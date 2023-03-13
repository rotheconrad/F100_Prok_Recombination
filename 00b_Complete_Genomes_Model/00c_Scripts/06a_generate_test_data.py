 #!/usr/bin/env python

'''Generates genomes around ani values based on gene RBMs along ANI gradient

** Needs python vs. 3.9 plus for rand.sample count parameter

    - uses a gamma distribution to model RBM gene distribution around ANI
    - Genome evolution is neutral random based on RBM distriubtion.
    - Average gene length and stdev (default 1000 bp, 250 stdev)
    - * no deletions, insertions or rearrangements * (future?)
    - Gene synteny is conserved
    - Genes start with ATG and stop with TAG
    - 10 bps between every gene (see future additions)
    - Genome starts and ends with 10 bps of space before/after genes
    - Pangenome Gene classes randomly distributed (Core or Accessory.
    - Creates draft genomes choosing n contigs from gaussian mu=10,sigma=5.
    - Contigs number stdev = average contigs / 2
    - creats 15% * source genes additional auxiliary genes.
    - shuffles auxillary genome choosing shared fraction from gaussion
      where mu=y=0.02(ANI)-1.05 and sigma=0.02

No input files are needed.

Writes a genome.fna, a cds.fnn for each genome.
Writes a gene cluster file with pangenome classifications.

The random.seed is set at 42 for reproducibility

# future additions
    - operons - sizes are randomly sampled from gaussian mu=10genes, std=5
    - intergenic regions are randomly sampled from gaussian mu=5bp, std=2
    - interoperon regions are randomly sampled from gaussian mu=100bp, std=25
    - recombination rate default 10%? simulate recombination event in 1/10
      of the genomes. How to keep ANI the same distance with increased F100?
      if some genes are more similar some genes need to be less similar.
      * also what fraction of genes are recombined?
    - accept input genome and gene fasta as seed genome
    - add mutation rate per generation per number of generations feature
        - no ANI model but keep gene and pangenome function
    - use daughter genomes as new source genomes for future generations

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: October 2022
License :: GNU GPLv3
Copyright 2022 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse, random, os
from pathlib import Path
import matplotlib.pyplot as plt


################################################################################
##### Generate source genome ###################################################
################################################################################

def get_pangenome_dist(gnum, crat):
    # calculates pangenome class numbers based on ratios and gene number

    arat = 1 - crat # auxilliary genes ratio
    hrat = 0.05 # highly conserved 5% of core genes
    crat = crat - hrat # new core ratio after highly conserved

    hnum = int(gnum*hrat) # number of highly conserved genes
    cnum = int(gnum*crat) # number of core genes
    anum = int(gnum*arat) # number of auxiliary genes

    pandist = [hnum, cnum, anum] # pangenome distribution array
    # adjust for rounding errors to equal total number of genes
    pandist_total = sum(pandist)
    if pandist_total < gnum:
        pandist[2] += gnum - pandist_total

    return pandist


def generate_source(gnum, pandist, tgenomes, mu_gene, stdev_gene):
    # creates the source genome and genes with pangenome classes

    ends = [] # 10bp ATGC string for genome beginning and end
    inters = [] # 10 bp ATGC string for between each gene
    genes = [] # list of lists: [gname, start, stop, pancat, sequence]
    axgenes = [] # additional auxillary genes to include with some genomes
    #glen = [20] # gene + intergene length. starts at 20 for 10 bp end caps
    # axgenes not included in source genome. These are extra auxiliary genes
    # to scatter into additional genomes.

    # random sample weighted pangenome classes
    classes = ['Conserved', 'Core', 'Accessory']
    pancats = random.sample(classes, k = gnum, counts = pandist)

    # generate 10 bp end caps
    for i in range(2):
        cap = ''.join(random.choice('CGTA') for _ in range(10))
        ends.append(cap)

    # generate gene and intergene sequences
    running_length = 10
    for i in range(gnum):
        gname = f'{i+1:03}' # gene name
        # generate gene length
        # start sequence ATG stop sequence TAG
        # random sample gene length avg gene length 1000 w/ stdev 250
        avg_gene_len = mu_gene - 6 # remove 6 bp for ATG and TAG 
        gl = int(random.gauss(avg_gene_len, stdev_gene))
        # implement minimum gene length 300 bp - ATG - TAG
        while gl < 294:
            gl = int(random.gauss(avg_gene_len, stdev_gene))
        strt = running_length + 1 # start position
        stop = running_length + gl + 6 # stop position
        pancat = pancats[i] # grab pancat
        cluster = f"Cluster_{i+1:03}"
        # generate random gene sequence
        rnd = ''.join(random.choice('CGTA') for _ in range(gl))
        gnseq = "ATG" + rnd + "TAG" # gene sequence
        genes.append([gname, strt, stop, pancat, gnseq, cluster])
        # create intergene sequence
        inters.append(''.join(random.choice('CGTA') for _ in range(10)))
        running_length += len(gnseq)+10

    # generate additional auxillary genes to sample from to simulate the shared
    # genome fraction.
    # get the number of accessory genes.
    accs = [i for i in pancats if i == 'Accessory']
    for i in range(len(accs)):
        gl = int(random.gauss(994, 250))
        # generate random gene sequence
        rnd = ''.join(random.choice('CGTA') for _ in range(gl))
        gnseq = "ATG" + rnd + "TAG" # gene sequence
        cluster = f"axCluster_{i:03}"
        axgenes.append([gnseq, cluster])

    # compile source genome dict to pass on
    source_genome = { 
                "genes": genes, "inters": inters,
                "ends": ends, "axgenes": axgenes,
                "glen": running_length #sum(glen)
                }

    return source_genome

################################################################################
##### Write output files #######################################################
################################################################################

def find_contig(cut_sites, strt):
    # finds the contig a gene is on based on the start site
    # returns the contig number and the length to contig start
    for i, v in enumerate(cut_sites):
        if strt >= cut_sites[-1]: return len(cut_sites) + 1, cut_sites[-1]
        elif strt < v: return i + 1, cut_sites[i-1]


def write_genome(gdict, genome_name, gcon, outpre, tPanCats):
    # writes genome.fna and cds.fnn fasta files from gdict

    # Build the full genome sequence string as seq_out
    seq_out = gdict["ends"][0] # add the first 10 bp start cap.

    ### divide genome sequence into random contig lengths
    # pick contigs number of random numbers in range len(seq_out)
    # I want the first ontig to be at least 5000bp
    # I want the last contig to be at least 5000bp
    # I want all contigs to be at least 5000
    # need the sum of gene lengths
    tmp_length = gdict["glen"] - 5000
    # randomly select the number of contigs
    contigs = int(random.gauss(gcon, gcon/4))
    while contigs < 2: contigs = int(random.gauss(gcon, gcon/4))
    # randomly select cut sites for contigs
    cuts = random.sample(range(5000, tmp_length, 5000), contigs)
    cut_sites = sorted(cuts)
    # write gene sequences to cds fasta and build genome sequence
    cdsout = f'{outpre}_FNNs/{outpre}_{genome_name}.fnn'
    with open(cdsout, 'w') as cout, open(tPanCats, 'a') as pout:
        for gene, inter in zip(gdict["genes"], gdict["inters"]):
            strt = gene[1]
            stop = gene[2]
            pancat = gene[3]
            gseq = gene[4]
            cluster = gene[5]
            # get contig number based on gene start position
            # and get the cutsite abrv: ctst
            try: contig, adjlen = find_contig(cut_sites, strt)
            except: print(f"find_contig error line 184", cut_sites, strt)
            if strt < cut_sites[0]: adjlen = 0
            # adjust strt and stop relative to the contig
            # adjlen is the position of the cut site
            # subtract this from start site to make it relative to contig
            strt = strt - adjlen
            stop = stop - adjlen
            # construct the gene name
            try: gname = f"{genome_name}_{contig:03}_{gene[0]}"
            except: print(f"gname error line 165 {genome_name} {contig} {gene[0]}")
            # write cds fasta sequence out in Prodigal style format
            cout.write(f">{gname} # {strt} # {stop} # {len(gseq)}\n{gseq}\n")
            # write pancat out
            pout.write(f"{gname}\t{cluster}\t{pancat}\tn/N\n")
            seq_out += gseq # add gene sequence to the full genome sequence
            seq_out += inter # add intergene regoin to the full genome sequence

    # add the last 10 bp end cap
    seq_out += gdict["ends"][1]

    # slice genome (seq_out) into contigs and write fasta
    with open(f'{outpre}_FNAs/{outpre}_{genome_name}.fna', 'w') as gout:
        j = 0 # previoius position for adjustment
        for i, cut in enumerate(cut_sites):
            pos = cut - j
            conout = seq_out[:pos]
            seq_out = seq_out[pos:]
            lineout = f">{genome_name}_{i+1:03}\n{conout}\n"
            gout.write(lineout)
            j = cut # set previous position.
        # write out last contig
        lineout = f">{genome_name}_contig_{len(cut_sites)+1:03}\n{seq_out}\n"
        gout.write(lineout)

    return True


def plot_hist(x, xmin, xmax, xstep, outpre, gname):

    m = sum(x) / len(x)

    fig, ax = plt.subplots(figsize=(7, 5))

    ax.hist(x, bins=100, edgecolor='k', color='#08519c', alpha=0.6,)
    ax.axvline(m, color='#ae017e', linestyle='dashed', linewidth=1)
    ax.set_xlabel('Average sequence distance of cluster', fontsize=14)
    ax.set_ylabel('Count', fontsize=14)

    # make adjustments for any floats in xmin, xmax, or xstep with range()
    hstep = xstep/10
    imin, imax, istep = int(xmin*100), int(xmax*100), int(xstep*100)
    xticks = [i/100 for i in range(imin,imax+istep,istep)]
    ax.set_xticks(xticks)
    ax.set_xlim(xmin-hstep, xmax+hstep)
    ax.tick_params(
     which='minor', axis='both', left=False, bottom=True
     )
    ax.tick_params(
             which='major', axis='both',
             left=True, bottom=True,
             size=6, width=2, tickdir='inout',
             labelsize=12, zorder=10
             )
    ax.yaxis.grid(
     which="major", color='#bdbdbd', linestyle='--',
     linewidth=1, alpha=0.4, zorder=1
     )
    ax.set_axisbelow(True)
    for spine in ax.spines.values(): spine.set_linewidth(2)

    fig.set_tight_layout(True)
    plt.savefig(f'{outpre}_PDFs/{gname}_rbms.pdf')
    plt.close()

################################################################################
##### Generate draft genomes ###################################################
################################################################################

def generate_genome(gnum, SG, ani, outpre, gname):
    # generates genome variants from source genome
    # uses gene list from SG["genes"]
    # which contains: [gene, start, stop, pancat, seq, cluster]
    # start by generating the gene reciprocal best matches around the ANI
    rbms = [get_gamma(ani) for i in range(gnum)]
    _ = plot_hist(rbms, 70, 100, 5, outpre, gname)
    # set variables for ordered lists for genes
    gene_names, pancats, seqs, clusters = [], [], [], []
    for g in SG["genes"]:
        gene_names.append(g[0])
        pancats.append(g[3])
        seqs.append(g[4])
        clusters.append(g[5])
    # initialize lists for evolved sequences
    new_ends = [] # 10bp sequence for genome beginning and end
    new_inters = [] # 10 bp ATGC string for between each gene
    new_genes = [] # # list of lists: [gname, start, stop, pancat, sequence]
    # set variables for additional auxiliary genes
    axgeneseq, axgeneclust = [], []
    for g in SG["axgenes"]:
        axgeneseq.append(g[0])
        axgeneclust.append(g[1])
    # reorder rbms for conserved genes
    rbms = reorder_rbms_conserved(pancats, rbms)
    # randomly select the shared fraction based on ANI value
    # at 95% ANI we swap 100% of accessory genes
    mu =  (-15*ani) + 1525 # (0.03*ani) - 2
    shared_frac = random.gauss(mu=mu, sigma=5)
    # iterate through source genes and create new genes based on rbms
    # iterate through the ordered gene and other lists to create new genome
    # initialize variables to store new gene sequences
    running_length = 10
    data = zip(gene_names, pancats, seqs, clusters, rbms)
    for gname, pancat, seq, clust, rbm in data:
        # switch auxilliary genes shared_frac percent of the time
        # set tempory axgeneseq and axgeneclust to sample without replacement
        tmpgseq = axgeneseq
        tmpgclst = axgeneclust
        if pancat == 'Accessory':
            event = random.randint(1,10000)/100
            if event >= shared_frac:
                randv = random.randint(0, len(axgeneseq)-1)
                seq = tmpgseq.pop(randv)
                clust = tmpgclst.pop(randv)
        # evolve the gene sequence based on the rbm and build new genome
        new_gene = evolve_seq(seq, rbm, True)
        strt = running_length + 1 # start position
        stop = running_length + len(new_gene) # stop position

        new_genes.append([gname, strt, stop, pancat, new_gene, clust])
        # increase sequence counter
        running_length += len(new_gene) + 10
        
    # Evolve the intergenic regions
    # join the start, end and inter lists to one string
    ends = SG["ends"]
    inter_seq = ends[0] + ''.join(SG["inters"]) + ends[1]
    # evolve joined sequence
    evo_seq = evolve_seq(inter_seq, ani, False)
    # split evo_seq back up into bits
    # everything is 10 currently because intergene segments are static at 10.
    new_ends.append(evo_seq[:10])
    new_ends.append(evo_seq[-10:])
    for i in range(10, len(evo_seq)-10, 10):
        new_inters.append(evo_seq[i:i+10])



    # compile source genome dict to pass on
    draft_genome = { 
                "genes": new_genes, "inters": new_inters,
                "ends": new_ends, "glen": running_length #sum(glen)
                }

    return draft_genome


def get_gamma(ani):
    # samples from gamma distribution around ani as mean.
    # uses python's random package
    # shape parameter k = alpha = 1.5
    # scale parameter = beta = (100-ani)/2
    # that sets the mean ~ equal to ani
    a = 1.5 # alpha or k for gamma distribution
    adjANI = 100.15 # adjust upper bound slightly to fit data better
    g = round(adjANI - random.gammavariate(a, (adjANI-ani)/a), 1)
    return g

def reorder_rbms_conserved(pancats, rbms):
    # all highly conserved genes need to be F100 if possible
    # get array positions where rbms == 100
    # asign 10% of these as highly conserved genes
    conserved_label_pos = [i for i,x in enumerate(pancats) if x == 'Conserved']
    # make sure there are more conserved rbm values than labels
    # if not, lower the conserved value
    conval=100
    conserved_rbm_pos = [i for i,x in enumerate(rbms) if x >= conval]
    while len(conserved_label_pos) > len(conserved_rbm_pos):
        conval -= 0.1
        conserved_rbm_pos = [i for i,x in enumerate(rbms) if x >= conval]
    # switch around the rbm array positions to reflect the conserved gene posi
    for lab_pos in conserved_label_pos:
        corresponding_rbm = rbms[lab_pos]
        if corresponding_rbm < conval:
            # pick random conserved rbm value from conserved rbm list
            rando = random.choice(range(len(conserved_rbm_pos)))
            swap = conserved_rbm_pos.pop(rando)
            new_val = rbms[swap]
            old_val = corresponding_rbm
            rbms[lab_pos] = new_val
            rbms[swap] = old_val

    return rbms


def evolve_seq(seq, rbm, remove_ends=False):
    # modifies positions in seq to equal rbm difference
    # returns new evolved sequence
    if remove_ends:
        seq = seq[3:-3]
    # get sequence length
    slen = len(seq)
    # calculate number of base pair changes to equal rbm
    if rbm > 100: rbm = 100
    frac = 1 - rbm/100
    changes = int(slen * frac)
    # select random positions in the sequence
    try: rnd = random.sample(range(slen), k=changes)
    except: print(f'\n\nLine 352: {slen} {changes} {rbm}\n{seq}\n\n')
    # iterate over positions and change base call
    for i in rnd:
        nucs = ['A', 'T', 'C', 'G']
        # get the source base at the rnd position
        source_base = seq[i]
        # remove the source base from the random choices
        nucs.pop(nucs.index(source_base))
        # pick random base to replace it with
        new_base = random.choice(nucs)
        # substitute new base
        seq = seq[:i] + new_base + seq[i+1:]

    if remove_ends:
        seq = "ATG" + seq + "TAG"

    return seq

################################################################################
################################################################################

def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-o', '--output_file_prefix',
        help='Please specify a prefix for the output file names!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-n', '--number_of_genomes',
        help='(OPTIONAL) Specify number of genomes to generate (default 10)!',
        metavar='',
        type=int,
        required=False,
        default=10
        )
    parser.add_argument(
        '-a', '--ani_range',
        help='(OPTIONAL) Specify ANI range min max step (default: 95 100 0.1)!',
        metavar='',
        type=float,
        nargs=3,
        required=False,
        default=[95, 100, 0.1]
        )
    parser.add_argument(
        '-g', '--number_source_genes',
        help='(OPTIONAL) Specify the number of source genes (default 3000)!',
        metavar='',
        type=int,
        required=False,
        default=3000
        )
    parser.add_argument(
        '-c', '--average_number_contigs',
        help='(OPTIONAL) Specify the number of contigs (default 20)!',
        metavar='',
        type=int,
        required=False,
        default=10
        )
    parser.add_argument(
        '-cr', '--core_gene_ratio',
        help='(OPTIONAL) Specify the core gene ratio (default 0.85)!',
        metavar='',
        type=float,
        required=False,
        default=0.85
        )
    parser.add_argument(
        '-mu', '--avg_gene_length',
        help='(OPTIONAL) Specify average gene length (default 1000)!',
        metavar='',
        type=int,
        required=False,
        default=1000
        )
    parser.add_argument(
        '-sd', '--gene_standard_deviation',
        help='(OPTIONAL) Specify gene length standard deviation (default 250)!',
        metavar='',
        type=int,
        required=False,
        default=250
        )
    parser.add_argument(
        '-rs', '--random_seed',
        help='(OPTIONAL) Set random seed (default 42)!',
        metavar='',
        type=int,
        required=False,
        default=42
        )
    args=vars(parser.parse_args())

    # Do what you came here to do:
    print('\n\nRunning Script...\n')

    # define parameters
    outpre = args['output_file_prefix']
    n = args['number_of_genomes']
    ani = args['ani_range']
    gnum = args['number_source_genes']
    gcon = args['average_number_contigs']
    crat = args['core_gene_ratio'] # core gene ratio
    mu_gene = args['avg_gene_length'] # average gene length
    stdev_gene = args['gene_standard_deviation'] # gene standard deviation

    # set the random seed
    #random.seed(args['random_seed'])

    # createout put directories
    # get outpre path
    p = outpre.split('/')
    if len(p) > 1:
        path = '/'.join(p[:-1])
    else: path = os.getcwd()
    # create directory for RBM hist plots
    Path(f"{path}/{outpre}_PDFs").mkdir(parents=True, exist_ok=True)
    # create directory for genomes
    Path(f"{path}/{outpre}_FNAs").mkdir(parents=True, exist_ok=True)
    # create directory for genes
    Path(f"{path}/{outpre}_FNNs").mkdir(parents=True, exist_ok=True)

    # convert floats (ani percent) to ints
    anix = [int(i*100) for i in ani]
    # generate the range and convert back to floats (ani percent)
    # ani_range = [ i/100 for i in range(anix[0], anix[1]+anix[2], anix[2])]
    ani_range = [ i/100 for i in range(anix[0], anix[1], anix[2])]
    # number of genomes to generate
    tgenomes = len(ani_range)

    #####################################################################
    ##### Generate Source ###############################################
    #####################################################################

    # pangenome distribution ratio
    pandist = get_pangenome_dist(gnum, crat)
    # initialize RBM and PanCat for concatenated results
    tRBMs = f"{outpre}_total_rbms.tsv"
    tPanCats = f"{outpre}_total_pancats.tsv"
    # write tPanCats file header
    with open(tPanCats, 'w') as out:
        out.write("Gene_Name\tCluster_Name\tPangenome_Category\tn/N\n")

    # generate source genome SG is a dict of 
    # { "genes": genes, "inters": inters, "ends": ends,
    #   "axgenes": axgenes, "glen": sum(glen) }
    SG = generate_source(gnum, pandist, tgenomes, mu_gene, stdev_gene)
    _ = write_genome(SG, "0SourceGenome", gcon, outpre, tPanCats)

    #####################################################################
    ##### Generate Genomes ##############################################
    #####################################################################

    # iterate over n for each ani to create n genomes for each ani in range
    for a in ani_range:
        for i in range(n):
            print(f'\t\tGenerating genomes {a} ANI genome {i+1}')
            gname = f"ANI{int(a*100)}Genome{i+1:03}"
            Gs = generate_genome(gnum, SG, a, outpre, gname)
            _ = write_genome(Gs, gname, gcon, outpre, tPanCats)

    #####################################################################

    print(f'\n\nComplete success space cadet!! Finished without errors.\n\n')
    
if __name__ == "__main__":
    main()
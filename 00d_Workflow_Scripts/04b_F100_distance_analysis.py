 #!/usr/bin/env python

'''Tests distance between recombing vs non-recombing genes.

This scripts requires 6 input files:

1) A reciprocal best match (RBM) file between two genomes from the AAI.rb
script of the enveomics collection.

2) Pangenome categories tsv file from 04d_Get_Genes_Clusters_PanCat.py.

3-4) Prodigal CDS in fasta format for each genome.
     Inputs cA and cB in the same order as they were given to AAI.rb.

5-6) The corresponding genome fasta fasta files.
     Inputs gA and gB in the same order as they were given to AAI.rb

This script returns data, graphics and statistics concerning the genomic
positions between genes categorized by highly conserved core genes,
recombining core genes, recombing accessory genes, and non-recombining
genes. A gene is said to be recombining if the RBM value equals 100 and
if it is not highly conserved.

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: July 2022
License :: GNU GPLv3
Copyright 2022 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse
from collections import defaultdict
import pandas as pd; import numpy as np
pd.options.mode.chained_assignment = None
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy import stats
import seaborn as sns
import statsmodels.api as sm


###############################################################################
##### SECTION 01: PARSE THE INPUT FILES #######################################
###############################################################################

def read_fasta(fp):
    '''This generator function reads fasta input and returns name, seq'''
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))


def parse_genome_fasta(gA, gB, switch):
    # parses the genome fasta files and retrieves the length of each
    # contig/chromosome. Creates nested dictionaries for both genomes
    # containing dictionaries of {contigID: len(contig)} for each genome.

    genomes = {
            'A': defaultdict(list),
            'B': defaultdict(list)
            }

    for i, G in enumerate([gA, gB]):
        with open(G, 'r') as file:
            for name, seq in read_fasta(file):
                contigID = name[1:] # remove the ">" from the fasta defline
                seqlen = len(seq)
                genomes[switch[i]][contigID] = seqlen

    # get the unique genome ID used in the fasta file for each genome
    # used to match genome A or genome B in the RBM file.
    gAid = '_'.join(list(genomes['A'].keys())[0].split('_')[:-1])
    gBid = '_'.join(list(genomes['B'].keys())[0].split('_')[:-1])

    return genomes, gAid, gBid


def parse_prodigal_CDS(cA, cB, switch):
    # parses the prodigal fasta and retrieves gene start, stop, and strand
    # from the fasta deflines. Creates nested dictionaries for both genomes
    # containing dictionaries of {GeneID: [values]} for each genome.
    CDS = {
        'A': defaultdict(lambda: defaultdict(list)),
        'B': defaultdict(lambda: defaultdict(list))
        }

    for i, C in enumerate([cA, cB]):
        with open(C, 'r') as file:
            for name, seq in read_fasta(file):
                X = name.split(' # ')
                ID = X[0][1:].split('_') # gene ID, '>' removed from fasta defline.
                contigID = '_'.join(ID[:-1])
                gID = ID[-1]
                start = int(X[1]) # gene start
                stop = int(X[2]) # gene stop
                strand = int(X[3]) # gene coding strand

                CDS[switch[i]][contigID][gID] = [start, stop, strand]

                if start > stop:
                    print(
                        f'parse_prodigal_CDS ERROR: start position greater '
                        f'than stop positoin for: \n{name}'
                        )

    return CDS


def parse_aai_rbm(rbm, gAid, gBid):

    # Parses the rbm file and returns a dictionary with lists as values.
    # Dictionary with a list of gene IDs for genome A and genome B
    # F100 is genes with 100% sequence similarity are set = 1
    # else F100 is set to 0 indicates < 100% RBMs

    RBM = {
            'A': defaultdict(lambda: defaultdict(int)),
            'B': defaultdict(lambda: defaultdict(int))
            }

    with open(rbm, 'r') as file:
        for line in file:
            X = line.rstrip().split('\t')
            A = X[0].split('_')
            contigIDA = '_'.join(A[:-1])
            genomeIDA = '_'.join(A[:-2])
            geneIDA = A[-1]
            B = X[1].split('_')
            contigIDB = '_'.join(B[:-1])
            genomeIDB = '_'.join(B[:-2])
            geneIDB = B[-1]
            ID = float(X[2])

            F100 = 1 if ID == 100 else 0

            if genomeIDA == gAid and genomeIDB == gBid:
                RBM['A'][contigIDA][geneIDA] = F100
                RBM['B'][contigIDB][geneIDB] = F100

    return RBM


def parse_pangenome_categories(PC):

    # parses the pangenome category file to dictionary of {gene: pancat}
    # where pancat equals Core, Accessory, or Specific.
    # Returns the dictionary

    pancats = {}

    with open(PC, 'r') as file:
        header = file.readline()
        for line in file:
            X = line.rstrip().split('\t')
            gene = X[0]
            pancat = X[2]
            pancats[gene] = pancat

    return pancats

###############################################################################
##### SECTION 02: COMPUTATIONS AND PLOTS ######################################
###############################################################################

def combine_input_data(RBM, CDS, genomes, pancats):
    # combines the parsed input into a single dataframe.
    # concatenates contigs and adjusts gene start and stop positions
    # returns pandas dataframe with columns of dict D keys and
    # returns contig positions for each genome.

    D = {
        'Genome': [], 'Gene': [], 'F100': [],
        'PanCat': [], 'Start': [], 'Stop': [], 'Strand': [], 
        }

    contig_positions = {'A': [], 'B': []}

    for genome, contigs in genomes.items():
        genome_length = 0
        for contig, contig_length in contigs.items():
            contig_genes = CDS[genome][contig]
            for gene, geneInfo in contig_genes.items():
                # 100% RBM = 1 or else = 0
                F100 = RBM[genome][contig][gene]
                # occasionally (< 1/1000) some genes do not make the pancat file
                # this seems to be from the mmseqs align step some sequences are
                # dropped from clusters failing an evalue cutoff. Since we 
                # have the evalue set super high and the occasional gene is
                # still removed. I label these genes and genome specific as in
                # they are not belonging to a gene cluster. This mmseqs issue
                # is noted in Part 03, Step 05 of the github readme.
                pc = pancats.get(f'{contig}_{gene}', 'Specific')
                start = geneInfo[0] + genome_length
                stop = geneInfo[1] + genome_length
                strand = geneInfo[2]

                D['Genome'].append(genome)
                D['Gene'].append(f'{contig}_{gene}')
                D['F100'].append(F100)
                D['PanCat'].append(pc)
                D['Start'].append(start)
                D['Stop'].append(stop)
                D['Strand'].append(strand)

            # increment genome_length. used to adjust gene start positions
            # and keep track of contig positions to add markers to plot
            genome_length += contig_length
            contig_positions[genome].append(genome_length)

    # convert dictionary to dataframe
    df = pd.DataFrame(D) #.sort_values(by=['Genome', 'Start'])
    # add gene lengths to Width column
    df['Width'] = df['Stop'] - df['Start']

    return df, contig_positions


def build_some_plots(df, genomes, pancats, outpre, cpos):

    # cpos is contig position array to mark them on the plot.
    
    colors = [
            '#fee391', # yellow for highly conserved genes (HC)
            '#dd3497', # pink recombining core (RC) 
            '#41ab5d', # green recombining accessory (RA)
            '#6baed6', # light blue non-recombining (NR)
            '#f03b20', # red genome specific genes (GS)
            '#000000', # neutral gray non-coding genome
            ]

    genome_out = f'{outpre}_genomes'

    fig, ax = plt.subplots(figsize=(7, 6))

    ax.set_xlabel('Gene location in genome (bp)')

    # Define label order
    labord = [
            "gA-HC", "gA-RC", "gA-RA", "gA-NR", "gA-GS",
            "gB-HC", "gB-RC", "gB-RA", "gB-NR", "gB-GS"
            ]

    # For each label add pc label and blank barplot point.
    for lab in labord:
        genome = lab[1]
        length = sum(genomes[genome].values())
        # neutral gray for genome length
        ax.barh(lab, length, left=0, color=colors[5], height=0.75)
        #ax.barh(lab, 1, left=0, color='w', height=0.75, alpha=0)

    pc_switch = {'Core': 'RC', 'Conserved': 'HC', 'Accessory': 'RA', 'Specific': 'GS'}

    for genome in genomes.keys():

        print(f'Computing and building plots for genome {genome} ...')

        # plot genes on the genome - select current genome (A or B)
        dfX = df[(df['Genome'] == genome)]
        dfX.sort_values(by='Start', inplace=True, ignore_index=True)
        dfX.to_csv(f'{outpre}_{genome}_gene_table.tsv', sep='\t', index=False)
        
        ## compute gene distances, run poisson simulation, KS test and plot
        dist_title = f'Genome {genome}'
        dist_out = f'{outpre}_{genome}_distance'
        length = sum(genomes[genome].values())
        _ = distance_plots(dfX, length, colors, dist_title, dist_out)

        genes = dfX['Gene'].to_numpy() # gene number
        pcat = dfX['PanCat'].to_numpy() # pancat
        Sta = dfX['Start'].to_numpy() # start position
        Wid = dfX['Width'].to_numpy() # Length or width of gene
        F10 = dfX['F100'].to_numpy() # F100 status
        Str = dfX['Strand'].to_numpy() # Strand
        #Pan = dfX['PanClass'].to_numpy() # Future - Pangenome category

        label = f'g{genome}'
    
        for G, P, S, W, F, Z in zip(genes, pcat, Sta, Wid, F10, Str):

            # get pancat abbreviation
            pc = pc_switch[P]

            # prep genome spceific genes for plot
            if pc == 'GS': c = colors[4]
            # prep non-recombinant genes for plot
            elif F < 1:
                c = colors[3]
                pc = 'NR'
            # prep recombinant genes for plot
            elif F == 1:
                if pc == 'RC': c = colors[1]
                elif pc == 'HC': c = colors[0]
                elif pc == 'RA': c = colors[2]
                else:
                    print('!!Panick!! How is a specific gene recombinant?')

            else:
                print('!!Panick!! something wrong lines 253-259 ish!')

            ylab = f'{label}-{pc}'
            ax.barh(ylab, W, left=S, color=c, height=0.75)


    # add contig markers
    mark_bars = {"gA-HC": ["v", -0.5], "gB-GS": ["^", 9.5]}
    for lab, marker in mark_bars.items():
        genome = lab[1]
        mark_pos = cpos[genome][:-1]
        ypos = [marker[1]] * len(mark_pos)
        ax.plot(mark_pos, ypos, marker="|", linestyle="", color='#969696')


    # Plot aesthetics
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.invert_yaxis()

    # add contig position markers

    # Build the legend
    legend_labels = [
                    'Highly Conserved (HC)',
                    'Recombant Core (RC)',
                    'Recombant Accessory (RA)',
                    'Non Recombinant (NR)',
                    'Genome Specific (GS)',
                    ]

    legend_elements = []

    for i, x in enumerate(legend_labels):
        n = Line2D(
            [0], [0], color='w', label=x, marker='s',
            markersize=15, markerfacecolor=colors[i]
            )
        legend_elements.append(n)
    marker_legend = Line2D(
            [0], [0], color='#969696', label='Contig Marker', marker="|",
            markersize=15, linestyle='None', #markerfacecolor='#969696'
            )
    legend_elements.append(marker_legend)

    ax.legend(
        handles=legend_elements,
        fontsize=12,
        fancybox=True,
        framealpha=0.0,
        frameon=False,
        loc='lower center',
        bbox_to_anchor=(0, 0.95, 1, 0.2),
        ncol=2
        )

    fig.set_tight_layout(True)
    plt.tick_params(left=False)
    plt.savefig(f'{genome_out}.pdf')
    plt.close()

    return True


###############################################################################
##### SECTION 03: STATISTICS AND PLOTS ########################################
###############################################################################

def poisson_simulation(array, n):

    ## 1st idea
    ## null hypothesis if genes are evenly spaced events every x base pairs
    ## the distance between them should follow a poisson distribution
    ## sig p value rejects the hypothesis and the genes are not evenly spaced.

    ## 2nd idea
    ## we tried 1st idea with distance between gene start sites. But since
    ## gene lengths are variable this was not a good way to test.
    ## We need to test distance between start and stop sites, or the number
    ## of genes between recombination events.

    # current version we test the number of genes between recombination or non-
    # recombination events. Null hypothesis is that if the distance between
    # events is evenly distributed they will follow a poisson distribution
    # high k and low p indicate we reject null hypothesis and the emperical
    # data do not follow a poisson distribution and are not evenly spaced.
    # hence they are randomly spaced events across the genome.

    emperical_mean_dist = np.diff(array)

    mu =  np.mean(emperical_mean_dist) # genes per kilo base pair
    poisson_array = stats.poisson.rvs(mu=mu, size=10000)
    # decided on a sample size of 10,000 instead of array length or gene count

    # Kolmogorov-Smirnov test emperical data against poisson distribution
    # returns ks_statistic, p_value
    k, p = stats.kstest(array, poisson_array)

    # Unit test / sanity check use poisson array instead of input array
    # k is small and p is large when the distributions are likely the same
    # fail to reject null hypothesis that the distributions are different.
    # poisson_A = stats.poisson.rvs(mu=mu, size=10000)
    # poisson_B = stats.poisson.rvs(mu=mu, size=10000)
    # k, p = stats.kstest(poisson_A, poisson_B)
    # emperical_mean_dist = poisson_A

    return emperical_mean_dist, poisson_array, mu, k, p


def qqplots(qqdata, distance_title, distance_out):
    # builds qq plots of distance between genes(emp) vs poisson simulation (psn)

    axlabsz = 12 # set ax title font size
    # initialize the plot object
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(
                                                2, 2, figsize=(7, 5),
                                                )
    # set the title
    #tt = f'Q-Q Plots: Genes Between Events\n{distance_title}'
    #fig.suptitle(tt, fontsize=14)

    if "nc_genes" in qqdata:
        emp = qqdata["nc_genes"]
        qdist = stats.poisson(np.mean(emp))
        sm.qqplot(emp, dist=qdist, line="45", ax=ax1)
    else:
        _ = ax1.text(
                0.5, 0.5, 'No Data',
                fontsize=18, ha='center',transform=ax1.transAxes
                )
    ax1.set_title('Non-recombinant genes', fontsize=axlabsz)
    ax1.set_xlabel('')
    if "conserved" in qqdata:
        emp = qqdata["conserved"]
        qdist = stats.poisson(np.mean(emp))
        sm.qqplot(emp, dist=qdist, line="45", ax=ax2)
    else:
        _ = ax2.text(
                0.5, 0.5, 'No Data',
                fontsize=18, ha='center',transform=ax2.transAxes
                )
    ax2.set_title('Conserved genes', fontsize=axlabsz)
    ax2.set_xlabel('')
    ax2.set_ylabel('')

    if "core" in qqdata:
        emp = qqdata["core"]
        qdist = stats.poisson(np.mean(emp))
        sm.qqplot(emp, dist=qdist, line="45", ax=ax3)
    else:
        _ = ax3.text(
                0.5, 0.5, 'No Data',
                fontsize=18, ha='center',transform=ax3.transAxes
                )
    ax3.set_title('Recombinant core genes', fontsize=axlabsz)

    if "accessory" in qqdata:
        emp = qqdata["accessory"]
        qdist = stats.poisson(np.mean(emp))
        sm.qqplot(emp, dist=qdist, line="45", ax=ax4)
    else:
        _ = ax4.text(
                0.5, 0.5, 'No Data',
                fontsize=18, ha='center',transform=ax4.transAxes
                )
    ax4.set_title('Recombinant accessory genes', fontsize=axlabsz)
    ax4.set_ylabel('')

    fig.set_tight_layout(True)
    plt.savefig(f'{distance_out}-qq.pdf')
    plt.close()

    '''
    # Test example plot for what a "good" Q-Q plot looks like
    print(qqdata)
    mu = np.mean(qqdata["nc_genes"])
    poisson_array = stats.poisson.rvs(mu=mu, size=10000) # emp analog
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.set_title("Poisson vs Poisson Example Q-Q Plot", fontsize=14)
    qdist = stats.poisson(np.mean(poisson_array))
    sm.qqplot(poisson_array, dist=qdist, line="45", ax=ax)

    fig.set_tight_layout(True)
    plt.savefig(f'{distance_out}-qqex.pdf')
    plt.close()
    '''

    return True


def distance_plots(df, length, colors, distance_title, distance_out):

    # Calculate distance between genes of four classes
    # Non-recombinant, conserved, recombinant core, and recombinant accessory
    # simulate poisson distribution
    # run kolmogorov-smirnov test
    # plot results
    # changed from contig length in bp to simply the length of the df
    length = len(df)
    # the dataframe index is the gene position in the genome
    # distance between genes is distance between index and the preceeding index 
    # subset conserved genes and select index as array
    conserved = df[df['PanCat'] == 'Conserved'].index
    # subset non conserved genes
    non_conserved = df[df['PanCat'] != 'Conserved']
    # subset non conserved and non recombining genes and select index as array
    nc_genes = non_conserved[non_conserved['F100'] < 1].index
    # subset non conserved but recombining genes
    recomb_genes = non_conserved[non_conserved['F100'] == 1]
    # subset recombinant core genes and select index as array
    core = recomb_genes[recomb_genes['PanCat'] == 'Core'].index
    # subset recombinant accessory genes and select index as array
    accessory = recomb_genes[recomb_genes['PanCat'] == 'Accessory'].index

    tx, ty = 0.72, 0.82 # stats test text position
    pcol = 'r' # poisson model kde line plot color
    mlw = 2 # kdeplot line width
    axlabsz = 12 # ax label sizes
    ts = 8 # stat test text size
    bw, ct = 3, 0 # kdeplot bw_adjus t and cut params
    minimum_genes = 5 # need at least 5 genes to test and plot

    # added this bit after the fact because I wanted to look at Q-Q plots
    # dicitonary of key nc_genes, conserved, core, accessory with
    # values of emp (the emperical distance between events array).
    qqdata = {} 

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(
                                                2, 2, figsize=(7, 5),
                                                #sharex=True, sharey=True
                                                )

    #tt = f'Genes Between Events\n{distance_title}'
    #fig.suptitle(tt, fontsize=14)

    if len(nc_genes) > minimum_genes:
        # calcuate distances, run poisson simulation, and ks test.
        emp, psn, mu, k, p = poisson_simulation(nc_genes, length)
        qqdata["nc_genes"] = emp
        _ = sns.histplot(x=emp, color=colors[5], ax=ax1, stat='density')
        _ = ax1.axvline(mu, color=colors[5], linestyle='dashed', linewidth=mlw)
        _ = sns.kdeplot(x=psn, color=pcol, ax=ax1, cut=ct, bw_adjust=bw)
        line = f'k = {k:.4f}\np = {p:.4f}'
        _ = ax1.text(tx, ty, line, transform=ax1.transAxes, fontsize=ts)
    else:
        _ = ax1.text(
                0.5, 0.5, 'No Data',
                fontsize=18, ha='center',transform=ax1.transAxes
                )
    ax1.set_ylabel('Density', fontsize=axlabsz)
    ax1.set_xlabel('Non-recombinant genes', fontsize=axlabsz)

    if len(conserved) > minimum_genes:
        emp, psn, mu, k, p = poisson_simulation(conserved, length)
        qqdata["conserved"] = emp
        _ = sns.histplot(x=emp, color=colors[0], ax=ax2, stat='density')
        _ = ax2.axvline(mu, color=colors[0], linestyle='dashed', linewidth=mlw)
        _ = sns.kdeplot(x=psn, color=pcol, ax=ax2, cut=ct, bw_adjust=bw)
        line = f'k = {k:.4f}\np = {p:.4f}'
        _ = ax2.text(tx, ty, line, transform=ax2.transAxes, fontsize=ts)
    else:
        _ = ax2.text(
                0.5, 0.5, 'No Data',
                fontsize=18, ha='center',transform=ax2.transAxes
                )
    ax2.set_ylabel('Density', fontsize=axlabsz)
    ax2.set_xlabel('Highly conserved genes', fontsize=axlabsz)

    if len(core) > minimum_genes:
        emp, psn, mu, k, p = poisson_simulation(core, length)
        qqdata["core"] = emp
        _ = sns.histplot(x=emp, color=colors[1], ax=ax3, stat='density')
        _ = ax3.axvline(mu, color=colors[1], linestyle='dashed', linewidth=mlw)
        _ = sns.kdeplot(x=psn, color =pcol, ax=ax3, cut=ct, bw_adjust=bw)
        line = f'k = {k:.4f}\np = {p:.4f}'
        _ = ax3.text(tx, ty, line, transform=ax3.transAxes, fontsize=ts)
    else:
        _ = ax3.text(
                0.5, 0.5, 'No Data',
                fontsize=18, ha='center',transform=ax3.transAxes
                )
    ax3.set_ylabel('Density', fontsize=axlabsz)
    ax3.set_xlabel('Recombinant core genes', fontsize=axlabsz)

    if len(accessory) > minimum_genes:
        emp, psn, mu, k, p = poisson_simulation(accessory, length)
        qqdata["accessory"] = emp
        _ = sns.histplot(x=emp, color=colors[2], ax=ax4, stat='density')
        _ = ax4.axvline(mu, color=colors[2], linestyle='dashed', linewidth=mlw)
        _ = sns.kdeplot(x=psn, color=pcol, ax=ax4, cut=ct, bw_adjust=bw)
        line = f'k = {k:.4f}\np = {p:.4f}'
        _ = ax4.text(tx, ty, line, transform=ax4.transAxes, fontsize=ts)
    else:
        _ = ax4.text(
                0.5, 0.5, 'No Data',
                fontsize=18, ha='center',transform=ax4.transAxes
                )
    ax4.set_ylabel('Density', fontsize=axlabsz)
    ax4.set_xlabel('Recombinant accessory genes', fontsize=axlabsz)

    #ax1.set_yscale('log')
    fig.set_tight_layout(True)
    plt.savefig(f'{distance_out}.pdf')
    plt.close()

    # build qq plots
    _ = qqplots(qqdata, distance_title, distance_out)

    return True

###############################################################################
###############################################################################

def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-rbm', '--RBM_from_AAI',
        help='Please specify the RBM file from AAI.rb!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-PC', '--pangenome_categories',
        help='Please specify the tsv from 04d_Get_Genes_Clusters_PanCat.py',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-cA', '--input_CDS_A',
        help='Please specify the first prodigal CDS in fasta format!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-cB', '--input_CDS_B',
        help='Please specify the second prodigal CDS in fasta format!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-gA', '--input_genome_A',
        help='Please specify the first genome file in fasta format!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-gB', '--input_genome_B',
        help='Please specify the second genome file in fasta format!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-m', '--minimum_contig_length',
        help='OPTIONAL: specify minimum contig length (default 2000)',
        metavar='',
        type=int,
        required=False,
        default=2000
        )
    parser.add_argument(
        '-o', '--output_file_prefix',
        help='Please specify a prefix for the output file names!',
        metavar='',
        type=str,
        required=True
        )
    args=vars(parser.parse_args())

    # Do what you came here to do:
    print('\n\nRunning Script...\n')

    # define parameters
    rbm = args['RBM_from_AAI']
    cA = args['input_CDS_A']
    cB = args['input_CDS_B']
    gA = args['input_genome_A']
    gB = args['input_genome_B']
    PC = args['pangenome_categories']
    m = args['minimum_contig_length']
    outpre = args['output_file_prefix']

    # setup a switch
    switch = {0: 'A', 1: 'B'}

    ## SECTION 01: Parse input data
    # read in the data
    # Get genome contig lengths
    genomes, gAid, gBid = parse_genome_fasta(gA, gB, switch)
    # get CDS info for each genome
    CDS = parse_prodigal_CDS(cA, cB, switch) 
    # get values from the rbm file
    RBM = parse_aai_rbm(rbm, gAid, gBid)
    # get the pancat info
    pancats = parse_pangenome_categories(PC)

    

    ## SECTION 02: Computations and Plots
    # process the data and build some plots
    # cpos is contig positions to use to mark them on the plot
    # cpos = {'A': [position array], 'B': [position array]}
    df, cpos = combine_input_data(RBM, CDS, genomes, pancats)
    _ = build_some_plots(df, genomes, pancats, outpre, cpos)
    ## SECTION 03 RUNS FROM INSIDE build_some_plots


    print(f'\n\nComplete success space cadet!! Finished without errors.\n\n')
    
if __name__ == "__main__":
    main()


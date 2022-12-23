 #!/usr/bin/env python

'''Tests distance between recombing (F100) vs non-recombing genes in genomes.

This script requires an reciprocal best match (RBM) file between two genomes
from the AAI.rb script of the enveomics collection.

It also requires a Prodigal CDS in fasta format for each genome.
Input f1 and f2 in the same order as they were given to the AAI.rb script.

And lastly, this script requires the two genome fasta files as well.

It returns statistics about the gene distances.

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
import numpy as np; import pandas as pd
pd.options.mode.chained_assignment = None
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import seaborn as sns
from scipy import stats

def read_fasta(fp):
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
    # contig/chromosome (C). Returns dictionary of {contigID: len(contig)}
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

    return genomes


def parse_prodigal_CDS(cA, cB, switch):
    # parses the prodigal fasta and retrieves gene start, stop, and strand
    # from the fasta deflines. Returns dictionary of {GeneID: [values]}
    CDS = {
        'A': defaultdict(lambda: defaultdict(list)),
        'B': defaultdict(lambda: defaultdict(list))
        }

    for i, C in enumerate([cA, cB]):
        with open(C, 'r') as file:
            for name, seq in read_fasta(file):
                X = name.split(' # ')
                ID = X[0][1:].split('_') # gene ID, '>' removed from fasta defline.
                contigID = '_'.join(ID[:2])
                gID = ID[2]
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


def parse_aai_rbm(rbm):

    # Parses the rbm file and returns a dictionary withs lists as values.
    # Dictionary with a list of gene IDs for genome A and genome B
    # F100 is genes with 100% sequence similarity
    # Fno is genes without 100% sequence similarity
    RBM = {'A': defaultdict(list), 'B': defaultdict(list)}

    with open(rbm, 'r') as file:
        for line in file:
            X = line.rstrip().split('\t')
            A = X[0].split('_')
            contigIDA = '_'.join(A[:2])
            geneIDA = A[2]
            B = X[1].split('_')
            contigIDB = '_'.join(B[:2])
            geneIDB = B[2]
            ID = float(X[2])

            F100 = 1 if ID == 100 else 0

            RBM['A'][contigIDA].append([geneIDA, F100])
            RBM['B'][contigIDB].append([geneIDB, F100])

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


def compute_gene_distances(RBM, CDS, pancats):
    # puts all the data together, runs statistics, and returns a dataframe.

    D = {
        'Genome': [], 'Contig': [], 'Gene': [], 'F100': [],
        'PanCat': [], 'Start': [], 'Stop': [], 'Strand': [], 
        }

    for genome, contigs in RBM.items():
        for contig, genes in contigs.items():
            for gene in genes:
                geneID = gene[0]
                F100 = gene[1]
                try: pc = pancats[f'{contig}_{geneID}']
                except: pc = 'NA'
                geneInfo = CDS[genome][contig][geneID]
                start = geneInfo[0]
                stop = geneInfo[1]
                strand = geneInfo[2]

                D['Genome'].append(genome)
                D['Contig'].append(contig)
                D['Gene'].append(geneID)
                D['F100'].append(F100)
                D['PanCat'].append(pc)
                D['Start'].append(start)
                D['Stop'].append(stop)
                D['Strand'].append(strand)

    df = pd.DataFrame(D).sort_values(by=['Genome', 'Contig', 'Start'])
    df['Width'] = df['Stop'] - df['Start']

    return df


def poisson_simulation(array, length, k=10000):

    ## The input array is the start position of genes on the contig/genome.
    ## The length is the length of contig/genome is base pairs.
    ## null hypothesis if genes are evenly spaced events every x base pairs
    ## the distance between them should follow a poisson distribution
    ## sig p value rejects the hypothesis and the genes are not evenly spaced.

    n = len(array)
    kbp_array = array/1000
    emperical_mean_dist = np.diff(kbp_array)
    mu =  np.mean(emperical_mean_dist) # genes per kilo base pair
    poisson_array = stats.poisson.rvs(mu=mu, size=n)

    # Kolmogorov-Smirnov test emperical data against poisson distribution
    # returns ks_statistic, p_value
    cdf = stats.poisson.cdf(n, mu)
    k, p = stats.kstest(kbp_array, 'poisson', args=(n, mu))

    ''' TO DO: build distribution of poisson samples to plot c.i.
    dist_array = defaultdict(list)

    for i in range(k):
        x = stats.poisson.rvs(mu=mu, size=n)
        dist = np.mean(np.diff(x))
        dist_array[i].append(dist)
    '''

    return emperical_mean_dist, poisson_array, mu, k, p


def distance_plots(df, length, colors, distance_title, distance_out):

    # simulate data for all genes
    all_genes = df['Start']
    conserved = df[df['PanCat'] == 'Conserved']['Start']
    core = df[df['PanCat'] == 'Core']['Start']
    accessory = df[df['PanCat'] == 'Accessory']['Start']

    tx, ty = 0.72, 0.82
    pcol = 'r'
    mlw = 2
    axlabsz = 12
    ts = 8
    bw, ct = 2, 3
    minimum_genes = 5

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(
                                                2, 2, figsize=(7, 5),
                                                #sharex=True, sharey=True
                                                )

    tt = f'Probability Density and KS test for gene distance\n{distance_title}'
    fig.suptitle(tt, fontsize=14)

    if len(all_genes) > minimum_genes:
        emp, psn, mu, k, p = poisson_simulation(all_genes, length)
        _ = sns.kdeplot(x=emp, color=colors[4], ax=ax1, cut=ct, bw_adjust=bw)
        _ = sns.kdeplot(x=psn, color=pcol, ax=ax1, cut=ct, bw_adjust=bw)
        _ = ax1.axvline(mu, color=colors[4], linestyle='dashed', linewidth=mlw)
        line = f'k = {k:.4f}\np = {p:.4f}'
        _ = ax1.text(tx, ty, line, transform=ax1.transAxes, fontsize=ts)
    else:
        _ = ax1.text(
                0.5, 0.5, 'No Data',
                fontsize=18, ha='center',transform=ax1.transAxes
                )
    ax1.set_ylabel('Density', fontsize=axlabsz)
    ax1.set_xlabel('All genes (kbp)', fontsize=axlabsz)

    if len(conserved) > minimum_genes:
        emp, psn, mu, k, p = poisson_simulation(conserved, length)
        _ = sns.kdeplot(x=emp, color=colors[0], ax=ax2, cut=ct, bw_adjust=bw)
        _ = sns.kdeplot(x=psn, color=pcol, ax=ax2, cut=ct, bw_adjust=bw)
        _ = ax2.axvline(mu, color=colors[0], linestyle='dashed', linewidth=mlw)
        line = f'k = {k:.4f}\np = {p:.4f}'
        _ = ax2.text(tx, ty, line, transform=ax2.transAxes, fontsize=ts)
    else:
        _ = ax2.text(
                0.5, 0.5, 'No Data',
                fontsize=18, ha='center',transform=ax2.transAxes
                )
    ax2.set_ylabel('Density', fontsize=axlabsz)
    ax2.set_xlabel('Highly conserved genes (kbp)', fontsize=axlabsz)

    if len(core) > minimum_genes:
        emp, psn, mu, k, p = poisson_simulation(core, length)
        _ = sns.kdeplot(x=emp, color=colors[1], ax=ax3, cut=ct, bw_adjust=bw)
        _ = sns.kdeplot(x=psn, color =pcol, ax=ax3, cut=ct, bw_adjust=bw)
        _ = ax3.axvline(mu, color=colors[1], linestyle='dashed', linewidth=mlw)
        line = f'k = {k:.4f}\np = {p:.4f}'
        _ = ax3.text(tx, ty, line, transform=ax3.transAxes, fontsize=ts)
    else:
        _ = ax3.text(
                0.5, 0.5, 'No Data',
                fontsize=18, ha='center',transform=ax3.transAxes
                )
    ax3.set_ylabel('Density', fontsize=axlabsz)
    ax3.set_xlabel('Recombinant core genes (kbp)', fontsize=axlabsz)

    if len(accessory) > minimum_genes:
        emp, psn, mu, k, p = poisson_simulation(accessory, length)
        _ = sns.kdeplot(x=emp, color=colors[2], ax=ax4, cut=ct, bw_adjust=bw)
        _ = sns.kdeplot(x=psn, color=pcol, ax=ax4, cut=ct, bw_adjust=bw)
        _ = ax4.axvline(mu, color=colors[2], linestyle='dashed', linewidth=mlw)
        line = f'k = {k:.4f}\np = {p:.4f}'
        _ = ax4.text(tx, ty, line, transform=ax4.transAxes, fontsize=ts)
    else:
        _ = ax4.text(
                0.5, 0.5, 'No Data',
                fontsize=18, ha='center',transform=ax4.transAxes
                )
    ax4.set_ylabel('Density', fontsize=axlabsz)
    ax4.set_xlabel('Recombinant accessory genes (kbp)', fontsize=axlabsz)

    #ax1.set_yscale('log')
    fig.set_tight_layout(True)
    plt.savefig(distance_out)
    plt.close()

    return True


def build_some_plots(df, genomes, pancats, outpre, mcl):
    
    colors = [
            '#54278f', # purple for conserved genes
            '#c51b7d', # bright pink recombining core
            '#4d9221', # dark green recombining accessory
            '#969696', # dark gray non-recombining
            '#bdbdbd', # neutral gray non-coding genome
            ]

    genome_out = f'{outpre}_genomes.pdf'

    fig, ax = plt.subplots(figsize=(7, 5))

    ax.set_xlabel('Gene location on contig (bp)')

    # set order of y-axis labels
    ylabel_set = []
    # get genome A and B contigs
    for genome, contigs in genomes.items():
        for i, (contig, length) in enumerate(contigs.items()):
            if length > mcl:
                ylabel_set.append(f'C{i+1:02}-G{genome}')

    # Sort label order
    ylabel_set.sort()
    # For each label add pc label and blank barplot point.
    for pc in ['HC', 'RC', 'RA', 'NR']:
        for L in ylabel_set:
            first_lab = f'{L}-{pc}'
            ax.barh(first_lab, 1, left=0, color='w', height=0.75, alpha=0)

    pc_switch = {'Core': 'RC', 'Conserved': 'HC', 'Accessory': 'RA'}

    for genome, contigs in genomes.items():

        print(f'Computing and building plots for genome {genome} ...')

        # plot genes on contigs
        for i, (contig, length) in enumerate(contigs.items()):

            # don't plot short contigs
            if length < mcl: continue

            dfX = df[(df['Genome'] == genome) & (df['Contig'] == contig)]

            ## statics gene distance plot
            dist_title = f'Genome {genome} Contig {contig}'
            dist_out = f'{outpre}_{genome}_{contig}_distance.pdf'
            _ = distance_plots(dfX, length, colors, dist_title, dist_out)
            
            genes = dfX['Gene'].to_numpy() # gene number
            Sta = dfX['Start'].to_numpy() # start position
            Wid = dfX['Width'].to_numpy() # Length or width of gene
            F10 = dfX['F100'].to_numpy() # F100 status
            Str = dfX['Strand'].to_numpy() # Strand
            #Pan = dfX['PanClass'].to_numpy() # Future - Pangenome category

            label = f'C{i+1:02}-G{genome}'
        
            for G, S, W, F, Z in zip(genes, Sta, Wid, F10, Str):
                if F < 1:
                    c = colors[3]
                    pc = 'NR'
                elif F == 1:
                    try:
                        pc = pc_switch[pancats[f'{contig}_{G}']]
                        if pc == 'RC': c = colors[1]
                        elif pc == 'HC': c = colors[0]
                        elif pc == 'RA': c = colors[2]
                        else:
                            print('!!Panick!! How is a specific gene recombinant?')
                    except:
                        c = colors[3]
                        pc = 'NR'
                        #print(f'{contig}_{G}')
                else:
                    print('!!Panick!! something wrong lines 253-259 ish!')

                #l = pos if Z == 1 else neg#removed pos neg strand labels
                #ax.barh(l, W, left=S, color=c, height=0.5)
                ylab = f'{label}-{pc}'
                ax.barh(ylab, W, left=S, color=c, height=0.75)

    # Plot aesthetics
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.invert_yaxis()

    # Build the legend
    legend_labels = [
                    'Highly Conserved (HC)',
                    'Recombant Core (RC)',
                    'Recombant Accessory (RA)',
                    'Non Recombinant (NR)',
                    #'Non-CDS'
                    ]

    legend_elements = []

    for i, x in enumerate(legend_labels):
        n = Line2D(
            [0], [0], color='w', label=x, marker='s',
            markersize=15, markerfacecolor=colors[i]
            )
        legend_elements.append(n)

    ax.legend(
        handles=legend_elements,
        fontsize=12,
        fancybox=True,
        framealpha=0.0,
        frameon=False,
        loc='lower center',
        bbox_to_anchor=(0, 1.02, 1, 0.2),
        ncol=2
        )

    fig.set_tight_layout(True)
    plt.tick_params(left=False)
    plt.savefig(genome_out)
    plt.close()

    return True


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
        '-PC', '--pangenome_categories',
        help='Please specify the tsv from 04d_Get_Genes_Clusters_PanCat.py',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-o', '--output_file_prefix',
        help='Please specify a prefix for the output file names!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-mcl', '--minimum_contig_length',
        help='(OPTIONAL) Specify the minimum contig length (default=500000)!',
        metavar='',
        type=int,
        required=False,
        default=500000
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
    outpre = args['output_file_prefix']
    mcl = args['minimum_contig_length']

    # setup a switch
    switch = {0: 'A', 1: 'B'}

    # read in the data
    RBM = parse_aai_rbm(rbm) # values from the rbm file
    CDS = parse_prodigal_CDS(cA, cB, switch) # get CDS info for each genome
    genomes = parse_genome_fasta(gA, gB, switch) # Get genome contig lengths
    pancats = parse_pangenome_categories(PC)

    # process the data
    df = compute_gene_distances(RBM, CDS, pancats)

    # build some plots
    _ = build_some_plots(df, genomes, pancats, outpre, mcl)


    print(f'\n\nComplete success space cadet!! Finished without errors.\n\n')
    
if __name__ == "__main__":
    main()


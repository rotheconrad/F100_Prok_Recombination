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


def compute_gene_distances(RBM, CDS):
    # puts all the data together, runs statistics, and returns a dataframe.

    D = {
        'Genome': [], 'Contig': [], 'Gene': [], 'F100': [],
        'Start': [], 'Stop': [], 'Strand': [], 
        }

    for genome, contigs in RBM.items():
        for contig, genes in contigs.items():
            for gene in genes:
                geneID = gene[0]
                F100 = gene[1]
                geneInfo = CDS[genome][contig][geneID]
                start = geneInfo[0]
                stop = geneInfo[1]
                strand = geneInfo[2]

                D['Genome'].append(genome)
                D['Contig'].append(contig)
                D['Gene'].append(geneID)
                D['F100'].append(F100)
                D['Start'].append(start)
                D['Stop'].append(stop)
                D['Strand'].append(strand)

    df = pd.DataFrame(D).sort_values(by=['Genome', 'Contig', 'Start'])
    df['Width'] = df['Stop'] - df['Start']

    return df


def distance_plots(df, colors, distance_title, distance_out):

    F100 = df[df['F100'] == 1]
    F100_dist = F100['Start'].diff()
    Fno = df[df['F100'] == 0]
    Fno_dist = Fno['Start'].diff()

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7, 5))

    fig.suptitle(distance_title, fontsize=20)

    _ = ax1.hist(F100_dist, rwidth=0.9, color='#c51b7d', alpha=0.6,)
    ax1.set_ylabel('Gene count', fontsize=14)
    ax1.set_xlabel('Distance between recombinging genes (bp)', fontsize=14)

    _ = ax2.hist(Fno_dist, rwidth=0.9, color='#4d9221', alpha=0.6)
    ax2.set_ylabel('Gene count', fontsize=14)
    ax2.set_xlabel('Distance between non-recombinging genes (bp)', fontsize=14)

    for ax in [ax1, ax2]:
        ax.set_yscale('log')
        ax.minorticks_on()
        ax.tick_params(
            which='minor', axis='both', left=False, bottom=False
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
    plt.savefig(distance_out)
    plt.close()

    return True


def build_some_plots(df, genomes, pancats, outpre):
    
    colors = [
            '#c51b7d', # bright pink recombining core
            '#4d9221', # dark green recombining accessory
            '#969696', # dark gray non-recombining
            '#bdbdbd', # neutral gray non-coding genome
            ]

    genome_out = f'{outpre}_genomes.pdf'

    fig, ax = plt.subplots(figsize=(7, 5))

    ax.set_xlabel('Gene location on contig (bp)')

    for genome, contigs in genomes.items():
        for i, (contig, length) in enumerate(contigs.items()):
            print(genome, contig, length)
            dfX = df[(df['Genome'] == genome) & (df['Contig'] == contig)]
            #print(dfX)
            distance_title = f'Genome {genome} Contig {contig}'
            distance_out = f'{outpre}_{genome}_{contig}_distance.pdf'
            _ = distance_plots(dfX, colors, distance_title, distance_out)

            genes = dfX['Gene'].to_numpy() # gene number
            Sta = dfX['Start'].to_numpy() # start position
            Wid = dfX['Width'].to_numpy() # Length or width of gene
            F10 = dfX['F100'].to_numpy() # F100 status
            Str = dfX['Strand'].to_numpy() # Strand
            #Pan = dfX['PanClass'].to_numpy() # Future - Pangenome category

            label = f'{genome}-{i+1:03}'
            ax.barh(label, length, left=1, height=0.5, color=colors[3])

            ''' removed positive and negative strand labels
            pos = f'{label} (+)'
            neg = f'{label} (-)'
            ax.barh(pos, length, left=1, height=0.5, color=colors[3])
            ax.barh(neg, length, left=1, height=0.5, color=colors[3])
            '''

            for G, S, W, F, Z in zip(genes, Sta, Wid, F10, Str):
                if F < 1:
                    c = colors[2]
                elif F == 1:
                    pc = pancats[f'{contig}_{G}']
                    if pc == 'Core': c = colors[0]
                    elif pc == 'Accessory': c = colors[1]
                    else:
                        print('!!Panick!! How is a specific gene recombinant?')
                else:
                    print('!!Panick!! something wrong lines 253-259 ish!')

                #l = pos if Z == 1 else neg#removed pos neg strand labels
                #ax.barh(l, W, left=S, color=c, height=0.5)
                ax.barh(label, W, left=S, color=c, height=0.5)

    # Plot aesthetics
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.invert_yaxis()
    ax.minorticks_on()
    ax.tick_params(
        which='minor', axis='both', left=False, bottom=False
        )
    ax.tick_params(
                which='major', axis='both',
                left=False, bottom=True,
                size=6, width=2, tickdir='inout',
                labelsize=12, zorder=10
                )
    ax.xaxis.grid(
        which="major", color='#bdbdbd', linestyle='--',
        linewidth=1, alpha=0.4, zorder=1
        )
    ax.set_axisbelow(True)
    for spine in ax.spines.values(): spine.set_linewidth(2)

    # Build the legend
    legend_labels = [
                    'Recombing core',
                    'Recombing accessory',
                    'Non-recombining',
                    'Non-CDS'
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

    # setup a switch
    switch = {0: 'A', 1: 'B'}

    # read in the data
    RBM = parse_aai_rbm(rbm) # values from the rbm file
    CDS = parse_prodigal_CDS(cA, cB, switch) # get CDS info for each genome
    genomes = parse_genome_fasta(gA, gB, switch) # Get genome contig lengths
    pancats = parse_pangenome_categories(PC)

    # process the data
    df = compute_gene_distances(RBM, CDS)

    # build some plots
    _ = build_some_plots(df, genomes, pancats, outpre)


    print(f'\n\nComplete success space cadet!! Finished without errors.\n\n')
    
if __name__ == "__main__":
    main()


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
        pclass = lab.split('-')[1]
        glength = sum(genomes[genome].values())
        dfY = df[(df['Genome'] == genome)]
        # black for full genome length background
        ax.barh(lab, glength, left=0, color=colors[5], height=0.75)
        # Add percent of total genes to end of each bar
        # this is probably an ugly way to do this that i've added at the end.
        total_genes = len(dfY)
        # iterate over each bar and annotate the value 
        if pclass == 'HC':
            # subset conserved genes 
            conserved = dfY[dfY['PanCat'] == 'Conserved']
            tc = colors[0]
            tl = round(len(conserved) / total_genes * 100, 2)
        elif pclass == 'RC':
            # subset recombinant core genes
            core = dfY[(dfY['PanCat'] == 'Core') & (dfY['F100'] == 1)]
            tc = colors[1]
            tl = round(len(core) / total_genes * 100, 2)
        elif pclass == 'RA':
            # subset recombinant accessory genes
            accessory = dfY[(dfY['PanCat'] == 'Accessory') & (dfY['F100'] == 1)]
            tc = colors[2]
            tl = round(len(accessory) / total_genes * 100, 2)
        elif pclass == 'NR':
            # subset non conserved and non recombining genes
            nc_genes = dfY[(dfY['PanCat'] != 'Conserved') & (dfY['PanCat'] != 'Specific') & (dfY['F100'] != 1)]
            tc = colors[3]
            tl = round(len(nc_genes) / total_genes * 100, 2)
        elif pclass == 'GS':
            # subset genome specific genes
            specific = dfY[dfY['PanCat'] == 'Specific']
            tc = colors[4]
            tl = round(len(specific) / total_genes * 100, 2)
        # add it to the plot
        ax.text(glength+(glength/100), lab, f'{tl}%', color=tc, va='center')

    pc_switch = {'Core': 'RC', 'Conserved': 'HC', 'Accessory': 'RA', 'Specific': 'GS'}

    for genome in genomes.keys():

        print(f'Computing and building plots for genome {genome} ...')
        # genome length
        glength = sum(genomes[genome].values())
        # plot genes on the genome - select current genome (A or B)
        dfX = df[(df['Genome'] == genome)]
        dfX.sort_values(by='Start', inplace=True, ignore_index=True)
        dfX.to_csv(f'{outpre}_{genome}_gene_table.tsv', sep='\t', index=False)

        ## compute gene distances, run poisson simulation, KS test and plot
        dist_title = f'Genome {genome}'
        dist_out = f'{outpre}_{genome}_distance'
        # changed doing by pangenome category to doing 1 for each genome.
        # and collapsing consecutive F100 genes to 1 event
        # returns percent recently recombining genes
        _ = distance_plots_2(dfX, colors, dist_title, dist_out)

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

    # the input array is the sorted index position of the genes in the genome
    # for each pangenome category. So the distance between index positions for
    # recombinant core genes is the number of genes between recombinant events
    # Using the diff method below the minimum distance is 1 instead of zero
    # because index position 4 - index position 5 equals 1.

    # make sure the array is sorted in ascending order.
    sorted_array = np.sort(array)

    # get the difference between gene index positions (genes between events)
    # and subtract 1 so the minimum distance is 0 instead of 1.
    temp_dist = np.diff(sorted_array) - 1

    # now we drop all zeros because theoretically adjacent genes both with
    # 100% sequence identity to their RBMs are a single recombination event.
    #emperical_mean_dist = temp_dist
    emperical_mean_dist = temp_dist[temp_dist != 0]
    #emperical_mean_dist = temp_dist[np.logical_and(temp_dist != 0, temp_dist != 1)]

    # print out value counts of distances between events
    values, counts = np.unique(emperical_mean_dist,return_counts=True)
    for v, c in zip(values, counts):
        print(v, c)
    

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


def distance_plots_2(df, colors, distance_title, distance_out):

    # Calculate distance between F100 gene events in the genome
    # Collapse consecutive F100 genes to single event.
    # simulate poisson distribution
    # run kolmogorov-smirnov test
    # plot results
    # changed from contig length in bp to simply the length of the df
    total_genes = len(df) # total events for poisson simulation
    # the dataframe index is the gene position in the genome
    # distance between genes is distance between index and the preceeding index 
    # subset non conserved genes ie disregard conserved genes
    # subset recombinant all genes that are not "conserved"
    recomb_genes = df[(df['F100'] == 1) & (df['PanCat'] != 'Conserved')].index

    # percent of genes that are recently recombining
    # 100% RBMs / total genes

    # skip statistical tests if fewer than 5 genes.
    if len(recomb_genes) < 5:
        print('Fewer than 5 recombinant genes. Skipping statistical test.')
        return False

    # calcuate distances, run poisson simulation, and ks test.
    emp, psn, mu, k, p = poisson_simulation(recomb_genes, total_genes)

    tx, ty = 0.72, 0.82 # stats test text position
    pcol = 'r' # poisson model kde line plot color
    mlw = 2 # kdeplot line width
    axlabsz = 12 # ax label sizes
    ts = 8 # stat test text size
    bw, ct = 3, 0 # kdeplot bw_adjus t and cut params

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7, 7))

    # plot emp distribution, poissoin distribution and stats test results
    _ = sns.histplot(x=emp, color=colors[5], ax=ax1, stat='density')
    _ = ax1.axvline(mu, color=colors[5], linestyle='dashed', linewidth=mlw)
    _ = sns.kdeplot(x=psn, color=pcol, ax=ax1, cut=ct, bw_adjust=bw)
    line = f'k = {k:.4f}\np = {p:.4f}'
    _ = ax1.text(tx, ty, line, transform=ax1.transAxes, fontsize=ts)
    ax1.set_ylabel('Density', fontsize=axlabsz)
    ax1.set_xlabel('Genes between recombinant events', fontsize=axlabsz)

    # plot Q-Q plot emperical vs theoretical
    qdist = stats.poisson(np.mean(emp))
    sm.qqplot(emp, dist=qdist, line="45", ax=ax2)

    fig.set_tight_layout(True)
    plt.savefig(f'{distance_out}.pdf')
    plt.close()

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
        help='Please specify the tsv from 04d_Get_Genes_Clusters_PanCat.py!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-a', '--annotation_file',
        help='Please specify the representative gene annotation file!',
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


    ## SECTION 03: Annotations
    # match genes to their annotations
    # partition into recombinant genes vs non-recombinant
    # test distributions and functional category differences

    print(f'\n\nComplete success space cadet!! Finished without errors.\n\n')
    
if __name__ == "__main__":
    main()


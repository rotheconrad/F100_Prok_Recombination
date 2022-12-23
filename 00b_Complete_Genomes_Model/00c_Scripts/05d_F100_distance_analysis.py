 #!/usr/bin/env python

'''Tests distance between recombing (F100) vs non-recombing genes.

This scripts requires 6 input files:

1) A reciprocal best match (RBM) file between two genomes from the AAI.rb
script of the enveomics collection.

2) Pangenome categories tsv file from 04d_Get_Genes_Clusters_PanCat.py.

3-4) Prodigal CDS in fasta format for each genome.
     Inputs cA and cB in the same order as they were given to AAI.rb.

5-6) The corresponding genome fasta fasta files.
     Inputs gA and gB in the same order as they were given to AAI.rb

This script returns graphics and statistics concerning the genomic
distance between genes categorized by highly conserved core genes,
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
import numpy as np; import pandas as pd
pd.options.mode.chained_assignment = None
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import seaborn as sns
from scipy import stats
import statsmodels.api as sm

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
        length = 0
        with open(G, 'r') as file:
            for name, seq in read_fasta(file):
                contigID = name[1:] # remove the ">" from the fasta defline
                seqlen = len(seq)
                length += seqlen
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
            contigIDA = '_'.join(A[:-1])
            geneIDA = A[-1]
            B = X[1].split('_')
            contigIDB = '_'.join(B[:-1])
            geneIDB = B[-1]
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
    # puts all the data together and returns a dataframe.

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

    df = pd.DataFrame(D)#.sort_values(by=['Genome', 'Contig', 'Start'])
    df['Width'] = df['Stop'] - df['Start']

    return df


def compute_gene_distances_draftmode(RBM, CDS, genomes, pancats):
    # puts all the data together and returns a dataframe.

    D = {
        'Genome': [], 'Contig': [], 'Gene': [], 'F100': [],
        'PanCat': [], 'Start': [], 'Stop': [], 'Strand': [], 
        }

    for genome, contigs in RBM.items():
        draft_length = 0
        for contig, genes in contigs.items():
            for gene in genes:
                geneID = gene[0]
                F100 = gene[1]
                try: pc = pancats[f'{contig}_{geneID}']
                except: pc = 'NA'
                geneInfo = CDS[genome][contig][geneID]
                start = geneInfo[0] + draft_length
                stop = geneInfo[1] + draft_length
                strand = geneInfo[2]

                D['Genome'].append(genome)
                D['Contig'].append(genome) # contig name == genome name draftmode
                D['Gene'].append(f'{contig}_{geneID}')
                D['F100'].append(F100)
                D['PanCat'].append(pc)
                D['Start'].append(start)
                D['Stop'].append(stop)
                D['Strand'].append(strand)

            # draftmode keep track of draft genome length of concatenated contigs
            contig_length = genomes[genome][contig]
            draft_length += contig_length
            #print(f'{genome}\t{contig}\t{contig_length}\t{draft_length}')

    df = pd.DataFrame(D).sort_values(by=['Genome', 'Contig', 'Start'])
    df['Width'] = df['Stop'] - df['Start']

    return df


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
    tt = f'Q-Q Plots: Genes Between Events\n{distance_title}'
    fig.suptitle(tt, fontsize=14)

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

    fig.set_tight_layout(True)
    plt.savefig(f'{distance_out}-qq.pdf')
    plt.close()

    # Test example plot
    mu = np.mean(qqdata["core"])
    poisson_array = stats.poisson.rvs(mu=mu, size=10000) # emp analog
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.set_title("Poisson vs Poisson Example Q-Q Plot", fontsize=14)
    qdist = stats.poisson(np.mean(poisson_array))
    sm.qqplot(poisson_array, dist=qdist, line="45", ax=ax)

    fig.set_tight_layout(True)
    plt.savefig(f'{distance_out}-qqex.pdf')
    plt.close()

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

    tt = f'Genes Between Events\n{distance_title}'
    fig.suptitle(tt, fontsize=14)

    if len(nc_genes) > minimum_genes:
        # calcuate distances, run poisson simulation, and ks test.
        emp, psn, mu, k, p = poisson_simulation(nc_genes, length)
        qqdata["nc_genes"] = emp
        _ = sns.histplot(x=emp, color=colors[4], ax=ax1, stat='density')
        _ = ax1.axvline(mu, color=colors[4], linestyle='dashed', linewidth=mlw)
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


def build_some_plots(df, genomes, pancats, outpre):
    
    colors = [
            #'#54278f', # purple for conserved genes
            '#ef6548', # brighter orange for conserved genes
            '#c51b7d', # bright pink recombining core
            '#4d9221', # dark green recombining accessory
            '#969696', # dark gray non-recombining
            '#bdbdbd', # neutral gray non-coding genome
            ]

    genome_out = f'{outpre}_genomes'

    fig, ax = plt.subplots(figsize=(7, 5))

    ax.set_xlabel('Gene location on contig (bp)')

    # set order of y-axis labels
    ylabel_set = []
    # get genome A and B contigs
    for genome, contigs in genomes.items():
        for i, (contig, length) in enumerate(contigs.items()):
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

            dfX = df[(df['Genome'] == genome) & (df['Contig'] == contig)]
            dfX.sort_values(by='Start', inplace=True, ignore_index=True)
            #dfX.to_csv(f'TEST_DF_{genome}.tsv', sep='\t')

            ## statics gene distance plot
            dist_title = f'Genome {genome} Contig {contig}'
            dist_out = f'{outpre}_{genome}_{contig}_distance'
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
    plt.savefig(f'{genome_out}.pdf')
    plt.close()

    return True


def build_some_plots_draftmode(df, genomes, pancats, outpre):
    
    colors = [
            '#88419d', # purple for highly conserved genes (HC) - og '#54278f'
            '#dd3497', # pink recombining core (RC) - og '#c51b7d' 
            '#41ab5d', # green recombining accessory (RA) - og '#4d9221'
            '#6baed6', # light blue non-recombining (NR)
            '#000000', # neutral gray non-coding genome
            ]

    genome_out = f'{outpre}_genomes'

    fig, ax = plt.subplots(figsize=(7, 5))

    ax.set_xlabel('Gene location in genome (bp)')

    # Define label order
    labord = [
            "gA-HC", "gB-HC", "gA-RC", "gB-RC",
            "gA-RA", "gB-RA", "gA-NR", "gB-NR"
            ]

    # For each label add pc label and blank barplot point.
    for lab in labord:
        genome = lab[1]
        length = sum(genomes[genome].values())
        # neutral gray for genome length
        ax.barh(lab, length, left=0, color=colors[4], height=0.75)
        #ax.barh(lab, 1, left=0, color='w', height=0.75, alpha=0)

    pc_switch = {'Core': 'RC', 'Conserved': 'HC', 'Accessory': 'RA'}

    for genome in genomes.keys():

        print(f'Computing and building plots for genome {genome} ...')

        # plot genes on the genome - select current genome (A or B)
        dfX = df[(df['Genome'] == genome)]
        dfX.sort_values(by='Start', inplace=True, ignore_index=True)
        dfX.to_csv(f'TEST_DF_{genome}.tsv', sep='\t')
        
        ## compute gene distances, run poisson simulation, KS test and plot
        dist_title = f'Genome {genome}'
        dist_out = f'{outpre}_{genome}_distance'
        length = sum(genomes[genome].values())
        _ = distance_plots(dfX, length, colors, dist_title, dist_out)

        genes = dfX['Gene'].to_numpy() # gene number
        Sta = dfX['Start'].to_numpy() # start position
        Wid = dfX['Width'].to_numpy() # Length or width of gene
        F10 = dfX['F100'].to_numpy() # F100 status
        Str = dfX['Strand'].to_numpy() # Strand
        #Pan = dfX['PanClass'].to_numpy() # Future - Pangenome category

        label = f'g{genome}'
    
        for G, S, W, F, Z in zip(genes, Sta, Wid, F10, Str):
            if F < 1:
                c = colors[3]
                pc = 'NR'
            elif F == 1:
                try:
                    pc = pc_switch[pancats[G]]
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
    plt.savefig(f'{genome_out}.pdf')
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
        '-o', '--output_file_prefix',
        help='Please specify a prefix for the output file names!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-draft', '--draft_genome_mode',
        help='Concatenates contigs into single contig. (Default=False)',
        metavar='',
        type=str,
        default=None,
        required=False
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
    draft = args['draft_genome_mode']

    # setup a switch
    switch = {0: 'A', 1: 'B'}

    # read in the data
    RBM = parse_aai_rbm(rbm) # values from the rbm file
    CDS = parse_prodigal_CDS(cA, cB, switch) # get CDS info for each genome
    genomes = parse_genome_fasta(gA, gB, switch) # Get genome contig lengths
    pancats = parse_pangenome_categories(PC)

    # process the data and build some plots
    if draft:
        print('\nOperating in draft genome mode ...\n')
        df = compute_gene_distances_draftmode(RBM, CDS, genomes, pancats)
        _ = build_some_plots_draftmode(df, genomes, pancats, outpre)
    else:
        print('\nOperating in closed genome mode ...\n')
        df = compute_gene_distances(RBM, CDS, pancats)
        _ = build_some_plots(df, genomes, pancats, outpre)


    print(f'\n\nComplete success space cadet!! Finished without errors.\n\n')
    
if __name__ == "__main__":
    main()


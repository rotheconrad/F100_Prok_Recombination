 #!/usr/bin/env python

'''This script is similar to the 03f script for genome pairs except it looks
at one genome compared to many genomes.

The input is a two column tsv file with paths to the genome fastas in the first
column and paths to the gene fastas in the second column. See the example input
file group_1_list.tsv. The first genome in the list will be compared to the
remaining genomes in the list.

And then it takes the same rbm, annotation, and pancat files as above.

This script returns data, graphics and statistics concerning the genomic
positions between genes categorized by highly conserved core genes,
recombining core genes, recombing accessory genes, and non-recombining
genes. A gene is said to be recombining if the RBM value is ≥ the -rec param
(default 99.8) and if it is not highly conserved.

This script performs some statistical tests on the distance between
recombinant genes, and on the distributions of gene annotations.

The -rec parameter (default 99.8) sets the threshold for an RBM to be
considered recombinant or not. It affects the results of the gene annotations
plot and chi-square hypothesis tests and it affects the recombinant positions
plot and Poisson and Geometric distribution tests.

see the github for a detail description of output files.
https://github.com/rotheconrad/F100_Prok_Recombination

# COG one letter code descriptions
# http://www.sbg.bio.ic.ac.uk/~phunkee/html/old/COG_classes.html

* SPECIAL NOTE ABOUT SPECIFIC gene pangenome classification:
The pangenome classification comes from MMSeqs2 clusters which filter more
strictly than the BLAST+ based RBM. The result is that some genes labeled with
the pangenome class of "Specific" based on singleton MMSeqs2 clusters do find
a RBM gene partner in other genomes. In the cases we observed this would be
because one gene of the RBM pair is about the half the length of other but with
high sequence identity.

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: March 2023
License :: GNU GPLv3
Copyright 2023 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse
from collections import defaultdict
from itertools import groupby
import pandas as pd; import numpy as np
pd.options.mode.chained_assignment = None
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy import stats
import seaborn as sns
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests


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


def parse_input_list(infile):
    ''' reads the input list, opens the file and returns:
    - mgenome is the main genome dict of {contig name: length}
    - genes is a nested dict for the other genomes {genome: {gene: gene info}}
    - geneinfo is a list of [start, stop, strand]
    - group is a dict of gene names from the other genomes in the group.
    '''
    print('\n\tCollecting input genomes and genes ...')

    # initialize return dicts
    mgenome, mgenes, group = {}, defaultdict((lambda: defaultdict())), {}
    # intialize gene file list
    group_files = []
    # genome order of input list
    gorder = []

    with open(infile, 'r') as file:
        # get the first line input files for the main genome and genes
        genome1, genes1 = file.readline().rstrip().split('\t')
        # parse the main genome fasta
        with open(genome1, 'r') as gfasta:
            for name, seq in read_fasta(gfasta):
                contigID = name[1:] # remove the ">" from the fasta defline
                seqlen = len(seq)
                mgenome[contigID] = seqlen
        # parse the main genes fasta
        with open(genes1, 'r') as gfasta:
            for name, seq in read_fasta(gfasta):
                X = name.split(' # ')
                # gene ID, '>' removed from fasta defline.
                geneID = X[0][1:]
                contigID = '_'.join(geneID.split('_')[:-1])
                # gene info
                start = int(X[1]) # gene start
                stop = int(X[2]) # gene stop
                strand = int(X[3]) # gene coding strand
                # add gene info to dict
                mgenes[contigID][geneID] = [start, stop, strand]
                # check all the genes are in order
                if start > stop:
                    print(
                        f'parse_prodigal_CDS ERROR: start position greater '
                        f'than stop positoin for: \n{name}'
                        )
        # Get the gene files names from the other genomes. Thats all we need here.
        for line in file:
            gene_file = line.rstrip().split('\t')[1]
            group_files.append(gene_file)
            with open(gene_file, 'r') as file:
                X = file.readline().split(' # ')
                genome = '_'.join(X[0].split('_')[:-2])[1:]
                gorder.append(genome)

    # Get the gene names from the group files. thats all we need here.
    for gene_file in group_files:
        with open(gene_file, 'r') as file:
            for name, seq in read_fasta(file):
                X = name.split(' # ')
                # gene ID, '>' removed from fasta defline.
                geneID = X[0][1:]
                group[geneID] = ''

    return mgenome, mgenes, group, gorder


def parse_rbm_file(rbm, rec, mgenes, group):
    ''' read the rbm file into a dict. Keep all rbms with the main genome genes
    and any of the group input genes. output dict of  {geneID: [rbm info]}
    where rbm info = [match, REC, pID, mismatch]
    Genes ≥ rec sequence similarity are set REC = 1 else REC is set to 0.
    '''
    print('\n\tReading RBM file ...')
    
    # initialize dictionary to store rbm data
    RBM = defaultdict(list)

    with open(rbm, 'r') as file:
        for line in file:
            X = line.rstrip().split('\t')
            geneA = X[0]
            contigA = '_'.join(geneA.split('_')[:-1])
            geneB = X[1]
            contigB = '_'.join(geneB.split('_')[:-1])
            pID = float(X[2])
            mismatch = int(X[4])

            REC = 1 if pID >= rec else 0

            # the rbm file should be the concatenated file from all genomes
            # This block should select only the input genomes from the main
            # genome and the group genomes.It should accept them in either
            # order, gA vs gB or gB vs gA.
            #tempgenes = {}
            #for contig, genes
            if contigA in mgenes and geneB in group:
                RBM[geneA].append([geneB, REC, pID, mismatch])
            elif contigB in mgenes and geneA in group:
                RBM[geneB].append([geneA, REC, pID, mismatch])

    return RBM


def parse_pangenome_file(pc):
    # parses the pangenome category file to dictionary of {gene: pancat}
    # where pancat equals Core, Accessory, or Specific.
    # Returns the dictionary

    print('\n\tReading PanCat file ...')

    pancats, repgenes = {}, {}
    
    with open(pc, 'r') as file:
        header = file.readline()
        for line in file:
            X = line.rstrip().split('\t')
            gene = X[0]
            repg = X[1]
            pancat = X[2]
            pancats[gene] = pancat
            repgenes[gene] = repg

    return pancats, repgenes


def parse_annotation_file(ano):
    
    # parses the annotation file into dict of {geneName: annotation}
    # autodetects EggNog or COGclassifier file
    # aggregates annotations into broad categories of:
    # Mobile, Metabolism, Conserved Hypothetical, Hypothetical, Other

    # define dictionary to hold output
    annos = {}
    # define dictionary to assign category for COGclassifier
    COGclass = {
            'X': 'Mobile', 'C': 'Metabolism 1', 'G': 'Metabolism 1',
            'E': 'Metabolism 1', 'F': 'Metabolism 1', 'H': 'Metabolism 1',
            'I': 'Metabolism 1', 'P': 'Metabolism 2', 'Q': 'Metabolism 2',
            'J': 'Ribosomal', 'A': 'Information', 'K': 'Information',
            'L': 'Information', 'B': 'Information', 'D': 'Cellular',
            'Y': 'Cellular', 'V': 'Cellular', 'T': 'Cellular',
            'M': 'Cellular', 'N': 'Cellular', 'Z': 'Cellular',
            'W': 'Cellular', 'U': 'Cellular', 'O': 'Cellular',
            'S': 'Conserved Hypothetical', 'R': 'Conserved Hypothetical'
            }

    # 'Hypothetical' will be added later for genes not in the annotation file

    # keywords to assign category for EggNog
    mobile = [
                'transposase', 'phage', 'integrase', 'viral', 'plasmid',
                'integron', 'transposon'
                ]

    with open(ano, 'r') as file:
        header = file.readline()
        if header[0] == 'Q':
            print('\n\tReading COGclassifier annotation file ...')
            for line in file:
                X = line.rstrip().split('\t')
                name = X[0]
                gene = X[5]
                desc = X[6]
                cog = X[7]
                cat = COGclass.get(cog, 'Other')
                annos[name] = [cat, cog, gene, desc]
        elif header[0] == '#':
            print('\n\tReading EggNog annotation file ...')
            for line in file:
                if line.startswith('#'): continue
                X = line.rstrip().split('\t')
                name = X[0] # representitive predicted gene name
                gene = X[8] # annotation short gene name
                desc = X[7] # annotation long gene name (description)
                cog = X[6][0] # select only first letter
                #same behavior as COGclassifier for multiple letter assignments
                # EggNog doesn't have cog X - mobile gene category
                # so we build it from keywords and the annotation description
                # check annotions and assign categories
                if any(mbl in desc.lower() for mbl in mobile):
                    cat = 'Mobile'
                elif cog == '-': cat = 'Hypothetical'
                else: cat = COGclass.get(cog, 'Other')
                # assigne the annotation to the gene
                annos[name] = [cat, cog, gene, desc]
        else:
            print('Annotation input file format error')

    return annos

###############################################################################
##### SECTION 02: CREATE DATA FRAME ###########################################
###############################################################################

def combine_input_data(mgenome, mgenes, RBM, pancats, repgenes, annos):
    
    # create dictionary to store info for dataframe
    D = {
        'Genome': [], 'Contig': [], 'Gene': [],
        'Match Genome': [], 'Match Contig': [], 'Match Gene': [],
        'PanCat': [], 'pID': [], 'REC': [],
        'Recombinant': [], 'Start': [], 'Stop': [], 'Strand': [],
        'COG Category': [], 'COG': [], 'Gene Annotation': [],
        'Annotation Description': [], 'Mismatch': []
        }
    # contig positions to mark
    contig_position = []
    # running genome length
    genome_length = 0
    for contig, contig_length in mgenome.items():
        genome = '_'.join(contig.split('_')[:-1])
        contig_genes = mgenes[contig]
        for gene, geneInfo in contig_genes.items():
            # set gene info, pancat, repgene, and annotation
            start = geneInfo[0] + genome_length
            stop = geneInfo[1] + genome_length
            strand = geneInfo[2]
            pc = pancats.get(gene, 'Specific')
            repg = repgenes.get(gene, 'na')
            default = ['Hypothetical', '-', '-', 'Hypothetical']
            ano = annos.get(repg, default)
            cat, cog, gn, desc = ano[0], ano[1], ano[2], ano[3]
            # RBM is {geneName: [list of rbm results for group]}
            # iterate through rbm results for each gene
            # if gene not in rbm group default to '-', 0, 0, 'x'
            # each rbm match will be a new row.
            rbm_matches = RBM.get(gene, [['-', 0, 0, 'x']])
            for rbmInfo in rbm_matches:
                match_gene = rbmInfo[0]
                if match_gene == '-': match_gene = '-_-_-'
                match_contig = '_'.join(match_gene.split('_')[:-1])
                match_genome = '_'.join(match_gene.split('_')[:-2])
                REC = rbmInfo[1]
                pID = rbmInfo[2]
                mismatch = rbmInfo[3]
                # assign recombinant or non-recombinant
                if REC == 1 and pc != 'Conserved':
                    rec = 'Recombinant'
                else: 
                    rec = 'Non-recombinant'

                # append data
                D['Genome'].append(genome)
                D['Contig'].append(contig)
                D['Gene'].append(gene)
                D['Match Genome'].append(match_genome)
                D['Match Contig'].append(match_contig)
                D['Match Gene'].append(match_gene)
                D['REC'].append(REC)
                D['pID'].append(pID)
                D['Mismatch'].append(mismatch)
                D['PanCat'].append(pc)
                D['COG Category'].append(cat)
                D['Recombinant'].append(rec)
                D['Start'].append(start)
                D['Stop'].append(stop)
                D['Strand'].append(strand)
                D['COG'].append(cog)
                D['Gene Annotation'].append(gn)
                D['Annotation Description'].append(desc)
                
        # increment genome
        genome_length += contig_length
        contig_position.append(genome_length)

    # convert dictionary to dataframe
    df = pd.DataFrame(D) #.sort_values(by=['Genome', 'Start'])
    # add gene lengths to Width column
    df['Width'] = df['Stop'] - df['Start']

    return df, contig_position

###############################################################################
##### SECTION 03: ANNOTATIONS HYPOTHESIS TESTING ##############################
###############################################################################

def plot_annotation_barplot(df, outpre):

    '''Plots annotation categories for recombinant vs non-recombinant genes'''

    print('\n\tBuilding annotation plot and tests ...')

    # set colors
    #colors = ['#c51b7d', '#e9a3c9', '#fde0ef', '#e6f5d0', '#a1d76a', '#4d9221']
    colors = [
            '#8dd3c7', '#ffffb3', '#bebada', '#fb8072',
            '#80b1d3', '#fdb462', '#b3de69', '#fccde5'
            ]

    # set annotation category order
    if 'Other' in df['COG Category'].unique():
        aorder = [
                'Other', 'Hypothetical', 'Conserved Hypothetical',
                'Ribosomal', 'Information', 'Cellular',
                'Metabolism 1', 'Metabolism 2', 'Mobile'
                ]
    else:
        aorder = [
                'Hypothetical', 'Conserved Hypothetical',
                'Ribosomal', 'Information', 'Cellular',
                'Metabolism 1', 'Metabolism 2', 'Mobile'
                ]

    xorder = ['Recombinant', 'Non-recombinant']
    # subset df
    adf = df.groupby(['Recombinant', 'COG Category'])['Gene'].count().unstack()
    adf = adf[aorder].T
    adf = adf[xorder].fillna(0)

    # categorical hypothesis testing with raw counts
    chip, pvals_corrected = annotation_hypothesis_test(adf)

    # calculate the percents of total per category
    total = adf.sum().to_list()
    ptots = [f'({round(i/sum(total) * 100, 2)}%)' for i in total]
    adf = adf.div(adf.sum(axis=0), axis=1).T

    fig, ax = plt.subplots(figsize=(7,5))

    ax = adf.plot.bar(stacked=True, ax=ax, color=colors, width=.7)
    # change axis labels
    ax.set_xlabel('')
    ax.set_ylabel("Gene fraction", fontsize=12)
    # set the axis parameters / style
    ax.minorticks_on()
    ax.tick_params(axis='both', labelsize=12)
    ax.tick_params(axis='x', labelrotation=0)
    ax.tick_params(axis='x', which='minor', bottom=False)
    # set grid style
    ax.yaxis.grid(
        which="minor", color='#f0f0f0', linestyle='--', linewidth=.75
        )
    ax.yaxis.grid(
        which="major", color='#d9d9d9', linestyle='--', linewidth=1
        )
    ax.set_axisbelow(True)

    # annotate percent of total gene difference
    ax.annotate(ptots[0], (0, 0.01), transform=ax.transAxes, ha='center')
    ax.annotate(ptots[1], (1, 0.01), transform=ax.transAxes, ha='center')

    # annotate individual percents and add asterisk if significant post hoc
    # double the corrected p value array since our plot has two columns
    sig = [i for i in pvals_corrected for _ in range(2)]
    for i, p in enumerate(ax.patches):
        width, height = p.get_width(), p.get_height()
        x, y = p.get_xy()
        # don't print categories with 0 genes on the bar plot
        if height == 0:
            continue
        # add asterisk if significant of alpha 0.05
        line = f'*{height:.2f}*' if sig[i] < 0.05 else f'{height:.2f}'
        ax.text(x+width/2, 
                y+height/2, 
                line, 
                horizontalalignment='center', 
                verticalalignment='center')

    # create legend to separate file
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(
                reversed(handles),
                reversed(labels),
                ncol=1,
                loc='center left',
                frameon=False,
                fontsize=12,
                bbox_to_anchor=(1, 0.5)
                )

    # adjust layout, save, and close
    fig.set_tight_layout(True)
    fig.savefig(f'{outpre}_annotations_bar.pdf')
    plt.close() 

    return True

def annotation_hypothesis_test(adf):
    ''' performs chi square and post hoc tests between annotation
    categorgies for recombining vs non-recombining genes.
    '''
    
    # create and print contingency table
    ctab = adf.copy(deep=True)
    ctab['Total'] = ctab.sum(axis=1)
    ctab.loc['Total'] = ctab.sum(axis=0)
    print(f'\nInitial Chi2 test contingency table:\n\n{ctab}')

    # chi2 test on the full data
    chi2, chip, dof, ex = stats.chi2_contingency(adf, correction=True)
    # create expected frequency table
    efreq = pd.DataFrame(ex, index=adf.index, columns=adf.columns)
    efreq['Total'] = efreq.sum(axis=1)
    efreq.loc['Total'] = efreq.sum(axis=0)
    print(f'\nChi2 expected frequency table:\n\n{efreq}')
    # print test statitic and p value
    print(f'\nchi2 statistic: {chi2:.4f}, dof: {dof}, chi2 pvalue: {chip:.6f}')

    # perform post hoc test on combinations if significant (< 0.05)
    if chip < 0.05:
        pvals, pvals_corrected = post_hoc_test(adf)
        print('\nPost hoc p values:\n', pvals)
        print('\nBenjamini/Hochberg corrected p values:\n', pvals_corrected)

    else:
        pvals_corrected = [1] * len(adf)

    return chip, pvals_corrected


def post_hoc_test(adf):
    ''' loops through individual rows and performs chi2 post hoc '''
    pvals = []

    bdf = adf.T
    for name in bdf.columns:
        xdf = bdf.drop(name, axis=1)
        xdf['OTHERs'] = xdf.sum(axis=1)
        xdf[name] = bdf[name]
        ydf = xdf[['OTHERs', name]].T
        c, p, d, x = stats.chi2_contingency(ydf, correction=True)
        # create expected frequency table
        extab = pd.DataFrame(x, index=ydf.index, columns=ydf.columns)
        extab['Total'] = extab.sum(axis=1)
        extab.loc['Total'] = extab.sum(axis=0)
        # create post hoc test contingency table
        ydf['Total'] = ydf.sum(axis=1)
        ydf.loc['Total'] = ydf.sum(axis=0)
        # print post hoc info
        print(f'\nPost hoc Chi2 test contingency table for {name}:\n')
        print(ydf)
        print(f'\nChi2 expected frequency table:\n\n{extab}')
        print(f'\nchi2 statistic: {c:.4f}, dof: {d}, chi2 pvalue: {p:.6f}')

        pvals.append(p)

        reject_list, pvals_corrected = multipletests(pvals, method='fdr_bh')[:2]

    return pvals, pvals_corrected

###############################################################################
##### SECTION 04: GENE RBM IDENTITY VS GENOME POSITION ########################
###############################################################################

def ani_sort(dfG, gorder):
    ''' sorts the genomes by ANI. Gets the ANI from RBMs '''

    data = {'Match Genome': [], 'ANI': []}

    for genome in gorder:
        dfA = dfG[dfG['Match Genome'] == genome]
        rbm_ani = dfA[dfA['pID'] != 0]['pID'].mean()
        data['Match Genome'].append(genome)
        data['ANI'].append(rbm_ani)

    dfANI = pd.DataFrame(data)
    dfANI = dfANI.sort_values('ANI', ascending=False)
    gorder = dfANI['Match Genome'].tolist()

    return gorder


def build_pos_line_plot(df, mgenome, outpre, cpos, rec, gorder):
    ''' plots gene sequence identity on y axis vs genome coords on x axis.
        cpos is a list of contig lengths to add markers '''

    dfG = df[df['Match Genome'] != '-']
    dfG['Mid'] = dfG[['Start', 'Stop']].mean(axis=1)

    # Sort by ANI
    gorder = ani_sort(dfG, gorder)

    colors = {
            'Conserved': '#e41a1c', # red high conserved core gene
            'Core': '#377eb8', # blue core gene 
            'Accessory': '#4daf4a', # green accessory gene
            'Specific': '#984ea3', # purple genome specific gene
            }

    print(f'\n\tBuilding gene line plots ...')

    # set current genome length / 4, rounded up to whole integer
    glength = sum(mgenome.values())
    # set yaxis max
    ymax = 100
    
    group_genomes = gorder #dfG['Match Genome'].unique()
    subs = len(group_genomes)

    # initialize the figure
    fig, axs = plt.subplots(subs, 1, figsize=(70, subs*2))

    for i, (ax, genome) in enumerate(zip(axs, group_genomes)):
        print(f'\t\tPlotting genome {genome}')
        xmin, xmax = 0, glength
        dfS = dfG[dfG['Match Genome'] == genome]
        # calculate non zero mean. ANI of RBMs
        rbm_ani = dfS[dfS['pID'] != 0]['pID'].mean()
        # set variable yaxis min
        ymin = 95 if rbm_ani >= 97.5 else 90
        x = dfS['Mid'].to_list()
        # replace values < ymin with ymin
        ids = dfS['pID'].to_numpy()
        y = np.where(ids < ymin, ymin, ids).tolist()
        c = [colors[i] for i in dfS['PanCat'].to_list()]

        # plot the data
        ax.scatter(x, y, color=c, marker='|', linestyle='-')

        # add genome name to plot
        gname = f'genome: {genome} | ANI: {rbm_ani:.4f}'
        ax.text(xmin, ymin-1.2, gname, ha='left', va='top', fontsize=8)

        # set axis limits
        ax.set_xlim(xmin-0.5, xmax+0.5)
        ax.set_ylim(ymin, ymax)
        #ax.set_xticks([])
        #ax.set_xticklabels([])

        # add dashed lines for recombinant threshold and ANI

        ax.hlines(y=rec, xmin=xmin, xmax=xmax, lw=1, colors='k', ls=':')
        ax.hlines(y=rbm_ani, xmin=xmin, xmax=xmax, lw=1, colors='k', ls='--')

        # add contig markers
        for mark in cpos[:-1]:
            ax.text(mark, ymax,"l", color='#969696',
                    fontsize=8, ha='center', va='bottom')

    # build the legend
    legend_labels = ['Conserved', 'Core', 'Accessory', 'Specific']
    legend_elements = []

    for x in legend_labels:
        n = Line2D(
            [0], [0], color='w', label=x, marker='s',
            markersize=15, markerfacecolor=colors[x]
            )
        legend_elements.append(n)

    n = Line2D(
        [0], [0], color='k', label='Recombinant', marker=None, linestyle=':'
        )
    legend_elements.append(n)

    n = Line2D(
        [0], [0], color='k', label='RBM ANI', marker=None, linestyle='--'
        )
    legend_elements.append(n)

    axs[0].legend(
        handles=legend_elements,
        fontsize=12,
        fancybox=True,
        framealpha=0.0,
        frameon=False,
        loc='lower center',
        bbox_to_anchor=(0, 0.98, 1, 0.2),
        ncol=6
        )

    fig.set_tight_layout(True)
    fig.subplots_adjust(hspace=0)
    plt.savefig(f'{outpre}_posline.pdf')
    plt.close()

    return True

###############################################################################
##### SECTION 05: RECOMBINANT GENE POSITION BY PANGENOME CLASS ################
###############################################################################

def build_pos_bar_plots(df, mgenome, pancats, outpre, cpos):

    print(f'\n\tBuilding recombinant position plots ...')

    # Builds gene/genome position bar style plot with pangenome categories
    # Runs evenness of recombinant gene positions tests
    # cpos is contig position array to mark them on the plot.
    
    colors = [
            '#fee391', # yellow for highly conserved genes (HC)
            '#dd3497', # pink recombining core (RC) 
            '#41ab5d', # green recombining accessory (RA)
            '#6baed6', # light blue non-recombining (NR)
            '#f03b20', # red genome specific genes (GS)
            '#000000', # neutral gray non-coding genome
            ]

    # outfile
    genome_out = f'{outpre}_posbar'

    # genome params
    glength = sum(mgenome.values())
    total_genes = len(df)

    # initialize fig and ax
    fig, ax = plt.subplots(figsize=(7, 5))
    # xlabel
    ax.set_xlabel('Gene location in genome (bp)')
    # Define label order
    labord = ["HC", "RC", "RA", "NR", "GS"]
    # For each label add pc label and blank barplot point.
    for lab in labord:
        pclass = lab
        # black for full genome length background
        ax.barh(lab, glength, left=0, color=colors[5], height=0.75)
        # Add percent of total genes to end of each bar
        # this is probably an ugly way to do this that i've added at the end.
        total_genes = len(df)
        # iterate over each bar and annotate the value 
        if pclass == 'HC':
            # subset conserved genes 
            conserved = df[df['PanCat'] == 'Conserved']
            tc = colors[0]
            # percent of total genes
            tl = round(len(conserved) / total_genes * 100, 2)
            # mutation rate
            conserved = conserved[conserved['Mismatch'] != 'x']
            mr = round(conserved['Mismatch'].sum() / conserved['Width'].sum() * 100, 4)
        elif pclass == 'RC':
            # subset recombinant core genes
            core = df[(df['PanCat'] == 'Core') & (df['REC'] == 1)]
            tc = colors[1]
            tl = round(len(core) / total_genes * 100, 2)
            core = core[core['Mismatch'] != 'x']
            mr = round(core['Mismatch'].sum() / core['Width'].sum() * 100, 4)
        elif pclass == 'RA':
            # subset recombinant accessory genes
            accessory = df[(df['PanCat'] == 'Accessory') & (df['REC'] == 1)]
            tc = colors[2]
            tl = round(len(accessory) / total_genes * 100, 2)
            accessory = accessory[accessory['Mismatch'] != 'x']
            mr = round(accessory['Mismatch'].sum() / accessory['Width'].sum() * 100, 4)
        elif pclass == 'NR':
            # subset non conserved and non recombining genes
            nc_genes = df[
                            (df['PanCat'] != 'Conserved') & 
                            (df['PanCat'] != 'Specific') & 
                            (df['REC'] != 1)
                            ]
            tc = colors[3]
            tl = round(len(nc_genes) / total_genes * 100, 2)
            nc_genes = nc_genes[nc_genes['Mismatch'] != 'x']
            mr = round(nc_genes['Mismatch'].sum() / nc_genes['Width'].sum() * 100, 4)
        elif pclass == 'GS':
            # subset genome specific genes
            specific = df[df['PanCat'] == 'Specific']
            tc = colors[4]
            tl = round(len(specific) / total_genes * 100, 2)
            specific = specific[specific['Mismatch'] != 'x']
            if len(specific) >= 1:
                mr = round(specific['Mismatch'].sum() / specific['Width'].sum() * 100, 4)
            else:
                mr = 'n/a'
        # add it to the plot
        ax.text(glength+(glength/100), lab, f'{tl}%\n{mr}%', color=tc, va='center')

    # dict lookup for pangenome labels
    pc_switch = {'Core': 'RC', 'Conserved': 'HC', 'Accessory': 'RA', 'Specific': 'GS'}

    ## compute gene distances, run poisson simulation, KS test and plot
    dist_out = f'{outpre}_gene_distribution'
    # changed doing by pangenome category to doing 1 for each genome.
    # and collapsing consecutive REC genes to 1 event
    # returns percent recently recombining genes
    _ = distance_plots_2(df, colors, dist_out)

    genes = df['Gene'].to_numpy() # gene number
    pcat = df['PanCat'].to_numpy() # pancat
    Sta = df['Start'].to_numpy() # start position
    Wid = df['Width'].to_numpy() # Length or width of gene
    F10 = df['REC'].to_numpy() # REC status
    Str = df['Strand'].to_numpy() # Strand

    # deduplicate points - don't plot over positions already plotted
    dd = {'RC': {}, 'HC': {}, 'RA': {}, 'GS': {}, 'NR': {}}

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

        # don't plot over a position thats already been plotted.
        if S in dd[pc]:
            continue
        # plot position if it has not been plotted.
        else:
            ylab = pc
            ax.barh(ylab, W, left=S, color=c, height=0.75, alpha=0.3)
            # store positions that have already been plotted.
            dd[pc][S] = ''

    # add contig markers
    mark_pos = cpos[:-1]
    ypos = [-0.5] * len(mark_pos)

    ax.plot(mark_pos, ypos, marker="|", linestyle="", color='#969696')

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

    ''' A zero here indicates adjacent genes with 100% sequence identity to
    their respective RBMs which is likely a signal from a single recombination
    event. Likewise, a 1 indicates their is only 1 non 100% sequence identity
    gene between two events which could also indicate a single recombination
    event. The array show the number of genes between each recombination event,
    so if we drop 0's or 1's or both from the array, this is the same as 
    combining consecutive 0's and 1's into a single event. '''
    #emperical_dist = temp_dist # retain all distances
    # collapse consectuve 0's into a single event
    emperical_dist = temp_dist[temp_dist != 0]
    # or collapse consecutive 0's and 1's into a single event.
    #emperical_dist = temp_dist[np.logical_and(temp_dist != 0, temp_dist != 1)]

    ''' print out value counts of distances between events
    values, counts = np.unique(emperical_mean_dist,return_counts=True)
    for v, c in zip(values, counts):
        print(v, c)
    '''

    mu =  np.mean(emperical_dist) # genes per kilo base pair
    poisson_array = stats.poisson.rvs(mu=mu, size=10000)
    # decided on a sample size of 10,000 instead of array length or gene count

    # Kolmogorov-Smirnov test emperical data against poisson distribution
    # returns ks_statistic, p_value
    k, p = stats.kstest(emperical_dist, poisson_array)

    # Unit test / sanity check use poisson array instead of input array
    # k is small and p is large when the distributions are likely the same
    # fail to reject null hypothesis that the distributions are different.
    # poisson_A = stats.poisson.rvs(mu=mu, size=10000)
    # poisson_B = stats.poisson.rvs(mu=mu, size=10000)
    # k, p = stats.kstest(poisson_A, poisson_B)
    # emperical_mean_dist = poisson_A

    return emperical_dist, poisson_array, mu, k, p


def geom_dist(array, n):

    # each gene is an event. 100% RBM is a success
    # array should contain indexed positions of all 100% RBM genes.
    # n should be the total genes in the genome.

    # make sure the array is sorted in ascending order.
    sorted_array = np.sort(array)

    # get the difference between gene index positions (genes between events)
    # and subtract 1 so the minimum distance is 0 instead of 1.
    geom_emp = np.diff(sorted_array) - 1

    # Geometric distributions are built around the probability of success
    # in this case the input array is the number of successes (recombinations)
    # and n is the total number of events (genes in the genome)
    # so p = len(array) / n
    p = len(array) / n
    geom_array = stats.geom.rvs(p=p, size=10000)

    # Kolmogorov-Smirnov test emperical data against geometric distribution
    # returns ks_statistic, p_value
    geom_k, geom_p = stats.kstest(geom_emp, geom_array)
    
    return geom_emp, geom_array, p, geom_k, geom_p


def distance_plots_2(df, colors, distance_out):

    # Calculate distance between REC gene events in the genome
    # Collapse consecutive REC genes to single event.
    # simulate poisson distribution
    # run kolmogorov-smirnov test
    # plot results
    # changed from contig length in bp to simply the length of the df
    total_genes = len(df) # total events for poisson simulation (total genes)
    # the dataframe index is the gene position in the genome
    # distance between genes is distance between index and the preceeding index 
    # subset non conserved genes ie disregard conserved genes
    # subset recombinant all genes that are not "conserved"
    recomb_genes = df[(df['REC'] == 1) & (df['PanCat'] != 'Conserved')].index

    # percent of genes that are recently recombining
    # 100% RBMs / total genes

    # skip statistical tests if fewer than 5 genes.
    if len(recomb_genes) < 5:
        print('Fewer than 5 recombinant genes. Skipping statistical test.')
        return False

    # calcuate distances, run poisson simulation, and ks test.
    emp, psn, mu, k, p = poisson_simulation(recomb_genes, total_genes)

    geom_emp, gmt, p, geom_k, geom_p = geom_dist(recomb_genes, total_genes)

    tx, ty = 0.72, 0.82 # stats test text position
    pcol = 'r' # poisson model kde line plot color
    mlw = 1 # kdeplot line width
    axlabsz = 12 # ax label sizes
    ts = 8 # stat test text size
    bw, ct = 3, 0 # kdeplot bw_adjus t and cut params

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(7, 7))

    # plot emp distribution, poissoin distribution and stats test results
    _ = sns.histplot(x=emp, bins=30, color=colors[5], ax=ax1, stat='density')
    _ = ax1.axvline(mu, color=colors[5], linestyle='dashed', linewidth=mlw)
    _ = sns.kdeplot(x=psn, color=pcol, ax=ax1, cut=ct, bw_adjust=bw)
    line = f'k = {k:.4f}\np = {p:.4f}'
    _ = ax1.text(tx, ty, line, transform=ax1.transAxes, fontsize=ts)
    ax1.set_ylabel('Density', fontsize=axlabsz)
    ax1.set_xlabel(' ')

    # plot emp distribution, geometric distribution and stats test results
    _ = sns.histplot(x=geom_emp, bins=30, color=colors[5], ax=ax2, stat='density')
    _ = sns.kdeplot(x=gmt, color=pcol, ax=ax2, cut=ct, bw_adjust=bw)
    line = f'k = {geom_k:.4f}\np = {geom_p:.4f}'
    _ = ax2.text(tx, ty, line, transform=ax2.transAxes, fontsize=ts)
    ax2.set_ylabel('', fontsize=axlabsz)
    ax2.set_xlabel('')

    # plot Q-Q plot emperical vs theoretical poisson
    qdist = stats.poisson(mu)
    sm.qqplot(emp, dist=qdist, line="45", ax=ax3)
    ax3.set_ylabel('Sample quantiles', fontsize=axlabsz)
    ax3.set_xlabel('')

    # plot Q-Q plot emperical vs theoretical geomitric
    qdist = stats.geom(p)
    sm.qqplot(geom_emp, dist=qdist, line="45", ax=ax4)
    ax4.set_ylabel('', fontsize=axlabsz)
    ax4.set_xlabel(' ')

    fig.text(0.55, 0.51, 'Genes between recombinant events', ha='center', fontsize=axlabsz,)
    fig.text(0.55, 0.02, 'Theoretical quantiles', ha='center', fontsize=axlabsz,)

    fig.set_tight_layout(True)
    plt.savefig(f'{distance_out}.pdf')
    plt.close()

    return True

###############################################################################
##### SECTION 06: Build binary matrix for recombinant RBMs ####################
###############################################################################

def build_rbm_binary_matrix(df, outpre):
    ''' generate and write out a binary matrix for RBM rarefaction analysis
    RBMs ≥ rec are assigned a 1 and RBMs < rec are assigned a 0
    Conserved genes are assigned a 2.
    Columns are genomes in the input group
    Rows are RBMs in gene order of the main (1st) genome in the input group
    '''

    print(f'\n\tBuilding binary recombinant RBM matrix ...')

    # create ternary option for conserved genes assign 2
    df.loc[df['PanCat'] == 'Conserved', 'REC'] = 2

    # Due to retaining tied RBMs, we need to select the tied RBM with the
    # greatest pID and retain only one RBM for the heatmap.
    df['Start-Match'] = df['Start'].astype(str) + '-' + df['Match Genome'].astype(str)
    dfX = df[df.duplicated(subset='Start-Match', keep=False) == True]
    dfX = dfX[['Start', 'pID', 'Match Gene', 'Gene', 'REC', 'Start-Match']]
    dfX = dfX.sort_values(by='pID', axis=0, ascending=False)
    dfX = dfX[dfX.duplicated(subset='Start-Match', keep="first") == True]
    df = df.drop(dfX.index)

    # Select columns. Start position is unique for each gene because we
    # concatenate the contigs and adjust the start relative to the concatenated
    # genome length when we parse the genome and gene fasta files.
    selection = ['Start', 'Match Genome', 'REC']
    bdf = df[selection]
    # turn our three column dataframe into a matrix (not square)
    # rows are genes represented by start position
    # columns are genomes represented by match genome name
    # values are 0 or 1 based on the rec threshold (default 99.8)
    matrix = bdf.pivot(index='Start', columns='Match Genome', values='REC')
    # drop the multi-index name,
    # drop the '-' column that arises from genome specific genes
    # fill nan's with 0s. these are genes in the reference genome without RBMs
    # in the query genome so they are non recombinant with the query genome.
    # this causes an error if there are not any genome specific genes.
    # future update, add logic to check. For now use try: except:
    try:
        matrix = matrix.rename_axis(None, axis=0).drop('-', axis=1).fillna(0)
    except:
        print(
            "\n\n\t\tLine 1043 Note: No genome specific genes ('-') identified"
            ".\n\t\tThis is not an error. Only a note. Everything is good.\n\n"
            )
    # write to file.
    matrix.to_csv(f'{outpre}_rbm_matrix.tsv', sep='\t', index=False)

    return True

###############################################################################
##### SECTION 07: CALCULATE RECOMBINATION RATES ###############################
###############################################################################

def get_recombination_rate(df, window, outfile, rbm_ani):

    # Calculate recombination vs mutation rate

    # get sums, define variable
    rec_mismatch = df['Rec_Mismatch'].sum()
    rec_length = df[df['Rec_Mismatch'] > 0]['Width'].sum()
    norec_mismatch = df['NoRec_Mismatch'].sum()
    norec_length = df[df['NoRec_Mismatch'] > 0]['Width'].sum()
    total_length = df['Width'].sum()
    # skip windows with zero mismatches
    if rec_mismatch == 0 or norec_mismatch == 0: return 0, 0, 0
    # cacluate mutation rates as (mismatches / gene length)
    divergence_time = 0.2
    Rp_time = (100 - rbm_ani - divergence_time) / 100
    Rm_time = divergence_time / 100
    Rmr = rec_mismatch / rec_length # recent mutation rate
    Pmr = norec_mismatch / norec_length # past mutation rate
    # calculate recent mutations and recently removed mutations
    Rm = Rmr * total_length / Rm_time # total recent mutations
    Rp = Pmr * rec_length / Rp_time # recently removed mutations
    # get the ratio
    slide1 = Rp / Rm if Rm > 0 else 0

    # slide2
    pm2_rate = (100 - rbm_ani) / 100
    pm2 = pm2_rate * rec_length
    rm2_rate = divergence_time / 100
    rm2 = rm2_rate * total_length
    slide2 = pm2 / rm2 if rm2 > 0 else 0

    # slide 3
    df['drate'] = (divergence_time - (df['NoRec_Mismatch'] + df['Rec_Mismatch'] / df['Width'] ))
    pm3 =  df['drate'].mean() * rec_length
    tm3 = (divergence_time / 100) * total_length
    am3 = tm3 - pm3
    slide3 = pm3 / am3

    lineout = (
                f'{window}\t{rec_mismatch}\t{rec_length}\t{norec_mismatch}\t'
                f'{norec_length}\t{total_length}\t{Rmr}\t{Pmr}\t{Rm}\t{Rp}\t'
                f'{pm2_rate}\t{pm2}\t{rm2_rate}\t{rm2}\t{pm3}\t{tm3}\t{am3}\t'
                f'{slide1}\t{slide2}\t{slide3}\n'
                )
    outfile.write(lineout)

    return slide1, slide2, slide3


def recombination_rate_plots(df, outpre):

    # use only RBM genes, drop genes without RBM
    dfx = df[df['Mismatch'] != 'x'] # drop genome specific genes
    rbm_ani = dfx[dfx['pID'] != 0]['pID'].mean()
    # Convert one line per gene
    dfy = dfx[['REC', 'Start', 'Mismatch', 'Width']]
    dfy = dfy.groupby(['Start', 'REC', 'Width'], as_index=False)['Mismatch'].sum()
    df0 = dfy[dfy['REC'] == 0].set_index('Start')
    df0 = df0.rename(columns={'Mismatch':'NoRec_Mismatch'})
    df0 = df0.drop('REC', axis=1)
    df1 = dfy[dfy['REC'] == 1].set_index('Start')
    df1 = df1.rename(columns={'Mismatch':'Rec_Mismatch'})
    df1 = df1.drop('REC', axis=1)
    dfZ = df0.join(df1, lsuffix='_0', rsuffix='_1', how='outer').reset_index()
    dfZ['Width'] = dfZ[['Width_0', 'Width_1']].max(axis=1)
    # open output file for rec rate table
    recout = open(f'{outpre}_rec_rate_table.tsv', 'w')
    header = (
                'Window\tRecombinant Mismatches\tRecombinant Gene Length\t'
                'Non-recombinant Mismatches\tNon-recombinant Gene Length\t'
                'Total Gene Length\tRmr * Gl\tPmr * Rl\t'
                'Rm\tRp\tPm2 rate\tPm2\tRm2 rate\tRm2\tPm3\tTm3\tAm3\t'
                'slide1\tslide2\tslide3\n'
                )
    recout.write(header)
    # get full genome recombination rate
    # Rmr = Number of recent mutations across genomic length
    # Pmr = Number of recent mutations purged by recombination
    slide1, slide2, slide3 = get_recombination_rate(dfZ, 'Full genome', recout, rbm_ani)
    # store the title for the plot
    titles = {
            'slide1': f'Slide 1: {slide1:.2f}',
            'slide2': f'Slide 2: {slide2:.2f}',
            'slide3': f'Slide 3: {slide3:.2f}'
            }

    # get sliding window recombination rates
    # dict to store results from sliding window
    data = defaultdict(list)
    # wdinows: number of genes, to iterate
    windows = [100, 200, 300, 400, 500, 600, 700, 800]
    # full length of the data frame
    dflen = len(dfZ)
    # iterate
    for window in windows:
        for start in range(dflen-window+2):
            stop = start+window-1
            dfW = dfZ.loc[start:stop]
            dfM = dfx[(dfx['Start'] >= start) & (dfx['Stop'] <= stop)]
            rbm_ani = dfM[dfM['pID'] != 0]['pID'].mean()
            slide1, slide2, slide3 = get_recombination_rate(dfW, window, recout, rbm_ani)
            data[f'{window}_slide1'].append(slide1)
            data[f'{window}_slide2'].append(slide2)
            data[f'{window}_slide3'].append(slide3)

    # close recout data table
    recout.close()

    # Begin building the plots
    axnum = len(windows)
    # text position xy, text size, axis label size
    tx, ty, ts, ls = 0.01, 0.95, 8, 12

    fig, axes = plt.subplots(axnum, 3, figsize=(15,24))
    axs = axes.reshape(-1)

    for i, (lab, y) in enumerate(data.items()):
        x = range(len(y))
        axs[i].plot(x, y)
        window = lab.split('_')[0]
        if i in [0,1,2]:
            title = titles[lab.split('_')[1]]
            axs[i].set_title(f'{title}')

        axs[i].hlines(y=1, xmin=min(x), xmax=max(x), lw=1, colors='k', ls=':')
        axs[i].set_yscale('symlog')

        axs[i].set_xlabel(f'{window} gene step')
        axs[i].set_ylabel('')

    fig.set_tight_layout(True)
    plt.savefig(f'{outpre}_Rec_rate_windows.pdf')
    plt.close()

    return True


def length_of_recombination_events(df, outpre):

    # Calculate length of consecutive REC gene events in the genome
    # plot histogram of results
    # count number of 1's before a 0
    rec_list = df['REC'].to_list()
    rec = []
    norec = []
    # groupby groups consecutive values into value, group
    # where group contains a value of each consecutive number/character
    for value, group in groupby(rec_list):
        count = sum(1 for _ in group)
        rec.append(count) if value == 1 else norec.append(count)

    axlabsz = 12 # ax label sizes

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7, 7))

    # plot emp distribution, poissoin distribution and stats test results
    _ = sns.histplot(x=rec, bins=30, color='#000000', ax=ax1, stat='frequency')
    ax1.set_ylabel('Frequency', fontsize=axlabsz)
    ax1.set_xlabel('Consecutive recombinant genes')

    # plot emp distribution, poissoin distribution and stats test results
    _ = sns.histplot(x=norec, bins=30, color='#000000', ax=ax2, stat='frequency')
    ax2.set_ylabel('Frequency', fontsize=axlabsz)
    ax2.set_xlabel('Consecutive non-recombinant genes')

    fig.set_tight_layout(True)
    plt.savefig(f'{outpre}_Adjecent_Rec_Length.pdf')
    plt.close()

    # write data to output file to compare between species or groups
    outfile = f'{outpre}_Adjecent_Rec_Length.tsv'
    with open(outfile, 'w') as out:
        recout = [str(i) for i in rec]
        out.write('Consecutive Recombinant,' + ','.join(recout) + '\n')
        nreco = [str(i) for i in norec]
        out.write('NonConsecutive Non-recombinant,' + ','.join(nreco) + '\n')

    return True

###############################################################################
##### SECTION 08: FIGURE X ###########################################
###############################################################################

# empty space for future ideas.

###############################################################################
###############################################################################

def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-i', '--input_list',
        help='Please specify the input list (2 column tsv file)!',
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
        '-rbm', '--RBM_allvall_file',
        help='Please specify the all vs all RBM file!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-pc', '--pangenome_categories',
        help='Please specify the tsv from 04d_Get_Genes_Clusters_PanCat.py!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-ano', '--annotation_file',
        help='Specify the representative gene annotation file.',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-rec', '--recombination_cutoff',
        help='(OPTIONAL) Specify recombination cutoff (default = 99.8).',
        metavar='',
        type=float,
        default=99.8,
        required=False
        )
    args=vars(parser.parse_args())

    # Do what you came here to do:
    print('\n\nRunning Script...\n')

    # define parameters
    infile = args['input_list']
    rbm = args['RBM_allvall_file']
    pc = args['pangenome_categories']
    outpre = args['output_file_prefix']
    ano = args['annotation_file']
    rec = args['recombination_cutoff']

    ## SECTION 01: Parse input data
    # parse input list.
    # mgenome is the main genome dict of {contig name: length}
    # mgenes is the main geneome dict of {gene name: geneinfo}
    # geneinfo is a list of [start, stop, strand]
    # group is a dict of gene names from the other genomes in the group
    mgenome, mgenes, group, gorder = parse_input_list(infile)
    # parse RBM file
    RBM = parse_rbm_file(rbm, rec, mgenes, group)
    # parse pancat file
    pancats, repgenes = parse_pangenome_file(pc)
    # parse annotation file
    annos = parse_annotation_file(ano)

    ## SECTION 02: Create data frame
    # this step takes all the input we just parsed and creates a dataframe
    df, cpos = combine_input_data(mgenome, mgenes, RBM, pancats, repgenes, annos)
    # write df to file
    group_df_file = f'{outpre}_group_data.tsv'
    df.to_csv(group_df_file, sep='\t', index=False)

    ## SECTION 03: Annotations Hypothesis testing
    # partition into recombinant genes vs non-recombinant REC ≥ REC
    _ = plot_annotation_barplot(df, outpre)

    ## SECTION 04: gene RBM identity vs. genome position
    _ = build_pos_line_plot(df, mgenome, outpre, cpos, rec, gorder)

    ## SECTION 05: recombinat gene positions by pangenome class
    _ = build_pos_bar_plots(df, mgenome, pancats, outpre, cpos)

    ## SECTION 06: Build binary matrix for recombinant RBMs
    # generate and write out a binary matrix for RBM rarefaction analysis
    # RBMs ≥ rec are assigned a 1 and RBMs < rec are assigned a 0
    # Columns are genomes in the input group
    # Rows are RBMs in gene order of the main (1st) genome in the input group
    _ = build_rbm_binary_matrix(df, outpre)

    ## SECTION 07: Calculate recombination rates
    print('\n\tCalculating recombination lengths ...')
    _ = length_of_recombination_events(df, outpre)
    # did not finish this function for the group genomes
    #_ = recombination_rate_plots(df, outpre)

    ## SECTION 08: empty space for future ideas.
    # this step creates x figure and does x stats
    # mgenome for length of genome and contigs

    print(f'\n\nComplete success space cadet!! Finished without errors.\n\n')

if __name__ == "__main__":
    main()


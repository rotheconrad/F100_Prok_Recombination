 #!/usr/bin/env python

'''Create RBM position line plots with more control.

This script is inteded to make polished final figures from the posline RBM
data.

The code is optimized for 10 query genomes + the 1 reference genome and
a 1 Mbp genome section.

The input is a two column tsv file with paths to the genome fastas in the first
column and paths to the gene fastas in the second column. See the example input
file group_1_list.tsv. The first genome in the list will be compared to the
remaining genomes in the list.

And then it takes the same rbm, annotation, and pancat files as above.

x-axis coordinates are in base pairs but the figure axis is in mega bp.
If you want to plot 1Mbp - 2Mbp enter
-xmin 1000000 -max 2000000

y-axis coordinates are in ANI %. Any RBM genes outside of the ymin and ymax
are plotted at the ymin and ymax right on the line.


see the github for a detail description of output files.
https://github.com/rotheconrad/F100_Prok_Recombination


-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: Dec 2023
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


def parse_metadata_file(metadata):
    # reads metadata into md dict of {genome name: phylogroup}
    # return md dict
    md = {}

    with open(metadata, 'r') as file:
        for line in file:
            X = line.rstrip().split('\t')
            genome = X[0]
            phylogroup = X[1]
            genomovar = X[2]
            md[genome] = [phylogroup, genomovar]

    return md

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
##### SECTION 04: GENE RBM IDENTITY VS GENOME POSITION ########################
###############################################################################

def ani_sort(dfG, gorder):
    ''' sorts the genomes by ANI. Gets the ANI from RBMs '''

    data = {'Match Genome': [], 'ANI': []}
    ani = {}

    for genome in gorder:
        dfA = dfG[dfG['Match Genome'] == genome]
        rbm_ani = dfA[dfA['pID'] != 0]['pID'].mean()
        data['Match Genome'].append(genome)
        data['ANI'].append(rbm_ani)
        ani[genome] = rbm_ani

    dfANI = pd.DataFrame(data)
    dfANI = dfANI.sort_values('ANI', ascending=False)
    gorder = dfANI['Match Genome'].tolist()

    return gorder, ani


def build_pos_line_plot(
        df, mgenome, outpre, cpos, rec, gorder, xmin, xmax, ymin, ymax, md
        ):
    ''' plots gene sequence identity on y axis vs genome coords on x axis.
        cpos is a list of contig lengths to add markers '''

    dfG = df[df['Match Genome'] != '-']
    dfG['Mid'] = dfG[['Start', 'Stop']].mean(axis=1)
    dfG = dfG[(dfG['Mid'] >= xmin) & dfG['Mid'] <= xmax]

    # write df to file
    group_df_file = f'{outpre}_group_data.tsv'
    dfG.to_csv(group_df_file, sep='\t', index=False)

    # Sort by ANI
    gorder, ani = ani_sort(dfG, gorder)

    colors = {
            'Conserved': '#e41a1c', # red high conserved core gene
            'Core': '#377eb8', # blue core gene 
            'Accessory': '#4daf4a', # green accessory gene
            'Specific': '#4daf4a', # Grouped specific with accessory genes
            #'Specific': '#984ea3', # purple genome specific gene
            }

    print(f'\n\tBuilding gene line plots ...')

    # set current genome length / 4, rounded up to whole integer
    glength = sum(mgenome.values())

    group_genomes = gorder #dfG['Match Genome'].unique()
    subs = len(group_genomes)

    # initialize the figure
    fig, axs = plt.subplots(subs, 1, figsize=(14, subs), sharex=True)

    for i, (ax, genome) in enumerate(zip(axs, group_genomes)):
        print(f'\t\tPlotting genome {genome}')
        dfS = dfG[dfG['Match Genome'] == genome]
        x = dfS['Mid'].to_list()
        # replace values < ymin with ymin
        ids = dfS['pID'].to_numpy()
        ymn = np.where(ids < ymin, ymin, ids)
        y = np.where(ids > ymax, ymax, ymn).tolist()
        c = [colors[i] for i in dfS['PanCat'].to_list()]

        # plot the data
        ax.scatter(x, y, color=c, marker='|', linestyle='-')

        # add genome ANI to plot
        rbm_ani = round(ani[genome], 2)
        #gname = f'genome: {genome} | ANI: {rbm_ani:.4f}'
        ax.text(
                1.04, 0.5, rbm_ani, ha='center', va='center',
                fontsize=14, transform = ax.transAxes
                )
        # add phylogroup and genomovar label to plot
        ax.text(# genomovar
                -0.08, 0.5, md[genome][1], ha='center', va='center',
                fontsize=14, transform = ax.transAxes
                )
        ax.text(# phylogroup
                -0.15, 0.5, md[genome][0], ha='center', va='center',
                fontsize=14, transform = ax.transAxes
                )

        # set axis limits
        ax.set_xlim(xmin-0.5, xmax+0.5)
        ax.set_ylim(ymin, ymax)
        ax.set_yticks([96, 98, 100])
        #ax.set_yticks([95, 97, 99], minor=True)

        # add dashed lines for recombinant threshold and ANI
        ax.hlines(y=rec, xmin=xmin, xmax=xmax, lw=1, colors='k', ls=':')
        ax.hlines(y=rbm_ani, xmin=xmin, xmax=xmax, lw=1, colors='k', ls='--')

        # set grid style
        ax.minorticks_on()
        ax.yaxis.grid(
            which="minor", color='#f0f0f0', linestyle='--', linewidth=.75
            )
        ax.yaxis.grid(
            which="major", color='#d9d9d9', linestyle='--', linewidth=1
            )
        ax.xaxis.grid(
            which="minor", color='#f0f0f0', linestyle='--', linewidth=.8
            )
        ax.xaxis.grid(
            which="major", color='#d9d9d9', linestyle='--', linewidth=1
            )
        ax.set_axisbelow(True)


    # add contig markers
    for mark in cpos[:-1]:
        if mark > xmin and mark < xmax:
            ax.text(mark, ymin,"l", color='#969696',
                    fontsize=8, ha='center', va='bottom')
    # add column labels
    axs[0].text(
            1.04, 1.1, 'ANI (%)', ha='center', va='center',
            fontsize=10, transform = axs[0].transAxes
            )
    # add phylogroup and genomovar label to plot
    axs[0].text(# genomovar
            -0.08, 1.1, 'Genomovar', ha='center', va='center',
            fontsize=10, transform = axs[0].transAxes, rotation=35
            )
    axs[0].text(# phylogroup
            -0.15, 1.1, 'Phylogroup', ha='center', va='center',
            fontsize=10, transform = axs[0].transAxes, rotation=35
            )
    # add x axis label
    axs[-1].set_xlabel('Genome position (mega base pairs)', fontsize=12)
    # y axis label
    axs[-1].text(# y axis label
            -0.02, -0.35, 'Sequence\nidentity (%)', ha='center', va='top',
            fontsize=10, transform = axs[-1].transAxes
            )

    # build the legend
    legend_labels = ['Conserved', 'Core', 'Accessory']
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

    #fig.set_tight_layout(True)
    fig.subplots_adjust(hspace=0.1, wspace=0, top=0.95, bottom=0.07, left=0.15, right=0.94)
    #fig.subplots_adjust(hspace=-0.03)
    plt.savefig(f'{outpre}_posline.pdf')
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
        '-md', '--metadata_file',
        help='Please specify the metadata file!',
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
    parser.add_argument(
        '-xmin', '--xaxis_minimum',
        help='(OPTIONAL) Specify the x-axis minimum (default = 0).',
        metavar='',
        type=float,
        default=0,
        required=False
        )
    parser.add_argument(
        '-xmax', '--xaxis_maximum',
        help='(OPTIONAL) Specify the x-axis maximum (default = 1000000).',
        metavar='',
        type=float,
        default=1000000,
        required=False
        )
    parser.add_argument(
        '-ymin', '--yaxis_minimum',
        help='(OPTIONAL) Specify the y-axis minimum (default = 95.0).',
        metavar='',
        type=float,
        default=95.0,
        required=False
        )
    parser.add_argument(
        '-ymax', '--yaxis_maximum',
        help='(OPTIONAL) Specify the y-axis maximum (default = 100.0).',
        metavar='',
        type=float,
        default=100,
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
    metadata = args['metadata_file']
    rec = args['recombination_cutoff']
    xmin = args['xaxis_minimum']
    xmax = args['xaxis_maximum']
    ymin = args['yaxis_minimum']
    ymax = args['yaxis_maximum']

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
    # parse the metadata file
    md = parse_metadata_file(metadata)

    ## SECTION 02: Create data frame
    # this step takes all the input we just parsed and creates a dataframe
    df, cpos = combine_input_data(mgenome, mgenes, RBM, pancats, repgenes, annos)

    ## SECTION 04: gene RBM identity vs. genome position
    # y-axis min and max happens during the plotting
    _ = build_pos_line_plot(
            df, mgenome, outpre, cpos, rec, gorder, xmin, xmax, ymin, ymax, md
            )

    print(f'\n\nComplete success space cadet!! Finished without errors.\n\n')

if __name__ == "__main__":
    main()


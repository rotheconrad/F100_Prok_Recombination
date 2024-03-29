#!/usr/bin/env python

''' rRBM rarefaction curve by clade

This script takes 3 input files.

The first file is the tab separated RBM matrix generated by the
03g_Recombinant_group_analysis.py script. Each column represents a
genome and each row represents an RBM gene. Genome IDs should be
on the first line (row). Each row after that is a RBM gene with a value
for each genome of 0 (absent), 1 (present), or 2 (conserved).

The second file is generated by the user and is a two column, comma-
separated list of genomeID,userAssignedClade where genomeID matches the
genome IDs from line 1 of the rbm matrix file and the userAssignedClade
can be anything the user assigns. Only genomes included in this file are
selected from the rbm matrix Example:

genome1,clade1
genome2,clade1
genome3,clade2
genome4,clade3
genome5,clade1
genome6,clade3
genome7,clade3
genome8,clade2

The third file is a two column, comma-separated list of clade,color.

clade1,#e41a1c
clade2,#377eb8
clade3,#4daf4a

Color guide: https://colorbrewer2.org/#type=qualitative&scheme=Set1&n=9

* conserved sites are masked (not included)

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: May 2023
License :: GNU GPLv3
Copyright 2023 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


def parse_clade_list(clade_list):
    ''' Builds and returns a dict of {genomeID, clade} from csv file '''
    
    clades = {}

    with open(clade_list, 'r') as file:
        for line in file:
            X = line.rstrip().split(',')
            genomeID = X[0]
            clade = X[1]
            clades[genomeID] = clade

    return clades


def parse_color_list(color_list):
    ''' Builds and returns dict of {cladeID, color} from csv file '''

    colors = {}

    with open(color_list, 'r') as file:
        for line in file:
            X = line.rstrip().split(',')
            cladeID = X[0]
            color = X[1]
            colors[cladeID] = color

    return colors


def parse_rbm_matrix(rbm_matrix, clades):
    ''' Creates a pandas dataframe from input matrix and labels clades.
        Removes conserved genes (labeled with 2) '''

    # read in the matrix as a pandas data frame
    df = pd.read_csv(rbm_matrix, sep='\t', header=0)
    # collect clade labels and gene number
    clades = [clades[i] for i in df.columns]
    genes = [f'Gene_{i+1:03}' for i in range(len(df))]
    # add gene number as index
    df['Genes'] = genes
    df = df.set_index('Genes')
    df.columns.names = ['Genomes']
    df = df[df[df.columns[0]] != 2] # remove conserved genes
    # transform df and add clades
    df = df.T
    df['Clade'] = clades

    return df


def compute_rarefaction(df, dfT, clade_perm, colors, sort_order):
    ''' iterates the pandas dataframe and stores rarefaction data '''

    data = {
            'Genome Count': [], 'Recombining Fraction': [], 'Color': [],
            'New Sites': [], 'Clade': []
            }
    
    genome_count = 0
    site_count = 0
    glist = []
    total_genes = len(dfT)

    for clade in clade_perm:
        color = colors[clade]
        dfx = df[df['Clade'] == clade].copy()
        dfx['Sum'] = dfx.sum(axis=1, numeric_only=True)
        dfx = dfx.sort_values('Sum',  ascending=sort_order)
        dfx = dfx.drop(columns=['Clade', 'Sum']).T
        xgenomes = dfx.columns.tolist()
        for g in xgenomes:
            genome_count += 1
            glist.append(g)
            x = dfT[glist].sum(axis=1, numeric_only=True).tolist()
            rbm_count = len([i for i in x if i != 0])
            # percent new sites
            #new_sites = round((rbm_count - site_count)/total_genes * 100, 2)
            # or count of new sites
            new_sites = (rbm_count - site_count)
            site_count = rbm_count
            recombining_fraction = round(rbm_count/total_genes * 100, 2)
            data['Genome Count'].append(genome_count)
            data['Recombining Fraction'].append(recombining_fraction)
            data['Color'].append(color)
            data['New Sites'].append(new_sites)
            data['Clade'].append(clade)

    df2 = pd.DataFrame(data)

    return df2


def build_rarefaction_plot(df, out_pre):
    ''' reads the dataframe and plots rarefaction curve by clade '''

    # write out the data
    df.to_csv(f'{out_pre}_data.tsv', sep='\t', index=False)

    fig, (ax1, ax2) = plt.subplots(
                                2, 1, figsize=(7.5,8), sharex=True,
                                gridspec_kw={'height_ratios': [2, 1]},
                                )
    fz = 12 # label font size

    ax1.set_title('', fontsize=fz)
    ax2.set_title('', fontsize=fz)
    ax1.set_ylabel('Recombining Fraction (%)', fontsize=fz)
    #ax2.set_ylabel('New Sites (%)', fontsize=fz) # percent new sites
    ax2.set_ylabel('New Sites', fontsize=fz) # or count of new sites
    ax2.set_xlabel('', fontsize=fz)
    ax2.set_xlabel('Genomes', fontsize=fz)

    legend_elements = []

    for clade, color in zip(df['Clade'].unique(), df['Color'].unique()):
        dfx = df[df['Clade'] == clade]
        x = dfx['Genome Count'].tolist()
        y1 = dfx['Recombining Fraction'].tolist()
        y2 = dfx['New Sites'].tolist()

        ax1.plot(x, y1, color=color, linestyle='-', linewidth=5)
        ax2.bar(x, y2, color=color, width=0.8)

        lg = Line2D(
        [0],[0], color='w', label=clade,
        markerfacecolor=color, marker='s', markersize=fz
        )
        legend_elements.append(lg)

    ax1.legend(
        title='Clade',
        handles=legend_elements,
        loc='lower right',
        fontsize=fz,
        fancybox=True,
        framealpha=0.0,
        frameon=False
        )

    # set ax1 y axis
    ax1.set_ylim(ymin=0, ymax=100)
    ax1.set_yticks(range(0, 101, 10))

    for ax in [ax1, ax2]:
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

    plt.subplots_adjust(
        left = 0.09,
        right = 0.98,
        bottom = 0.07,
        top = 0.98,
        hspace = 0.05
        )

    plt.savefig(f'{out_pre}_plot.pdf')
    plt.close()

    return True


def build_bar_plot(df, out_pre):

    return True


def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-r', '--rbm_matrix',
        help='Please specify the rbm_matrix.tsv input file!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-l', '--clade_list',
        help='Please specify the clade list input file!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-c', '--color_list',
        help='Please specify the color list input file!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-o', '--output_prefix',
        help='What do you want to name the output files?',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-sort_order', '--sort_order',
        help='(OPTIONAL): sort order. default: ascending=False',
        metavar='',
        type=bool,
        required=False,
        default=False
        )
    args=vars(parser.parse_args())

    # Run this scripts main function
    print('\n\nCalculating the rarefaction curve ...\n')

    # define input params
    rbm_matrix = args['rbm_matrix']
    clade_list = args['clade_list']
    color_list = args['color_list']
    out_pre = args['output_prefix']
    sort_order = args['sort_order']

    # parse the input files
    clades = parse_clade_list(clade_list)
    colors = parse_color_list(color_list)
    df = parse_rbm_matrix(rbm_matrix, clades)

    # build the plots
    # get the clade list
    clades = df.Clade.unique().tolist()
    # get the number of genomes in each clade
    sizes = [len(df[df['Clade'] == clade]) for clade in clades]
    # get the clade order based on clade size
    order = np.argsort(sizes)
    # need a transformed version of the df as well
    dfT = df.drop(columns=['Clade']).T # transformed dataframe
    for i in order[::-1]:
        clade_perm = clades[i:] + clades[:i]
        df2 = compute_rarefaction(df, dfT, clade_perm, colors, sort_order)
        _ = build_rarefaction_plot(df2, f'{out_pre}_{i+1:02}')

    print(f'\n\nComplete success space cadet!! Finished without errors.\n\n')

if __name__ == "__main__":
    main()
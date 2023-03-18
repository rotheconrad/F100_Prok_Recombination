#!/usr/bin/env python

''' Plot clustermap from binary matrix with header row and index column


-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: March 20th, 2020
License :: GNU GPLv3
Copyright 2020 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse, sys
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

sys.setrecursionlimit(50000)

def plot_clustermap(
                binary_matrix, outfile, core_threshold, w, h, x, excore, exspec
                ):
    ''' Takes a binary matrix input file in tsv format that includes
        a header row and index column and builds a clustermap plot '''

    # set the outfile name
    #outfile = binary_matrix.split('.')[0] + 'clustermap.pdf'
    # Read in the tsv binary matrix to a pandas dataframe with head and index
    df = pd.read_csv(binary_matrix, sep='\t', header=0)

    # set total or length of genomes (entries) in the set
    n = len(df.columns)
    # set value to consider core
    c = n * core_threshold
    # count number of core
    core = len(df[(df.sum(axis=1) >= c) & (df.sum(axis=1) <= n)])
    # count number of genome specific
    specific = len(df[df.sum(axis=1) == 0])
    # Compose data line to print at top of plot
    data_line = (
                f'Total Genomes: {n} | Core rRBM Genes: {core} | '
                f'Genes without rRBM: {specific} '
                )
    print(data_line, '\n\n')
    # set colors
    colors = ['#f0f0f0', '#000000', '#fee391']
    # plot it

    # Exclude core genes if ec flag
    if excore:
        df = df[df.sum(axis=1) < c]

    # Exlude genome-specific if es flag
    if exspec:
        df = df[df.sum(axis=1) > 1]

    g = sns.clustermap(
                    df, figsize=(w,h),
                    metric="euclidean", method="ward",
                    cmap=colors, vmin=0, vmax=2, 
                    cbar_kws={"ticks":[0,1,2]},
                    row_cluster=False
                    )
    # Retrieve ax object to access axes features
    ax = g.ax_heatmap
    ax.legend(labels=['Non-recombinant', 'Recombinant', 'Conserved'])
    '''
    # add data text to top of figure
    plt.text(
            0.5, 1.2, data_line,
            fontsize=f, color='#980043',
            horizontalalignment='center', verticalalignment='center',
            transform=ax.transAxes
            )
    '''
    # turn off y-axis labels
    ax.set_yticks([])
    ax.set_ylabel('Genes (in genome order)')
    # rotate x-axis labels
    plt.setp(ax.get_xticklabels(), rotation=90, fontsize=x)
    # turn off x-axis labels
    #ax.set_xticks([])
    # adjust plot margins
    #plt.subplots_adjust()
    plt.tight_layout()
    # save figure and close
    plt.savefig(outfile)
    plt.close

    return True

def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-b', '--binary_matrix_tsv_file',
        help='Please specify the binary matrix input tsv file!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-o', '--output_file',
        help='What do you want to name the output file? (use .pdf)',
        #metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-c', '--set_core_threshold',
        help='(Optional) Set the cutoff for core gene class (Default = 1.0).',
        metavar='',
        type=float,
        required=False,
        default=1.0
        )
    parser.add_argument(
        '-x', '--set_figure_width',
        help='(Optional) Set the figure width (Default = 24).',
        metavar='',
        type=int,
        required=False,
        default=24
        )
    parser.add_argument(
        '-y', '--set_figure_height',
        help='(Optional) Set the figure height (Default = 14).',
        metavar='',
        type=int,
        required=False,
        default=14
        )
    parser.add_argument(
        '-a', '--xaxis_text_size',
        help='(Optional) Set the size of x-axis text of figure (Default = 18).',
        metavar='',
        type=int,
        required=False,
        default=18
        )
    parser.add_argument(
        '-ec', '--exclude_core_genes',
        help='(Optional) Cluster and plot without the core genes.',
        action='store_true',
        required=False,
        )
    parser.add_argument(
        '-es', '--exclude_genome_specific_genes',
        help='(Optional) Cluster and plot without the genome-specific genes.',
        action='store_true',
        required=False,
        )
    args=vars(parser.parse_args())

    print(
        '\nBuilding pangenome clustermap with seaborn.clustermap using '
        'metric: euclidean and method: ward. For more info see: '
        'https://seaborn.pydata.org/generated/seaborn.clustermap.html'
        )

    _ = plot_clustermap(
                    args['binary_matrix_tsv_file'],
                    args['output_file'],
                    args['set_core_threshold'],
                    args['set_figure_width'],
                    args['set_figure_height'],
                    args['xaxis_text_size'],
                    args['exclude_core_genes'],
                    args['exclude_genome_specific_genes']
                    )

if __name__ == "__main__":
    main()
#!/usr/bin/env python

''' Plot clustermap from binary matrix with header row and index column

This tool takes the _rbm_matrix.tsv output from the script
03g_Recombinant_group_analysis.py and returns a clustered
heatmap in .pdf format.

This script requires the following packages:

    * argparse
    * pandas
    * matplotlib
    * seaborn

This script can plot meta data with the clustermap if provided.
Requires two files:
    1) A tab separated metadata file with a row for each genome in the
    ANI file with the exact genome name and columns for each meta value.
    File should include a header.
            example:
                Genome\tSite\tNiche\tPhylogroup
                Genome1\tA\tSheep\tB1
                Genome2\tB\tCow\tE

    2) A comma separated file of unique meta values,color
            example:
                Site1,#ffffb3
                Site2,#377eb8
                Species1,#ff7f00
                Species2,#f781bf
                Temp1,#4daf4a
                Temp2,#3f1gy5


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


def parse_colors(metacolors):

    ''' reads meta colors file into a dict of {meta value: color} '''

    cdict = {}

    with open(metacolors, 'r') as file:
        for line in file:
            X = line.rstrip().split(',')
            mval = X[0]
            color = X[1]
            cdict[mval] = color

    return cdict


def plot_clustermap(
            df, outfile, core_threshold, w, h, x, excore, exspec, metadf, cdict
            ):

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

    # build the plot with metadata
    if not metadf.empty:
        # Build legend
        for meta in metadf.columns:
            labels = metadf[meta].unique()

            fig, ax = plt.subplots(figsize=(10,10))
            ax.xaxis.set_visible(False)
            ax.yaxis.set_visible(False)

            for label in labels:
                ax.bar(
                    0, 0, color=cdict[label], label=label, linewidth=0
                    )

            ax.legend(
                title=meta, title_fontsize='xx-large', loc="center",
                frameon=False, markerscale=5, fontsize='xx-large'
                )
            outpre = '.'.join(outfile.split('.')[:-1])
            plt.savefig(f'{outpre}_Legend_{meta}.pdf', dpi=300)
            plt.close()

        # Build the clustermap
        metadf = metadf.replace(cdict) # swap colors in for values
        g = sns.clustermap(
                        df, figsize=(w,h),
                        metric="euclidean", method="ward",
                        cmap=colors, vmin=0, vmax=2, 
                        cbar_kws={"ticks":[0,1,2]},
                        row_cluster=False,
                        col_colors=metadf
                        )

    # build the plot without metadata
    else:
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
        '-m', '--meta_data_file',
        help='Please specify a meta data file!',
        metavar=':',
        type=str,
        required=False
        )
    parser.add_argument(
        '-c', '--meta_colors_file',
        help='Please specify a meta colors file!',
        metavar=':',
        type=str,
        required=False
        )
    parser.add_argument(
        '-o', '--output_file',
        help='What do you want to name the output file? (use .pdf)',
        #metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-sct', '--set_core_threshold',
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
        '-ecg', '--exclude_core_genes',
        help='(Optional) Cluster and plot without the core genes.',
        action='store_true',
        required=False,
        )
    parser.add_argument(
        '-esg', '--exclude_genome_specific_genes',
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

    # define input params
    binary_matrix = args['binary_matrix_tsv_file']
    metadata = args['meta_data_file']
    metacolors = args['meta_colors_file']
    outfile = args['output_file']
    core_threshold = args['set_core_threshold']
    w = args['set_figure_width']
    h = args['set_figure_height']
    x = args['xaxis_text_size']
    excore = args['exclude_core_genes']
    exspec = args['exclude_genome_specific_genes']

    # read in the binary matrix as pandas df
    df = pd.read_csv(binary_matrix, sep='\t', header=0)

    if metadata and metacolors:
        cdict = parse_colors(metacolors)
        metadf = pd.read_csv(metadata, sep='\t', index_col=0, header=0)
        metadf = metadf.reindex(df.T.index).dropna()
        # select only columns in metadf
        df = df[metadf.index.tolist()]
        _ = plot_clustermap(
            df, outfile, core_threshold, w, h, x, excore, exspec, metadf, cdict
                )

    elif metadata:
        print('\n\nBoth meta data and meta colors are required to use them!')

    elif metacolors:
        print('\n\nBoth meta data and meta colors are required to use them!')

    else:
        # create the plot without meta data colors
            _ = plot_clustermap(
            df, outfile, core_threshold, w, h, x, excore, exspec, None, None
                )

    print('\n\nComplete success space cadet! Hold on to your boots.\n\n')

if __name__ == "__main__":
    main()
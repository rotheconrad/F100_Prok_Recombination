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
    

def plot_clustermap(df, outfile, core, W, H, ts, excore, exspec, metadf, cdict):
    ''' Takes a binary matrix input file in tsv format that includes
        a header row and index column and builds a clustermap plot '''

    # set total or length of genomes (entries) in the set
    n = len(df.columns)
    # set value to consider core
    c = n * core
    # count number of core
    core = df[df.sum(axis=1) >= c].shape[0]
    # count number of genome specific
    specific = df[df.sum(axis=1) == 1].shape[0]
    # count number of variable
    variable = df[(df.sum(axis=1) > 1) & (df.sum(axis=1) < c)].shape[0]
    # Compose data line to print at top of plot
    data_line = (
                f'Total Genomes: {n} | Core Genes: {core} | '
                f'Genome Specific Genes: {specific} | '
                f'Variable Genes: {variable}'
                )
    print('\n\nPangenome Summary:\n')
    print(data_line, '\n\n')

    # set colors
    colors = ['#f0f0f0', '#525252']

    # Exclude core genes if ec flag
    if excore:
        df = df[df.sum(axis=1) < c]

    # Exlude genome-specific if es flag
    if exspec:
        df = df[df.sum(axis=1) > 1]

    # build the plot with metadata
    if not metadf.empty:
        print('Building plot with metadata ...')
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
                        df, figsize=(W, H),
                        metric="euclidean", method="ward",
                        cmap=colors, vmin=0, vmax=1, 
                        cbar_kws={"ticks":[0,1]},
                        col_colors=metadf
                        )

    # build the plot without metadata
    else:
        print('Building plot without metadata ...')
        g = sns.clustermap(
                        df, figsize=(W, H),
                        metric="euclidean", method="ward",
                        cmap=colors, vmin=0, vmax=1, 
                        cbar_kws={"ticks":[0,1]}
                        )

    # Retrieve ax object to access axes features
    ax = g.ax_heatmap

    # turn off y-axis labels
    ax.set_yticks([])
    # rotate x-axis labels
    plt.setp(ax.get_xticklabels(), rotation=90, fontsize=ts)
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
        help='What do you want to name the output file? (Use .pdf)',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-c', '--set_core_threshold',
        help='(Optional) Set the cutoff for core gene class (Default = 0.9).',
        metavar='',
        type=float,
        required=False,
        default=0.9
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
    parser.add_argument(
        '-md', '--meta_data_file',
        help='Please specify a meta data file!',
        metavar=':',
        type=str,
        required=False
        )
    parser.add_argument(
        '-mc', '--meta_colors_file',
        help='Please specify a meta colors file!',
        metavar=':',
        type=str,
        required=False
        )
    args=vars(parser.parse_args())

    # define parameters
    binary_matrix = args['binary_matrix_tsv_file']
    outfile = args['output_file']
    core = args['set_core_threshold']
    W = args['set_figure_width']
    H = args['set_figure_height']
    ts = args['xaxis_text_size']
    excore = args['exclude_core_genes']
    exspec = args['exclude_genome_specific_genes']
    metadata = args['meta_data_file']
    metacolors = args['meta_colors_file']

    print(
        '\nBuilding pangenome clustermap with seaborn.clustermap using '
        'metric: euclidean and method: ward. For more info see: '
        'https://seaborn.pydata.org/generated/seaborn.clustermap.html'
        )

    # Read in the tsv binary matrix to a pandas dataframe with head and index
    print('\n\n', binary_matrix, '\n\n')
    df = pd.read_csv(binary_matrix, sep='\t', header=0, index_col=0)

    if metadata and metacolors:
        print('\n\nParsing metadata and metacolors ...')
        cdict = parse_colors(metacolors)
        metadf = pd.read_csv(metadata, sep='\t', index_col=0, header=0)
        print('\n\nBuilding the plot ...')
        _ = plot_clustermap(
                df, outfile, core, W, H, ts, excore, exspec, metadf, cdict
                )

    elif metadata:
        print('\n\nBoth meta data and meta colors are required to use them!')

    elif metacolors:
        print('\n\nBoth meta data and meta colors are required to use them!')

    else:
        print('\n\nBuilding the plot ...')
        # create the plot without meta data colors
        # use empty dataframe to skip meta data and colors
        E = pd.DataFrame()
        _ = plot_clustermap(df, outfile, core, W, H, ts, excore, exspec, E, E)

    print('\n\nComplete success space cadet! Hold on to your boots.\n\n')

if __name__ == "__main__":
    main()
#!/usr/bin/env python

'''Plots clustered heatmap of All vs. All ANI matrix.

This tool takes the fastANI_allV.ani output from Part 01, Step 02
as input, creates a distance matrix using 100-ANI and returns a
clustered heatmap in .pdf format

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
Date Created :: Feb 2023
License :: GNU GPLv3
Copyright 2023 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse
import pandas as pd
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
import seaborn as sns


def parse_ANI_file(infile):

    ''' creates a square distance matrix of 100 - ANI '''

    dist = {} # keep ANI distance of each genome pair
    gdict = {} # keep set of genome names
    data = {} # build dict to convert to dataframe (distance matrix)

    # read the file, compute distance 100 - ANI, store data in dist dict
    with open(infile, 'r') as file:
        # header = file.readline() # ANI output should not have a header
        for line in file:
            X = line.rstrip().split('\t')
            g1 = X[0].split('/')[-1].split('.')[0]
            g2 = X[1].split('/')[-1].split('.')[0]
            ani = float(X[2])
            d = 100 - ani
            dist['-'.join(sorted([g1, g2]))] = d
            # add all genome names to gdict to get the complete set of names
            gdict[g1] = ''
            gdict[g2] = ''

    # get the sorted genome name list
    glist = sorted(list(gdict.keys()))
    #data['rows'] = glist # add rows or index names to dict.
    # iterate through genome list and dictionary to create square matrix df
    for g1 in glist:
        data[g1] = []
        for g2 in glist:
            key = '-'.join(sorted([g1, g2]))
            d = dist.get(key, 0)
            data[g1].append(d)

    df = pd.DataFrame(data, index=glist, columns=glist)

    return df


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


def plot_clustred_heatmap(df, outfile, metadf, cdict):

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
                        df, figsize=(16,9),
                        xticklabels=True,
                        yticklabels=True,
                        col_colors=metadf
                        )


    # build the plot without metadata
    else:
        g = sns.clustermap(
                        df, figsize=(16,9),
                        xticklabels=True,
                        yticklabels=True,
                        )

    g.savefig(outfile)
    plt.close()

    return True


def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-i', '--input_file',
        help='Please specify the input file name!',
        metavar=':',
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
        help='Please specify the output file name!',
        metavar=':',
        type=str,
        required=True
        )
    args=vars(parser.parse_args())

    # Do what you came here to do:
    print('\n\nRunning Script...')
    
    # define parameters
    infile = args['input_file']
    metadata = args['meta_data_file']
    metacolors = args['meta_colors_file']
    outfile = args['output_file']

    # read in all vs. all fastANI file.
    df = parse_ANI_file(infile)
    
    if metadata and metacolors:
        cdict = parse_colors(metacolors)
        metadf = pd.read_csv(metadata, sep='\t', index_col=0, header=0)
        metadf = metadf.reindex(df.index).dropna()
        # select only rows in metadf
        df = df[df.index.isin(metadf.index)]
        # select only columns in metadf
        df = df[metadf.index.tolist()]
        _ = plot_clustred_heatmap(df, outfile, metadf, cdict)

    elif metadata:
        print('\n\nBoth meta data and meta colors are required to use them!')

    elif metacolors:
        print('\n\nBoth meta data and meta colors are required to use them!')

    else:
        # create the plot without meta data colors
        _ = plot_clustred_heatmap(df, outfile, pd.DataFrame(), pd.DataFrame())

    print('\n\nComplete success space cadet! Hold on to your boots.\n\n')


if __name__ == "__main__":
    main()

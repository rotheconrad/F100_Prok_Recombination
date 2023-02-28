#!/usr/bin/env python

'''Plots clustered heatmap of All vs. All F100 matrix.

This tool takes the ${my_species}_F100.tsv output from Part 02, Step 03
as input, creates a distance matrix using 1-F100 and returns a clustered
heatmap in .pdf format

This script requires the following packages:

    * argparse
    * pandas
    * matplotlib
    * seaborn

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


def parse_F100_file(infile):

    ''' creates a square distance matrix of 1 - F100 '''

    dist = {} # keep F100 distance of each genome pair
    gdict = {} # keep set of genome names
    data = {} # build dict to convert to dataframe (distance matrix)

    # read the file, compute distance 1 - F100, store data in dist dict
    with open(infile, 'r') as file:
        header = file.readline()
        for line in file:
            X = line.rstrip().split('\t')
            g1 = X[0]
            g2 = X[1]
            F1 = float(X[3])
            d = 1 - F1
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


def plot_clustred_heatmap(df, outfile):

    # build the plot
    g = sns.clustermap(
                    df, figsize=(16,9),
                    xticklabels=True,
                    yticklabels=True
                    )

    # adjust layout, save, and close
    #plt.gca().invert_yaxis()
    #fig.set_tight_layout(True)
    g.savefig(outfile)
    plt.close()


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
    outfile = args['output_file']

    # read in the *.pim file from EMBL-EBI simple phylogeny
    df = parse_F100_file(infile)

    # create the plot
    _ = plot_clustred_heatmap(df, outfile)

    print('\n\nComplete success space cadet! Hold on to your boots.\n\n')


if __name__ == "__main__":
    main()

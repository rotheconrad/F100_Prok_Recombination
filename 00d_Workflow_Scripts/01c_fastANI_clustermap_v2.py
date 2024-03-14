#!/usr/bin/env python

'''Predict subpopulation ANI clusters and build heatmap with metadata.

This tool takes the fastANI_allV.ani output from Part 01, Step 02
as input. It has a few options and ultimately returns an ANI matrix
as a hierarchical clustered heatmap.

Their are three initial options:

1) Predict ANI clusters for local minimumns in the ANI distribution (default).
2) Predict ANI clusters for user defined values (-user).
    e.g. -user 95 97.5 98.5 99.5
3) Skip cluster prediction and use user defined metadata.
    e.g. -m meta_data.tsv -c meta_colors.tsv

The end result is a publication ready figure of an ANI clustermap with 
ANI cluster and/or user defined metadata in vector PDF format.

It will also write out a tsv file with genome name and predicted clusters
if options 1 or 2 were used.

This script requires the following packages:

    * argparse
    * pandas
    * matplotlib
    * seaborn

For option 3, the user defined metadata requires the following two files:

    1) A tab separated metadata file with a row for each genome in the
    ANI file with the exact genome name, and columns for each meta value.
    File should include a header.
            example:
                Genome\tSite\tNiche\tPhylogroup
                Genome1\tA\tSheep\tB1
                Genome2\tB\tCow\tE

    2) A tab separated file of unique meta catagories and colors
            example:
                Site1\t#ffffb3
                Site2\t#377eb8
                Species1\t#ff7f00
                Species2\t#f781bf
                Temp1\t#4daf4a
                Temp2\t#3f1gy5

* NOTE: The color ranges for the heatmap accommodate to 90% ANI difference
        The code will need modification for more diverse datasets.

To set different cluster method choose from:
    single, complete, average, weighted, centroid, median, ward
    *See webpage/docs for scipy.cluster.hierarchy.linkage for more info

To set different cluster metrics choose from:
    euclidean, seuclidean, correlation, hamming, jaccard, braycurtis
    *See webpage/docs for scipy.spatial.distance.pdist for additional 
     options and more info

-------------------------------------------
Author :: Roth Conrad w/ thanks to Dorian Feistel
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: Feb 2023
License :: GNU GPLv3
Copyright 2023 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import gaussian_kde
from scipy.signal import find_peaks
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.cluster.hierarchy import fclusterdata
from scipy.spatial.distance import pdist, squareform

def parse_ANI_file(infile, ani_min, ani_max):

    ''' creates a square distance matrix of 100 - ANI '''

    dist = {} # keep ANI distance of each genome pair
    gdict = {} # keep set of genome names
    data = {} # build dict to convert to dataframe (distance matrix)
    ani_array = []

    # read the file, compute distance 100 - ANI, store data in dist dict
    with open(infile, 'r') as file:
        # header = file.readline() # ANI output should not have a header
        for line in file:
            X = line.rstrip().split('\t')
            # genomes in pair
            g1 = X[0].split('/')[-1].split('.')[0]
            g2 = X[1].split('/')[-1].split('.')[0]
            # get ANI
            ani = round(float(X[2]), 2)
            # get second genome length (total kmers)
            g2size = int(X[4])
            # filter min and max ANI values
            if ani < ani_min or ani > ani_max: continue
            # add ani to array
            ani_array.append(ani)
            # convert to ANI distance
            #d = 100 - ani
            d = ani
            # store distance by genome pair name
            gpair = '-'.join(sorted([g1, g2]))
            if gpair in dist:
                gsize = dist[gpair][1]
                if g2size > gsize:
                    dist[gpair] = [d, g2size]
            else:
                dist[gpair] = [d, g2size]
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
            d = dist.get(key, [ani_min,0])[0]
            data[g1].append(d)

    df = pd.DataFrame(data, index=glist, columns=glist)

    return df, ani_array


def get_valleys(ani_array, outpre):

    # use peak_finder algorithm to find local minimums (valleys)

    kde = gaussian_kde(ani_array, bw_method=0.15)
    x = np.linspace(min(ani_array), max(ani_array), 1000)
    y = kde(x)
    yp = find_peaks(y) # returns array (peak indices in x, {properties dict})
    yn = y * -1 # inverse of y to find the valleys
    yv = find_peaks(yn) # these are the valley position
    pks = yp[0] # peaks
    vls = yv[0] # valleys
    local_maximums = [round(i,1) for i in sorted(x[pks], reverse=True)]
    local_minimums = [round(i,1) for i in sorted(x[vls], reverse=True)]

    max_txt = 'Local maximums: ' + ', '.join([str(i) for i in local_maximums])
    min_txt = 'Local minimums: ' + ', '.join([str(i) for i in local_minimums])
    line = f'{max_txt}\n{min_txt}'

    # Plotting
    plt.figure(figsize=(10, 6))
    plt.plot(x, y, label='KDE')
    plt.scatter(x[pks], y[pks], color='red', label='Peaks')
    plt.scatter(x[vls], y[vls], color='blue', label='Valleys', marker='v')
    plt.title('ANI KDE with Peaks and Valleys')
    plt.xlabel('ANI')
    plt.ylabel('Density')
    ax = plt.gca()
    plt.text(0.01, 0.99, line, horizontalalignment='left',
     verticalalignment='top', transform=ax.transAxes)
    plt.legend()
    
    plt.savefig(f'{outpre}_Local_minimums.pdf')

    return local_minimums


def predict_clusters(df, ani_thresholds, metric, method):

    # initialize dicts for meta data and colors
    metadf = {}
    cdict = {}

    genomes = df.columns.tolist()
    # convert ANI and thresholds to ANI distance
    df = (100-df)/100
    sorted_threshold = sorted(ani_thresholds)
    distance_thresholds = [round((100 - i)/100, 4) for i in sorted_threshold]

    # convert to condensed matrix
    # for some reason the condensed matrix gives (triangle) gives different
    # results than the uncondensed matrix (square)
    # squareform converts between condensed and uncondensed
    cond_matrix = squareform(df.to_numpy())

    for n, d in enumerate(distance_thresholds):
        label = f'clade-{n} ({sorted_threshold[n]}% ANI; t={d})'
        
        Z = linkage(cond_matrix, metric=metric, method=method)
        clusters = fcluster(Z=Z, t=d, criterion='distance')
        #clusters = fclusterdata(
                            #X=df, t=d, criterion='distance',
                            #metric=metric, method=method
                            #)
        named = [f'c{n:01}-{i:03}' for i in clusters]
        metadf[label] = named

    # prepare metadf as dataframe
    metadf = pd.DataFrame(metadf, index=genomes)

    # populate cdict with color for each unique metadf value
    count = 0
    for m in metadf.columns:
        count += 1
        unique_clusters = sorted(metadf[m].unique())
        colors = sns.color_palette("Paired", len(unique_clusters)).as_hex()
        if count % 2 == 0: colors.reverse()
        for i,u in enumerate(unique_clusters):
            cdict[u] = colors[i]

    return metadf, cdict


def parse_colors(metacolors):

    ''' reads meta colors file into a dict of {meta value: color} '''

    cdict = {}

    with open(metacolors, 'r') as file:
        for line in file:
            X = line.rstrip().split('\t')
            mval = X[0]
            color = X[1]
            cdict[mval] = color

    return cdict


def custom_cmap(ani_min, ani_max):

    '''
    # two colors per 1% ANI (e.g. light to dark)
    ani_range = int(np.ceil(ani_max) - np.floor(ani_min) + 1) * 2

    # Read these as colors per whole ANI value
    # eg: [100-99 lt-dark red; 99-98 lt-dk blue; 98-97 lt-dk green]
    # eg: [97-96 orange-yellow; 96-95 lt-dk purple; 95-94 lt-dk teal]
    # eg: [94-93 lt-dk pink; 93-92 lt-dk brown; 92-91 gray-black]

    colors = [
              '#67000d', '#fc9272', '#08306b', '#6baed6', '#00441b', '#a1d99b',
              '#ec7014', '#ffeda0', '#3f007d', '#9e9ac8',  '#01665e', '#c7eae5',
              '#c51b7d', '#fde0ef', '#8c510a', '#f6e8c3', '#d9d9d9', '#252525'
              ]

    cmap = sns.blend_palette(reversed(colors[:ani_range-1]), as_cmap=True)
    '''

    # one color per 1% ANI range
    ani_range = int(np.ceil(ani_max) - np.floor(ani_min) + 1)

    # Read these as colors per whole ANI value
    # eg: [100 red, 99 orange, 98 yellow, 97 green, ]

    colors = [
              '#e41a1c', '#ff7f00', '#ffff33', '#006d2c', '#80cdc1',
              '#377eb8', '#e78ac3', '#984ea3', '#bf812d', '#bababa'
              ]

    cmap = sns.blend_palette(reversed(colors[:ani_range]), as_cmap=True)

    return cmap


def plot_clustred_heatmap(
            df, outpre, cmap, metadf, cdict, ani_min, ani_max, metric, method
            ):

    # build the plot with metadata
    if not metadf.empty:
        # Build legend
        for meta in metadf.columns:
            labels = metadf[meta].unique()

            # create columns for larger metadata sets
            cols = int(np.ceil(len(labels)/20))

            fig, ax = plt.subplots(figsize=(cols*4,10))
            ax.xaxis.set_visible(False)
            ax.yaxis.set_visible(False)

            for label in labels:
                ax.bar(
                    0, 0, color=cdict[label], label=label, linewidth=0
                    )

            ax.legend(
                title=meta, title_fontsize='xx-large', loc="center",
                frameon=False, markerscale=5, fontsize='xx-large', ncol=cols
                )
            plt.savefig(f'{outpre}_Legend_{meta}.pdf')
            plt.close()

        # Build the clustermap
        colordf = metadf.replace(cdict) # swap colors in for values
        g = sns.clustermap(
                        df, figsize=(16,9), cmap=cmap,
                        xticklabels=True,
                        yticklabels=True,
                        col_colors=colordf,
                        vmin=ani_min, vmax=ani_max,
                        metric=metric, method=method
                        )


    # build the plot without metadata
    else:
        g = sns.clustermap(
                        df, figsize=(16,9), cmap=cmap,
                        xticklabels=True,
                        yticklabels=True,
                        vmin=ani_min, vmax=ani_max,
                        metric=metric, method=method
                        )

    # reorder the distance matrix and write to file
    orows = df.columns[g.dendrogram_row.reordered_ind] # row order
    ocols = df.index[g.dendrogram_col.reordered_ind] # col order
    df = df.reindex(orows)
    df = df[ocols]
    df.to_csv(f'{outpre}_distmat.tsv', sep='\t')

    metadf = metadf.reindex(orows)
    metadf.to_csv(f'{outpre}_predicted_clusters.tsv', sep='\t')

    with open(f'{outpre}_meta_colors.tsv', 'w') as cout:
        for meta, color in cdict.items():
            cout.write(f'{meta}\t{color}\n')

    g.cax.set_xlabel('ANI (%)', rotation=0, labelpad=10, fontsize=12)
    # remove labels and tick marks
    g.ax_heatmap.set_xlabel('')
    g.ax_heatmap.set_ylabel('')
    g.ax_heatmap.set_xticklabels([])
    g.ax_heatmap.set_yticklabels([])
    # If you want to also remove the ticks
    g.ax_heatmap.tick_params(axis='x', which='both', length=0)
    g.ax_heatmap.tick_params(axis='y', which='both', length=0)
    # write file
    g.savefig(f'{outpre}_clustermap.pdf')
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
        '-o', '--output_file_prefix',
        help='Please specify the prefix to use for output file names!',
        metavar=':',
        type=str,
        required=True
        )
    parser.add_argument(
        '-user', '--user_clusters',
        help='Input ANI value list for cluster predction (eg 95, 97.5, 99.5)!',
        metavar=':',
        type=float,
        nargs='+',
        required=False,
        default=None
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
        '-no', '--no_meta_data',
        help='Set -no True to build basic seaborn clustermap of ANI distance!',
        metavar=':',
        type=str,
        required=False,
        default=None
        )
    parser.add_argument(
        '-min', '--minimum_ANI',
        help='Minimum ANI to include. Removes data below. (Default: 95.0)!',
        metavar=':',
        type=float,
        required=False,
        default=95.0
        )
    parser.add_argument(
        '-max', '--maximum_ANI',
        help='Maximum ANI to include. Removes data above. (Default: 100.0)!',
        metavar=':',
        type=float,
        required=False,
        default=100.0
        )
    parser.add_argument(
        '-metric', '--cluster_metric',
        help='Cluster metric (default: euclidean)',
        metavar=':',
        type=str,
        required=False,
        default='euclidean'
        )
    parser.add_argument(
        '-method', '--cluster_method',
        help='Cluster metric (default: average)',
        metavar=':',
        type=str,
        required=False,
        default='average'
        )
    args=vars(parser.parse_args())

    # Do what you came here to do:
    print('\n\nRunning Script...')
    
    # define parameters
    infile = args['input_file']
    outpre = args['output_file_prefix']
    user_clusters = args['user_clusters']
    metadata = args['meta_data_file']
    metacolors = args['meta_colors_file']
    no_meta = args['no_meta_data']
    ani_min = args['minimum_ANI']
    ani_max = args['maximum_ANI']
    metric = args['cluster_metric']
    method = args['cluster_method']

    # read in all vs. all fastANI file.
    df, ani_array = parse_ANI_file(infile, ani_min, ani_max)

    # create custom color map for the heatmap gradient
    cmap = custom_cmap(ani_min, ani_max)
    
    if no_meta:
        # create the plot without meta data colors
        metadf, cdict = pd.DataFrame(), pd.DataFrame()
        _ = plot_clustred_heatmap(
            df, outpre, cmap, metadf, cdict, ani_min, ani_max, metric, method
            )

    elif metadata and metacolors:
        cdict = parse_colors(metacolors)
        metadf = pd.read_csv(metadata, sep='\t', index_col=0, header=0)
        metadf = metadf.reindex(df.index).dropna()
        # select only rows in metadf
        df = df[df.index.isin(metadf.index)]
        # select only columns in metadf
        df = df[metadf.index.tolist()]
        _ = plot_clustred_heatmap(
            df, outpre, cmap, metadf, cdict, ani_min, ani_max, metric, method
            )

    elif metadata:
        print('\n\nBoth meta data and meta colors are required to use them!')

    elif metacolors:
        print('\n\nBoth meta data and meta colors are required to use them!')

    elif user_clusters:
        # run clustering algo on user defined values
        metadf, cdict = predict_clusters(df, user_clusters, metric, method)

        _ = plot_clustred_heatmap(
            df, outpre, cmap, metadf, cdict, ani_min, ani_max, metric, method
            )

    else:
        # default case. find_peaks, get local minimums, predict clusters, plot
        predicted_clusters = get_valleys(ani_array, outpre)
        metadf, cdict = predict_clusters(df, predicted_clusters, metric, method)

        _ = plot_clustred_heatmap(
            df, outpre, cmap, metadf, cdict, ani_min, ani_max, metric, method
            )

    print('\n\nComplete success space cadet! Hold on to your boots.\n\n')


if __name__ == "__main__":
    main()

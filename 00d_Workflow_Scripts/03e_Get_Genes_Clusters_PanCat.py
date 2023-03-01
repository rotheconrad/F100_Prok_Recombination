#!/usr/bin/env python

'''Generate .tsv file Genes Clusters Pangenome Categories

Takes the Binary Matrix .tsv file from 03 and the Mmseqs2 cluster file 
as input and outputs a .tsv file of:
Gene_Name Cluster_Name Pangenome_Category
Where Pangenome_Category is one of: (Core, Accessory, Specific)

This tool takes the following input parameters:

    * bnry - binary gene presence absence matrix .tsv file
    * clstr - MMSeqs2 tsva cluster file
    * out - output file name

This script returns the following files:

    * outputfilename.tsv

This script requires the following packages:

    * argparse

This file can also be imported as a module and contains the follwing 
functions:

    * get_category - reads the bnry file and assigns pangenome category
    * build_the_list - Orchestrates Parsing, computing, and output
    * main - the main function of the script

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: Thursday, September 5th, 2019
License :: GNU GPLv3
Copyright 2019 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse
from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt

def get_category(bnry):
    ''' Read bnry file and calculate the category of each cluster '''

    pan_category = {} # by cluster. Core, Accessory, Specific 

    with open(bnry, 'r') as f:
        header = f.readline().rstrip().split()
        genomes = header # genome names are the header
        n = len(genomes) # This is the number of genomes ie(100)
        # Define gene category cutoffs
        core, specific = 0.90, 1/n
        # each cluster is a line in the bnry file
        # compute gene category for each cluster
        for l in f:
            X = l.rstrip().split('\t')
            cluster = X[0] # Cluster # is the first column of bnry
            gs = [int(i) for i in X[1:]] # genes in cluster count
            g = sum(gs) # This is how many genomes cluster is in
            v = g / n # percent genomes with gene
            if v == specific: pan_category[cluster] = ['Specific', v]
            elif v < core: pan_category[cluster] = ['Accessory', v]
            elif v >= core: pan_category[cluster] = ['Core', v]
            else: print(l, v, 'error in assigning gene category')

    return pan_category


def parse_cluster_data(clstr, pan_category):

    cluster_data = defaultdict(lambda: defaultdict())
    data = {}

    # read the cluster and alignment data
    with open(clstr, 'r') as f:
        for l in f:
            X = l.rstrip().split('\t')
            cluster = X[0]
            gene = X[1]
            genome = '_'.join(gene.split('_')[:-1])
            pcat = pan_category[cluster] # [pancat, n/N]
            alnid = float(X[3])
            dist = 1 - alnid

            # look for highly conserved genes
            # dist = 1 - sequence alignmend id
            # store only one gene per genome for each cluster.
            # * gene with the smallest distance to cluster rep.
            if pcat[0] == 'Core':
                if genome in cluster_data[cluster]:
                    test_dist = cluster_data[cluster][genome][1]
                    if dist < test_dist:
                        cluster_data[cluster][genome] = [gene,dist]
                else:
                    cluster_data[cluster][genome] = [gene, dist]

            # store gene data to for outfile
            data[gene] = [gene, cluster, pcat[0], pcat[1], dist]

    return cluster_data, data


def get_avg_dists(cluster_data):

    cluster_dist = {'Cluster': [], 'Avg_dist': []}
    cluster_genes = defaultdict(list)

    # create cluster distance distributions.
    # each core cluster has one gene from each genome at this point.
    # core genes are 90% - 100% of the genomes.
    # average distance of genes in cluster = sum(dist) / len(dist)
    for cluster, genomes in cluster_data.items():
        distances = []
        for genome, values in genomes.items():
            gene = values[0]
            dist = values[1]
            distances.append(dist)
            cluster_genes[cluster].append(gene)
        cluster_dist['Cluster'].append(cluster)
        cluster_dist['Avg_dist'].append(sum(distances)/len(distances))

    return cluster_dist, cluster_genes


def plot_dist(input_array, out):

    plot_out = out.split('.')[0] + '.pdf'
    qarray = sorted(input_array)
    q05 = np.quantile(qarray, 0.05)

    fig, ax = plt.subplots(figsize=(7, 5))

    ax.hist(input_array, bins=75, edgecolor='k', color='#08519c', alpha=0.6,)
    ax.axvline(q05, color='#ae017e', linestyle='dashed', linewidth=1)
    ax.set_xlabel('Average sequence distance of cluster', fontsize=14)
    ax.set_ylabel('Count', fontsize=14)

    #ax.set_yscale('log')
    ax.minorticks_on()
    ax.tick_params(
        which='minor', axis='both', left=False, bottom=False
        )
    ax.tick_params(
                which='major', axis='both',
                left=True, bottom=True,
                size=6, width=2, tickdir='inout',
                labelsize=12, zorder=10
                )
    ax.yaxis.grid(
        which="major", color='#bdbdbd', linestyle='--',
        linewidth=1, alpha=0.4, zorder=1
        )
    ax.set_axisbelow(True)
    for spine in ax.spines.values(): spine.set_linewidth(2)

    fig.set_tight_layout(True)
    plt.savefig(plot_out)
    plt.close()

    return q05


def build_the_list(bnry, clstr, out):
    ''' Reads in the files, builds the list, and writes to out '''

    # read in the binary matrix file and find Pan Cat for each cluster
    pan_category = get_category(bnry)
    # read in the cluster file and get distances for each cluster
    cluster_data, data = parse_cluster_data(clstr, pan_category)
    # get average distance for each cluster and cluster, gene names
    cluster_dist, cluster_genes = get_avg_dists(cluster_data)
    # select the lowest 5% of clusters as highly conserved.
    q05 = plot_dist(cluster_dist['Avg_dist'], out)
    tmp_data = zip(cluster_dist['Cluster'], cluster_dist['Avg_dist'])
    for clust, dist in tmp_data:
        if dist <= q05:
            for gene in cluster_genes[clust]:
                data[gene][2] = 'Conserved' 
    # write out the data
    with open(out, 'w') as o:
        o.write('Gene_Name\tCluster_Name\tPangenome_Category\tn/N\tDistance\n')
        # n/N = number of genomes with gene in cluster / total genomes
        for gene, v in data.items():
            line_out = f'{v[0]}\t{v[1]}\t{v[2]}\t{v[3]}\t{v[4]}\n'
            o.write(line_out)

    return True


def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-b', '--binary_file',
        help='Please specify the binary.tsv input file!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-c', '--cluster_file',
        help='Please specify the Mmseqs2 Cluster TSV File!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-o', '--output_file',
        help='Please specify the name to us for the output .tsv file',
        metavar='',
        type=str,
        required=True
        )
    args=vars(parser.parse_args())

    # Run this scripts main function
    print('Generating Genes, Clusters, and PanCat list...')

    build_the_list(
        args['binary_file'],
        args['cluster_file'],
        args['output_file']
        )

if __name__ == "__main__":
    main()

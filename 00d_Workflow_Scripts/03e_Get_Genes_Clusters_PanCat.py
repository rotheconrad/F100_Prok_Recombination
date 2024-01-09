#!/usr/bin/env python

'''Generate .tsv file Genes Clusters Pangenome Categories

Takes the Binary Matrix .tsv file and the Mmseqs2 cluster file 
as input and outputs a .tsv file with columns:
Gene_Name, Cluster_Name, Pangenome_Category, n/N, Distance

* Gene name is the name of all individual genes in the pangenome
* Cluster_Name is the representative gene of the assigned cluster
* Pangenome_Category is one of: (Conserved, Core, Accessory, Specific)
* n/N is number of genomes with gene category over total genomes in analysis
* Distance is the average sequence distance from the representative gene.
* Conserved genes are a subset of the least divergent core genes.
* Core genes are found in ≥ 90% of genomes in the analysis.
* Specific genes are found in only 1 genome in the analysis.
* Accessory genes are all other genes.

This tool takes the following input parameters:

    * bmat - binary gene presence absence matrix .tsv file
    * mmsq - MMSeqs2 tsva cluster file
    * rbms - All vs all RBM file
    * out - output file name

(OPTIONAL) parameters for defining the most conserved gene clusters.
The average within cluster sequence distance is computed for each cluster.
The Core gene clusters with least sequence distance are labeled "Conserved"
Gene clusters with distance = 0 are labeled as conserved and removed from
the distribution. Lower quantile is computed on remaining distribution.

    * c - Fraction of genomes required to have a gene for that gene to be 
        ... classified as a core gene.

    * qt - Lower quantile to define additional "Conserved" gene clusters.
        ... this quantile is computed from the within cluster average
        ... sequence distance distribution after 0 dist genes are removed.
        ... default = 0.1; recommend values 0 - 0.5.

    * xmax - Use this to reduce the max x-axis value plotted to zoom into
        ... the histogram and get a closer look. Typically not necessary
        ... within species, but useful for closely related 90-95% ANI.

This script returns the following files:

    * outputfilename.tsv

This script requires the following packages:

    * argparse

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


def parse_rbm_file(rbms):
    ''' Read all vs all RBMs into dict.
    RBMs = {'-'.join(sort(gA,gB)): 100-pID}

    '''
    print('\n\tReading RBM file ...')
    
    # initialize dictionary to store rbm data
    # {gene-names: pID} * pID is 100 - sequence identity
    rbm_dict = {}

    with open(rbms, 'r') as file:
        for line in file:
            X = line.rstrip().split('\t')
            geneA = X[0]
            geneB = X[1]
            # sort the pair name for easy lookup later
            match = '-'.join(sorted([geneA, geneB]))
            # take the sequence identity as distance 100-pID
            pID = 100 - float(X[2])
            # store in dict
            rbm_dict[match] = pID

    return rbm_dict


def get_category(bmat, core_thresh):
    ''' Read bmat file and calculate the category of each cluster '''

    print('\n\tReading binary matrix file ...')

    pan_category = {} # by cluster. Core, Accessory, Specific 

    with open(bmat, 'r') as f:
        header = f.readline().rstrip().split()
        genomes = header # genome names are the header
        n = len(genomes) # This is the number of genomes ie(100)
        # Define gene category cutoffs
        core, specific = core_thresh, 1/n # default core = 0.90
        # each cluster is a line in the bmat file
        # compute gene category for each cluster
        for l in f:
            X = l.rstrip().split('\t')
            cluster = X[0] # Cluster # is the first column of bmat
            gs = [int(i) for i in X[1:]] # genes in cluster count
            g = sum(gs) # This is how many genomes cluster is in
            v = g / n # percent genomes with gene
            if v == specific: pan_category[cluster] = ['Specific', v]
            elif v < core: pan_category[cluster] = ['Accessory', v]
            elif v >= core: pan_category[cluster] = ['Core', v]
            else: print(l, v, 'error in assigning gene category')

    return pan_category


def parse_cluster_data(mmsq, pan_category, rbm_dict):

    print('\n\tReading mmseqs cluster file ...')

    # cluster data is a dict of dicts. Each cluster rep gene is a dict with
    # each gene in the cluster as a dict
    cluster_data = defaultdict(lambda: defaultdict())
    data = {}
    total_genes = 0
    mmseq_norbm = 0
    dup_count = 0

    # read the cluster and alignment data
    with open(mmsq, 'r') as f:
        for l in f:
            total_genes += 1
            X = l.rstrip().split('\t')
            cluster = X[0]
            gene = X[1]
            genome = '_'.join(gene.split('_')[:-1])
            pcat = pan_category[cluster] # [pancat, n/N]

            # look for highly conserved genes
            # dist = 1 - sequence alignmend id from RBM file
            # mmseqs seemed to be reporting low alignment scores.
            # for instance a multiple sequence alignment would show a gene
            # cluster is all 100% similar but the mmseq alignment would 0.6633.
            # The RBM from blast would say 100%. So we use the RBM score.
            # store only one gene per genome for each cluster.
            # * gene with the smallest distance to cluster rep.

            if cluster == gene:
                dist = 0
            else:
                match = '-'.join(sorted([cluster, gene]))
                dist = rbm_dict.get(match, -1)
                # if the RBM is not found between a gene in the cluster
                # and the cluster representative (cluster vs gene), then
                # the most likely scenario is that the genome the gene 
                # belongs to has a duplicate, and the duplicate is the 
                # better match. We keep only the distance from RBMs and
                # we discard any cluster-gene match that does not have
                # an RBM.

            if dist == -1:
                mmseq_norbm += 1
                continue

            # here, if a genome has more than one gene in a cluster, we are
            # using the gene with the smaller dist score.
            elif genome in cluster_data[cluster]:
                test_dist = cluster_data[cluster][genome][1]
                if dist < test_dist:
                    dup_count += 1
                    cluster_data[cluster][genome] = [gene, dist]
            else:
                cluster_data[cluster][genome] = [gene, dist]

            # store gene data to for outfile
            data[gene] = [gene, cluster, pcat[0], pcat[1], dist]

    print(f'\t\tTotal genes in mmseqs cluster file: {total_genes}')
    print(f'\t\tGenes in mmseqs clusters without RBM match: {mmseq_norbm}')
    print(f'\t\tPercent of total genes: {mmseq_norbm/total_genes*100:.2f}%')
    print('\t\tLower is better. If over 10% I might look into it further.')
    print(f'\t\tDuplicate genes from same genome in cluster: {dup_count}')
    print(f'\t\tPercent of total genes: {dup_count/total_genes*100:.2f}%')
    print('\t\tLower is better. We do not expect frequent duplicates.')

    return cluster_data, data


def get_avg_dists(cluster_data, pan_category):

    print('\n\tComputing average within cluster sequence distances ...')

    cluster_dist = {'Cluster': [], 'Avg_dist': []}
    cluster_genes = defaultdict(list)

    # create cluster distance distributions for core genes only.
    # each core cluster has one gene from each genome at this point.
    # core genes are 60% - 100% of the genomes.
    # average distance of genes in cluster = sum(dist) / len(dist)
    for cluster, genomes in cluster_data.items():
        distances = []
        pcat = pan_category[cluster]
        if pcat[0] == 'Core':
            for genome, values in genomes.items():
                gene = values[0]
                dist = values[1]
                distances.append(dist)
                cluster_genes[cluster].append(gene)
            cluster_dist['Cluster'].append(cluster)
            cluster_dist['Avg_dist'].append(sum(distances)/len(distances))

    return cluster_dist, cluster_genes


def plot_dist(input_array, qt, xmax, out):

    print('\n\tPlotting average within cluster sequence distances ...')

    # define output file name for the plot
    plot_out = out.split('.')[0] + '.pdf'

    # We want to get a quantile to grab the most conserved genes.
    # With highly similar genomes, the distribution may not be gaussian.
    # there can be a peak at 0 from super high sequence similarities.
    # before computing the quantile we remove values = 0.
    zarray = [i for i in input_array if i != 0]
    qarray = sorted(zarray)
    qtv = np.quantile(qarray, qt)

    fig, ax = plt.subplots(figsize=(7, 5))

    if xmax:
        input_array = [i for i in input_array if i <= xmax]
        
    ax.hist(input_array, bins=75, edgecolor='k', color='#08519c', alpha=0.6,)
    ax.axvline(qtv, color='#ae017e', linestyle='dashed', linewidth=1)
    ax.set_xlabel('Average sequence distance of cluster (%)', fontsize=14)
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

    return qtv


def build_the_list(bmat, mmsq, rbms, core_thresh, qt, xmax, outf):
    ''' Reads in the files, builds the list, and writes to out '''

    # read in the binary matrix file and find Pan Cat for each cluster
    pan_category = get_category(bmat, core_thresh)
    # read in the all vs all rbm data
    rbm_dict = parse_rbm_file(rbms)
    # read in the cluster file and get distances for each cluster
    cluster_data, data = parse_cluster_data(mmsq, pan_category, rbm_dict)
    # get average distance for each cluster and cluster, gene names
    cluster_dist, cluster_genes = get_avg_dists(cluster_data, pan_category)
    # define highly conserved genes lower qt quantile of sequence differences.
    qtv = plot_dist(cluster_dist['Avg_dist'], qt, xmax, outf)
    tmp_data = zip(cluster_dist['Cluster'], cluster_dist['Avg_dist'])
    for clust, dist in tmp_data:
        if dist <= qtv:
            for gene in cluster_genes[clust]:
                data[gene][2] = 'Conserved'
    # write out the data
    with open(outf, 'w') as o:
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
        '-b', '--binary_matrix_file',
        help='Please specify the binary.tsv input file!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-m', '--mmseqs_cluster_file',
        help='Please specify the Mmseqs2 Cluster TSV File!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-r', '--RBM_allvall_file',
        help='Please specify the all vs all RBM file!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-c', '--core_gene_threshold',
        help='(OPTIONAL) Core fraction of genomes with gene (default = 0.90)',
        metavar='',
        type=float,
        required=False,
        default=0.90
        )
    parser.add_argument(
        '-qt', '--quantile_cutoff',
        help='(OPTIONAL) Specify the quantile cutoff (default = 0.1).',
        metavar='',
        type=float,
        default=0.1,
        required=False
        )
    parser.add_argument(
        '-xmax', '--max_xaxis_value',
        help='(OPTIONAL) Specify the max x-axis value (default = None).',
        metavar='',
        type=float,
        default=None,
        required=False
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
    print('\n\tRunning the script ...')

    # define input params
    bmat = args['binary_matrix_file']
    mmsq = args['mmseqs_cluster_file']
    rbms = args['RBM_allvall_file']
    core_thresh = args['core_gene_threshold']
    qt = args['quantile_cutoff']
    xmax = args['max_xaxis_value']
    outf = args['output_file']

    _ = build_the_list(bmat, mmsq, rbms, core_thresh, qt, xmax, outf)

    print(f'\n\nComplete success space cadet!! Finished without errors.\n\n')

if __name__ == "__main__":
    main()

#!/usr/bin/env python

''' Convert MMSeqs2 Cluster TSV file to binary matrix

The binary matrix stores gene presence/absence data for each genome.
Gene clusters are the rows and genomes are the columns.
Zero indicates the gene cluster is absent from a genome.

input: MMSeqs2 Cluster TSV file.
output: Binary matrix as TSV file.

Genome names are inferred from the gene sequence names.
Modify line 93 genome variable to select the right thing.

The cluster key output file just links the cluster number assigned
to rows of the binary matrix output to the representative sequence
for that cluster.

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: December, 2021
License :: GNU GPLv3
Copyright 2021 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse
from collections import defaultdict
import pandas as pd


def generate_binary_matrix(byCluster, byGenome):

    '''
    This function takes two dictionaries with cluster data stored
    by cluster representative name and by genome name.
    It iterates through the set of cluster representative names
    for each genome entry and checks to see if the genome is present
    or absent in that gene cluster. It encodes a binary datatable with
    gene clusters as rows and genomes as columns. Zero indicates a 
    gene cluster is absent from a genome.
    '''

    print('\n\nGenerating binary matrix ...')
    #initialize variables
    # get list of all the cluster ids
    cid = list(byCluster.keys())
    # create new list to use for index. cleaner names: Cluster_#
    # used to name the rows Cluster_01 through n but changed my mind
    #index = [f'Cluster_{i}' for i, c in enumerate(cid)]
    # store binary encoding in dict of {genome: [binary]}
    matrix = defaultdict(list)
    # read through data and encode binary info
    for genome, clusters in byGenome.items():
        for c in cid:
            if c in clusters: matrix[genome].append(1)
            else: matrix[genome].append(0)
    # create dataframe from binary matrix encoding
    # used to rename the cluster name but changed my mind
    #binary_matrix = pd.DataFrame(matrix, index=index)
    binary_matrix = pd.DataFrame(matrix, index=cid)
    # create key dictionary of {index: cid}
    #key = dict(zip(index, cid))

    return binary_matrix#, key


def parse_mmseqs_cluster_tsv(infile):

    '''
    The mmseqs cluster tsv file is two columns.
    column 1 is the cluster representative sequence.
    This sequence is repeated for every sequence in the cluster.
    Column 2 are sequences in the cluster.
    The general idea with this function is to create 2 dictionaries.
    One dictionary stores data by the cluster represenative name as in
    {cluster: genome}. Which genomes are in each cluster.
    The other dicitonary stores data by the genome name.
    {genome: cluster}. Which clusters are in each genome.
    '''

    print('\n\nParsing MMSeqs2 Cluster TSV file ...')
    # initialize variables
    # Store data in two dictionaries.
    byCluster = defaultdict(lambda: defaultdict(int))
    byGenome = defaultdict(lambda: defaultdict(int))

    # read through the file and populate the dictionaries
    with open(infile, 'r') as file:
        for line in file:
            X = line.rstrip().split('\t')
            cluster = X[0]
            genome = '_'.join(X[1].split('_')[:-2])
            # store data in dict of dict
            byCluster[cluster][genome] += 1
            byGenome[genome][cluster] += 1

    return byCluster, byGenome


def main():
    #byGenome = defaultdict(lambda: defaultdict(int))

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-i', '--mmseqs_cluster_tsv_input_file',
        help='Please specify the mmseqs cluster tsv input file!',
        metavar=':',
        type=str,
        required=True
        )
    parser.add_argument(
        '-o', '--binary_matrix_tsv_output_file',
        help='Please specify the name to use for the output file!',
        metavar=':',
        type=str,
        required=True
        )

    args=vars(parser.parse_args())

    # define the input parameters
    infile = args['mmseqs_cluster_tsv_input_file']
    outfile = args['binary_matrix_tsv_output_file']
    #keyout = args['GeneCluster_RepGeneName_key']

    # parse the input file
    byCluster, byGenome = parse_mmseqs_cluster_tsv(infile)
    # generate the binary matrix
    #binary_matrix, key = generate_binary_matrix(byCluster, byGenome)
    binary_matrix = generate_binary_matrix(byCluster, byGenome)
    # write the binary matrix to a tsv file
    binary_matrix.to_csv(outfile, sep='\t', index=True, header=True)
    

    print('\n\nComplete success space cadet! Finished without errors.\n\n')

if __name__ == "__main__":
    main()


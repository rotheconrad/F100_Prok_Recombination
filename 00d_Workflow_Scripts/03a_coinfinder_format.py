#!/usr/bin/env python

''' Reads binary pangenome matrix and writes out Coinfinder tsv file.

Coinfinder accepts a tab-delimetered list of genes present in each
strain as input. An example of a tab-delimited list of genes:

gene_1  genome_1
gene_1  genome_2
gene_1  genome_3
gene_2  genome_2
gene_2  genome_3
gene_3  genome_1
gene_3  genome_2

https://github.com/fwhelan/coinfinder

The Recombinant Genes Analysis from PART 03, Step 01 outputs a binary
pangenome matrix. This script turns that into the tsv Coinfinder 
input file. Simple.

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: June 2023
License :: GNU GPLv3
Copyright 2023 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse


def generate_coinfinder_file(infile, outfile):

    with open(infile, 'r') as inf, open(outfile, 'w') as ouf:
        # genome names are in the header column
        genome_names = inf.readline().rstrip().split('\t')[1:]

        # each line is gene. 1 gene is in the genome. 0 gene is absent
        for line in inf:
            X = line.rstrip().split('\t')
            gene = X[0]
            genomes = [int(i) for i in X[1:]]
            for i, value in enumerate(genomes):
                gname = genome_names[i]
                if value == 1:
                    ouf.write(f'{gene}\t{gname}\n')

    return True


def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-i', '--input_binary_matrix',
        help='Please specify the binary_pangenome_matrix.tsv input file!',
        #metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-o', '--output_file_name',
        help='What do you want to name the output file?',
        #metavar='',
        type=str,
        required=True
        )
    args=vars(parser.parse_args())

    # Do what you came here to do:
    print('\n\nRunning Script...')
    
    # define parameters
    infile = args['input_binary_matrix']
    outfile = args['output_file_name']

    _ = generate_coinfinder_file(infile, outfile)

    print('\n\nComplete success space cadet! Hold on to your boots.\n\n')

if __name__ == "__main__":
    main()

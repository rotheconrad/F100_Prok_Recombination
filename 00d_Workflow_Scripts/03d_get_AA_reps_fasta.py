 #!/usr/bin/env python

'''Renames fasta deflines sequentially using the filename.

Genomes downloaded from NCBI have a typical naming scheme like this:
example filename: GCF_000007105.1_ASM710v1_genomic.fna

The default behavior of this script is cut the second underscore position
("_") and use it as a prefix for renaming the fasta deflines in numeric
consecutive order.

So for a fasta file with three contigs or chromosomes the script cuts
"ASM710v1" from filename and renames fasta deflines as:

>ASM710v1_1
>ASM710v1_2
>ASM710v1_n

Alternatively, the user can input their own desired prefix.

This script effectively renames the file inplace.
eg. Writes temporary file while renaming and then replaces the original.

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: April 2022
License :: GNU GPLv3
Copyright 2022 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse, subprocess

def read_fasta(fp):
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))

def parse_repseqs_fasta(repseqs):
    ''' reads repseq fasta and returns dict of {GeneName: ''} '''

    repNames = {}

    with open(repseqs, 'r') as file:
        for name, seq in read_fasta(file):
            n = name[1:].split(' ')[0]
            repNames[n] = ''

    return repNames


def get_rep_gene_fasta(aaseqs, repNames, outfile):
    ''' read the amino acid sequence fasta and write out genes in repNames '''

    with open(aaseqs, 'r') as file, open(outfile, 'w') as out:
        for name, seq in read_fasta(file):
            n = name[1:].split(' ')[0]
            if n in repNames:
                out.write(f'{name}\n{seq}\n')

    return True


def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-r', '--representative_gene_fasta',
        help='Please specify the representative gene fasta!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-a', '--amino_acid_fasta',
        help='Please specify the amino acid sequence fasta!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-o', '--output_file',
        help='Please specify the output fasta file name!',
        metavar='',
        type=str,
        required=True
        )
    args=vars(parser.parse_args())

    # define input parameters
    repseqs = args['representative_gene_fasta']
    aaseqs = args['amino_acid_fasta']
    outfile = args['output_file']

    print(f'\n\nRenaming fasta file deflines ...')

    # get the gene names of the representative genes
    repNames = parse_repseqs_fasta(repseqs)
    # read the amino acid sequence fasta and write out genes in repNames
    _ = get_rep_gene_fasta(aaseqs, repNames, outfile)

    print(f'\n\nFasta file deflines renamed successfully!\n\n')
    
if __name__ == "__main__":
    main()


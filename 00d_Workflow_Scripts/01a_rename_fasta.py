 #!/usr/bin/env python

'''Renames fasta deflines sequentially using the filename.

This script renames the deflines of fasta files with a genome prefix
and the number of the contig in sequential order.

This script renames the file inplace ie. it overwrites the orginal file.

This script has three modes:

1) Default mode:

Renames all fasta deflines with the file name as the genome prefix.
Any underscores in the file name are replaced with a period.

So My_genome_file.fasta will get renamed to:

>My.genome.file_1
AACGAAGTTGCTGACGGCGGAAGCGACATAGGGATCTGTCAGTTGTCATTCGCGAAAAACATCCGTCCCCGA
>My.genome.file_2
GAAGCCTAGGGGAACAGGTTAGTTTGAGTAGCTTAAGAATGTAAATTCTGGGATTATAGTGTAGTAATCTCT
>My.genome.file_3
AATTAACGGTGACGGTTTTAAGACAGGTCTTCGCAAAATCAAGCGGGGTGATTTCAACAGATTTTGCTGATG

#######################

2) NCBI mode:

Genomes downloaded from NCBI have a typical naming scheme like this:
example filename: GCF_000007105.1_ASM710v1_genomic.fna

The shortest unique description for these files is the "ASM710v1" at the
third underscore position. This mode cuts the third underscore position
and uses it as the genome prefix.

To use this mode set -ncbi True.

So GCF_000007105.1_ASM710v1_genomic.fna fasta will get renamed to:

>ASM710v1_1
AACGAAGTTGCTGACGGCGGAAGCGACATAGGGATCTGTCAGTTGTCATTCGCGAAAAACATCCGTCCCCGA
>ASM710v1_2
GAAGCCTAGGGGAACAGGTTAGTTTGAGTAGCTTAAGAATGTAAATTCTGGGATTATAGTGTAGTAATCTCT
>ASM710v1_n
AATTAACGGTGACGGTTTTAAGACAGGTCTTCGCAAAATCAAGCGGGGTGATTTCAACAGATTTTGCTGATG

#######################

3) User custom mode:

The user can input their own desired genome prefix using the -user flag.

So My_genome_file.fasta using -user myNameScheme will get renamed to:

>myNameScheme_1
AACGAAGTTGCTGACGGCGGAAGCGACATAGGGATCTGTCAGTTGTCATTCGCGAAAAACATCCGTCCCCGA
>myNameScheme_2
GAAGCCTAGGGGAACAGGTTAGTTTGAGTAGCTTAAGAATGTAAATTCTGGGATTATAGTGTAGTAATCTCT
>myNameScheme_n
AATTAACGGTGACGGTTTTAAGACAGGTCTTCGCAAAATCAAGCGGGGTGATTTCAACAGATTTTGCTGATG

* ATTENTION!! DO NOT USE UNDERSCORES IN THE USER DEFINED GENOME PREFIX

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

def Fasta_rename_sequences(infile, prefix):

    outfile = infile + '.rename'

    with open(infile, 'r') as f, open(outfile, 'w') as o:
        i = 1
        for name, seq in read_fasta(f):
            newName = '>%s_%d\n%s\n' % (prefix, i, seq)
            o.write(newName)
            i += 1

    _ = subprocess.run(['mv', outfile, infile])

    return True

def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-i', '--input_file',
        help='Please specify the input fasta file!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-ncbi', '--NCBI_mode',
        help='(OPTIONAL) Set -ncbi True to activate NCBI mode!',
        metavar='',
        type=str,
        required=False,
        default=None
        )
    parser.add_argument(
        '-user', '--user_prefix',
        help='(OPTIONAL) Please specify a naming prefix!',
        metavar='',
        type=str,
        required=False,
        default=None
        )
    args=vars(parser.parse_args())

    # define input parameters
    infile = args['input_file']
    ncbi_mode = args['NCBI_mode']
    user_prefix = args['user_prefix']

    print(f'\n\nRenaming fasta file deflines ...')
    # check mode
    if user_prefix:
        prefix = user_prefix
    elif ncbi_mode:
        prefix = infile.split('/')[-1].split('_')[2]
    else:
        prefix = infile.split('/')[-1].split('.')[0].replace("_", ".")
        
    # rename the fasta deflines
    _ = Fasta_rename_sequences(infile, prefix)

    print(f'\n\nFasta file deflines renamed successfully!\n\n')
    
if __name__ == "__main__":
    main()


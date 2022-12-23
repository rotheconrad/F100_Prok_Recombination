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
        '-p', '--prefix',
        help='(OPTIONAL) Please specify a naming prefix!',
        metavar='',
        type=str,
        required=False,
        default=None
        )
    args=vars(parser.parse_args())

    # define input parameters
    infile = args['input_file']
    prefix = args['prefix']

    print(f'\n\nRenaming fasta file deflines ...')
    # check for prefix
    if not prefix:
        prefix = infile.split('/')[-1].split('_')[2]
    # rename the fasta deflines
    _ = Fasta_rename_sequences(infile, prefix)

    print(f'\n\nFasta file deflines renamed successfully!\n\n')
    
if __name__ == "__main__":
    main()


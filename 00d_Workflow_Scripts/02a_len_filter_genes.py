#!/usr/local/pacerepov1/python/2.7/bin/python

 #!/usr/bin/env python

'''Filters CDS fasta and removes short and long sequences

Intended for gene prediction fasta files such as Prodigal output.

Removes predicted CDS sequences shorter than 100 base pairs or longer
than 40,000 base pairs. Default minimum and maximum length parameters
may be changed by the user.

Caution:
This script effectively filters the file inplace.
eg. Writes temporary file then replaces the original.

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

def Fasta_filter_sequences(infile, min_len, max_len):

    outfile = infile + '.lenfilter'
    small = 0
    large = 0
    c = 0
    i = 0

    with open(infile, 'r') as f, open(outfile, 'w') as o:
        for name, seq in read_fasta(f):
            c += 1
            slen = len(seq)
            if slen < min_len:
                small += 1
            elif slen > max_len:
                large += 1
            else:
                o.write(f'{name}\n{seq}\n')
                i += 1

    _ = subprocess.run(['mv', outfile, infile])

    filtered = c - i
    print(
        f'\n\nWrote {i} of {c} genes to file.\n'
        f'Genes filtered: {filtered}.\n'
        f'Genes < {min_len} bp: {small}\n'
        f'Genes > {max_len} bp: {large}'
        )

    return True

def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-i', '--input_file',
        help='Please specify the input file!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-min', '--min_gene_length',
        help='(OPTIONAL) Specify minimum gene length (default: 100)',
        metavar='',
        type=int,
        required=False,
        default=100
        )
    parser.add_argument(
        '-max', '--max_gene_length',
        help='(OPTIONAL) Specify minimum gene length (default: 40000)',
        metavar='',
        type=int,
        required=False,
        default=40000
        )
    args=vars(parser.parse_args())

    infile = args['input_file']
    min_len = args['min_gene_length']
    max_len = args['max_gene_length']
    _ = Fasta_filter_sequences(infile, min_len, max_len)

    print(f'\nFasta file filtered successfully!\n\n')
    
if __name__ == "__main__":
    main()


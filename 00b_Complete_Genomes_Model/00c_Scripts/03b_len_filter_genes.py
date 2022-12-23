#!/usr/local/pacerepov1/python/2.7/bin/python

 #!/usr/bin/env python

'''Filters CDS fasta and removes short and long sequences

Intended for gene prediction fasta files.

Removes sequences < 100bp and >40000 bp.

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

def Fasta_filter_sequences(infile):

    outfile = infile + '.lenfilter'
    small = 0
    large = 0
    c = 0
    i = 0

    with open(infile, 'r') as f, open(outfile, 'w') as o:
        for name, seq in read_fasta(f):
            c += 1
            slen = len(seq)
            if slen < 100:
                small += 1
            elif slen > 40000:
                large += 1
            else:
                o.write(f'{name}\n{seq}\n')
                i += 1

    _ = subprocess.run(['mv', outfile, infile])

    filtered = c - i
    print(
        f'\n\nWrote {i} of {c} genes to file.\n'
        f'Genes filtered: {filtered}.\n'
        f'Genes < 100bp: {small}\n'
        f'Genes > 40000bp: {large}'
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
    args=vars(parser.parse_args())

    infile = args['input_file']
    _ = Fasta_filter_sequences(infile)

    print(f'\nFasta file filtered successfully!\n\n')
    
if __name__ == "__main__":
    main()


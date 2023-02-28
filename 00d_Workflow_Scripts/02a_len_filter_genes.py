#!/usr/local/pacerepov1/python/2.7/bin/python

 #!/usr/bin/env python

'''Filters CDS fasta and removes short and long sequences

Intended for gene prediction fasta files such as Prodigal output.

Looks at predicted CDS length distribution from input file and removes
the top and bottom 0.5% of gene lengths. Use -qnt to change the default of
0.01. Since the length filter is determined from the data, this setting
works with nucleotide or amino acid sequences.

Writes out histogram of the gene length distribution. Notes the x-axis
units will be in amino acid residue length for amino acid sequence fasta
files rather than the default label of base pairs (bp).

(OPTIONAL) set -ul True to remove predicted CDS sequences based on
a minimum and maximum sequence length in base pairs. default min of 100
base pairs and max of 40,000 base pairs. Use -min and -max to change
the default lengths. If the input is amino acid sequence, divide the min
and max bp filters by 3.

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
import numpy as np
import matplotlib.pyplot as plt


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


def parse_gene_fasta(infile):

    data = {}

    with open(infile, 'r') as f:
        for name, seq in read_fasta(f):
            gene = f'{name}\n{seq}\n'
            glen = len(seq)
            data[gene] = glen

    return data


def plot_dist(data, out, min_qnt, max_qnt):

    plot_out = out.split('.')[0] + '_geneLenDist.pdf'
    input_array = list(data.values())
    qarray = sorted(input_array)
    qmn, qmx = np.quantile(qarray, min_qnt), np.quantile(qarray, max_qnt)

    fig, ax = plt.subplots(figsize=(7, 5))

    ax.hist(input_array, bins=75, edgecolor='k', color='#08519c', alpha=0.6,)
    ax.axvline(qmn, color='#ae017e', linestyle='dashed', linewidth=1)
    ax.axvline(qmx, color='#ae017e', linestyle='dashed', linewidth=1)
    ax.set_xlabel('Gene Length (bp)', fontsize=14)
    ax.set_ylabel('Count', fontsize=14)
    mline = f'Min: {qmn:.2f}\nMax: {qmx:.2f}'
    ax.text(0.99, 0.99, mline, transform=ax.transAxes, ha='right', va='top')

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

    return qmn, qmx


def Fasta_filter_sequences(data, infile, min_len, max_len):

    outfile = infile + '.lenfilter'
    small = 0
    large = 0
    c = 0
    i = 0

    with open(outfile, 'w') as o:
        for line, glen in data.items():
            c += 1
            gname = line.split('\n')[0].split(' ')[0][1:]
            if glen < min_len:
                print(f'{gname} to short: {glen}')
                small += 1
            elif glen > max_len:
                print(f'{gname} to long: {glen}')
                large += 1
            else:
                o.write(line)
                i += 1

    _ = subprocess.run(['mv', outfile, infile])

    # write summary out
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
        '-qnt', '--quantile_gene_length',
        help='(OPTIONAL) Specify filter quantile (default: 0.005)',
        metavar='',
        type=float,
        required=False,
        default=0.005
        )
    parser.add_argument(
        '-ul', '--use_min_max',
        help='(OPTIONAL) use min and max gene length (default: False)',
        metavar='',
        type=str,
        required=False,
        default=False
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

    # define input params
    infile = args['input_file']
    min_qnt = args['quantile_gene_length']
    max_qnt = 1 - min_qnt
    ul = args['use_min_max']
    min_len = args['min_gene_length']
    max_len = args['max_gene_length']

    # read in the fasta sequences
    data = parse_gene_fasta(infile)

    # plot the gene length distribution
    # calculates and returns quantile lengths
    lmin, lmax = plot_dist(data, infile, min_qnt, max_qnt)

    # if ul is true use min and max length for the filter instead of quantiles
    if ul:
        lmin, lmax = min_len, max_len

    # filter the sequences and write the output file.
    _ = Fasta_filter_sequences(data, infile, lmin, lmax)

    print(f'\nFasta file filtered successfully!\n\n')
    
if __name__ == "__main__":
    main()


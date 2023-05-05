#!/usr/local/pacerepov1/python/2.7/bin/python

 #!/usr/bin/env python

'''Filters CDS fasta and removes short and long sequences

Intended for gene prediction fasta files such as Prodigal output.

This script filters a fasta file by minimum and maximum gene length in 
base pairs for nucleotide sequences with a min of 100bp and a max of
8000bp by default. You can use the -min and -max parameter flags to
change these defaults.

(OPTIONAL) set -aa True for amino acid input which will divide the min
and max by 3 before filtering.

(CAUTION) This script filters inplace e.g. it replaces the input file.

This script writes gene name, gene length, and filter summary for the
genes that are removed by the filter. This script also writes out a
histogram of the gene length distribution and labels the top and bottom
0.5% quantile lengths. This plot can used to inform custom length
selection. Use -qnt to change the 0.005 default.

(OPTIONAL) set -uq True to use the quantile lengths for the min and max
gene lengths of the filter.

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


def plot_dist(data, out, min_qnt, max_qnt, aa):

    plot_out = out.split('.')[0] + '_geneLenDist.pdf'
    input_array = list(data.values())
    qarray = sorted(input_array)
    qmn, qmx = np.quantile(qarray, min_qnt), np.quantile(qarray, max_qnt)

    fig, ax = plt.subplots(figsize=(7, 5))

    ax.hist(input_array, bins=75, edgecolor='k', color='#08519c', alpha=0.6,)
    ax.axvline(qmn, color='#ae017e', linestyle='dashed', linewidth=1)
    ax.axvline(qmx, color='#ae017e', linestyle='dashed', linewidth=1)
    xlab = 'Gene length (bp)' if aa == False else 'Gene length (aa)'
    ax.set_xlabel(xlab, fontsize=14)
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
            if glen <= min_len:
                print(f'{gname} to short: {glen}')
                small += 1
            elif glen >= max_len:
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
        '-plot', '--plot_gene_length_distribution',
        help='(OPTIONAL) Set true to plot histogram of gene lengths (default: False)',
        metavar='',
        type=str,
        required=False,
        default=False
        )
    parser.add_argument(
        '-aa', '--amino_acid_input',
        help='(OPTIONAL) For amino acid input, divide legnth filter by 3 (default: False)',
        metavar='',
        type=str,
        required=False,
        default=False
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
        '-uq', '--use_qnt_min_max',
        help='(OPTIONAL) use quantiles for min and max gene length (default: False)',
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
        help='(OPTIONAL) Specify minimum gene length (default: 8000)',
        metavar='',
        type=int,
        required=False,
        default=8000
        )
    args=vars(parser.parse_args())

    # define input params
    infile = args['input_file']
    plot = args['plot_gene_length_distribution']
    aa = args['amino_acid_input']
    min_qnt = args['quantile_gene_length']
    max_qnt = 1 - min_qnt
    uq = args['use_qnt_min_max']
    min_len = args['min_gene_length']
    max_len = args['max_gene_length']

    # read in the fasta sequences
    data = parse_gene_fasta(infile)

    # plot the gene length distribution
    # calculates and returns quantile lengths
    if plot:
        lmin, lmax = plot_dist(data, infile, min_qnt, max_qnt, aa)

    # if uq is true use min and max quantile length for the filter
    if uq and plot:
        min_len, max_len = lmin, lmax
    elif uq and not plot:
        input_array = list(data.values())
        qarray = sorted(input_array)
        lmin, lmax = np.quantile(qarray, min_qnt), np.quantile(qarray, max_qnt)
        min_len, max_len = lmin, lmax
    # if aa is true divide min max length by 3.
    if aa:
        min_len, max_len = int(min_len / 3), int(max_len / 3)

    # filter the sequences and write the output file.
    _ = Fasta_filter_sequences(data, infile, min_len, max_len)

    print(f'\nFasta file filtered successfully!\n\n')
    
if __name__ == "__main__":
    main()


 #!/usr/bin/env python

'''This script reads the pancat_file.tsv file and outputs a boxplot
displaying the fraction of total genes per pangenome class per genome.
It also writes a data table to file of the boxplot data.

Each gene cluster is scored in four pangenome classes:
    - genome specific (n = 1)
    - rare (n/N ≤ 0.25)
    - common (0.25 < n/N < 0.9)
    - core (n/N ≥ 0.9)

    where n = number of genomes is set with that gene cluster
    and N = total number of genomes in the set

The fraction of each pangenome class for each genome is computed:
    - genome specific in genome / total genes in genome
    - rare in genome / total genes in genome
    - common in genome / total genes in genome
    - core in genome / total genes in genome    

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: Oct 2023
License :: GNU GPLv3
Copyright 2023 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse
from collections import defaultdict
from collections import Counter
import pandas as pd
import matplotlib.pyplot as plt


def parse_pancat_file(infile, outpre, rgc, cgc):

    # first we parse the pangenome categories for all genomes into
    # a list of all genes for each genome and store in data dict
    # data = {'genome name': [gene class labels]}
    data1 = defaultdict(list)

    with open(infile, 'r') as file:
        header = file.readline()
        for line in file:
            X = line.rstrip().split()
            genome = X[0].split('_')[0]
            label = X[2]
            nN = float(X[3])

            if label == 'Accessory' or label == 'Core' or label == 'Conserved':
                if nN >= cgc: label = 'Core'
                elif nN > rgc: label = 'Common'
                else: label = 'Rare'

            data1[genome].append(label)

    # next we compute the percent each gene class represents in each genome
    data2 = {
                'Genome Name': [],
                'Gene Counts': [],
                'Specific': [],
                'Rare': [],
                'Common': [],
                'Core': [],
                }

    for genome, genes in data1.items():
        data2['Genome Name'].append(genome)
        n = len(genes) # gene count per genome
        data2['Gene Counts'].append(n)
        counts = Counter(genes)
        # Counter returns dict of unique value counts
        # {'Core': 2528, 'Common': 345, 'Rare': 232, 'Specific': 100}
        # divide count by genome gene count and store in data2

        data2['Specific'].append(round(counts.get('Specific', 0)/n * 100, 2))
        data2['Rare'].append(round(counts.get('Rare', 0)/n * 100, 2))
        data2['Common'].append(round(counts.get('Common', 0)/n * 100, 2))
        data2['Core'].append(round(counts.get('Core', 0)/n * 100, 2))

    df = pd.DataFrame(data2)

    df.to_csv(f'{outpre}_data.tsv', sep='\t', index=False)

    return df


def build_boxplot(df, outpre):

    fig, (ax1, ax2) = plt.subplots(
                1, 2, figsize=(7, 5), gridspec_kw={'width_ratios': [2, 1]}
                )
    gclass = ['Specific', 'Rare', 'Common', 'Core']
    b1 = df.boxplot(column=gclass, ax=ax1, grid=False)
    b2 = df.boxplot(column=['Gene Counts'], ax=ax2, grid=False) 

    fz = 12 # label font size

    ax1.set_title('Genes per class per genome', fontsize=fz)
    ax2.set_title('Genes per genome', fontsize=fz)
    ax1.set_ylabel('Genes (%)', fontsize=fz)
    ax2.set_ylabel('Genes', fontsize=fz)

    # Print Averages on plot
    avg_line = (
                "Averages:\n"
                f"Genes: {df['Gene Counts'].mean().round(2)}\n"
                f"Core: {df['Core'].mean().round(2)}%\n"
                f"Common: {df['Common'].mean().round(2)}%\n"
                f"Rare: {df['Rare'].mean().round(2)}%\n"
                f"Specific: {df['Specific'].mean().round(2)}%\n"
                )

    ax1.text(
        0.04, 0.98, avg_line, verticalalignment='top', transform=ax1.transAxes
        )


    # ax style
    for ax in [ax1, ax2]:
        # set the axis parameters / style
        ax.minorticks_on()
        ax.tick_params(axis='both', labelsize=12)
        ax.tick_params(axis='x', labelrotation=0)
        ax.tick_params(axis='x', which='minor', bottom=False)
        # set grid style
        ax.yaxis.grid(
            which="minor", color='#f0f0f0', linestyle='--', linewidth=.75
            )
        ax.yaxis.grid(
            which="major", color='#d9d9d9', linestyle='--', linewidth=1
            )
        ax.set_axisbelow(True)

    plt.tight_layout()
    plt.savefig(f'{outpre}_boxplot.pdf')
    plt.close()

    return True


def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-i', '--input_file_name',
        help='Please specify the pancat file!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-o', '--output_file_prefix',
        help='Please specify a prefix for the output file names!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-rgc', '--rare_gene_cutoff',
        help='(OPTIONAL) Specify rare gene cutoff (default = 0.25).',
        metavar='',
        type=float,
        default=0.25,
        required=False
        )
    parser.add_argument(
        '-cgc', '--core_gene_cutoff',
        help='(OPTIONAL) Specify core gene cutoff (default = 0.9).',
        metavar='',
        type=float,
        default=0.9,
        required=False
        )
    args=vars(parser.parse_args())

    # Do what you came here to do:
    print('\n\nRunning Script...\n')

    # define parameters
    infile = args['input_file_name']
    outpre = args['output_file_prefix']
    rgc = args['rare_gene_cutoff']
    cgc = args['core_gene_cutoff']

    df = parse_pancat_file(infile, outpre, rgc, cgc)

    _ = build_boxplot(df, outpre)

    print(f'\n\nComplete success space cadet!! Finished without errors.\n\n')

if __name__ == "__main__":
    main()


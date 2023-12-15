 #!/usr/bin/env python

'''Barplot by metadata group.

Unfortunately the GAM model upper and lower confidence intervals reach
above 1 (y-axis) and 100 (x-axis) right around the genomovar group level
of ≥99.5% ANI and so we don't get significant pairs at this level.

But, for metadata groups below this level, this script can count the
number of genome pairs outside the upper and lower confidence intervals
within and between groups. I'm using phylogroups from the 1st column of
my metadata files for this first draft.

* the phylogroup column can contain any user defined groupings. It does
  not have to be strictly phylogroup. It can be clades or niches or
  sites etc. Any way the user wants to look at different groups. The
  structure is a two column tsv file with genome name and group name.

Input:

- metadata file with genome name and phylogroup column
- sig-pairs file from 02d_f100_scatter_pyGAM.py

Outputs: bar plots as vectorized pdf

- {outpre}_F100_sigpairs-01.pdf high level 4 category overview figure
- {outpre}_F100_sigpairs-02-A.pdf category A by metadata groupings
- {outpre}_F100_sigpairs-03-B.pdf category B by metadata groupings
- {outpre}_F100_sigpairs-04-C.pdf category C by metadata groupings
- {outpre}_F100_sigpairs-05-D.pdf category D by metadata groupings
- f'{outpre}_F100_sigpairs-06-GroupCounts.pdf'

Categories:
    A - Same phylogroup above 97.5% c.i.
    B - Different phylogroup above 97.5% c.i.
    C - Same phylogroup below 2.5% c.i.
    D - Different phylogroup below 2.5% c.i.

Options:

    -norm True: Divides each bar by the total number of genome pairs for that
    grouping and multiplies by 100 so that each bar shows the percent of total
    genome pairs for that grouping that are significant instead of raw counts.
    e.g. If two groups A vs. B1 has a total of 100 genome pairs and the 20 of
    those are found to be significant above the 97.5% c.i. then with -norm True
    the bar would show 20%.

    * The denominators for normalization are computed from the genome counts in
      the metadata file so it is important the metadata file includes all 
      genomes used in the analysis.
    
    * Figure 01 is normalized by the total number of n choose 2 genome
      combinations where is n is the total genomes in the data set. To remove
      self matches this is calculated as: combinations(n, 2) – n.

    * Figures 02 – 05 are normalized by the total number of n choose 2 genome
      combinations where n is the total genomes within each grouping (i.e. A-A,
      B1-A, A-G etc.). To remove self matches this is calculated as:
      combinations(n,2) – n.

    * Figure 06 is normalized by the total number of genomes in the dataset.
      (this figure should be the same for all GAM models for each dataset).

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: Dec 2023
License :: GNU GPLv3
Copyright 2023 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse
from collections import defaultdict
from math import comb
import pandas as pd; import numpy as np
import matplotlib.pyplot as plt


def parse_metadata_file(metadata):
    # reads metadata into md dict of {genome name: phylogroup}
    # return md dict
    md = {}
    gpg_counts = defaultdict(int) # phylogroup genome counts
    # count of genome pairs in each dict used to normalize.
    data03 = defaultdict(int)

    with open(metadata, 'r') as file:
        header = file.readline() # skip header
        for line in file:
            X = line.rstrip().split('\t')
            genome = X[0]
            phylogroup = X[1]
            md[genome] = phylogroup
            gpg_counts[phylogroup] += 1

    # for norm == 'True'
    # count all vs all genome pair combinations total
    genome_total = sum(gpg_counts.values())
    # compute n choose 2 combinations - genome total to remove self matches
    gpair_combinations = comb(genome_total, 2) - genome_total
    data03['total'] = gpair_combinations
    # count all vs all genome pair combinations for each grouping
    for group1, count1 in gpg_counts.items():
        for group2, count2 in gpg_counts.items():
            pgpair = '-'.join(sorted([group1, group2]))
            genome_total = count1 + count2
            gpair_combinations = comb(genome_total, 2) - genome_total
            data03[pgpair] = gpair_combinations

    return md, gpg_counts, data03


def parse_sig_pairs_file(sigpairs, md):
    # reads sig-pairs.tsv file, matches to phylogroups, write dict
    # spg = same phylogroup, dpg = different phylogroup
    # above 97.5% c.i., below 2.5% c.i.
    
    data01 = {'spg_above': 0, 'spg_below': 0, 'dpg_above': 0, 'dpg_below': 0}
    data02 = {
                'spg_above': defaultdict(int),
                'dpg_above': defaultdict(int),
                'spg_below': defaultdict(int),
                'dpg_below': defaultdict(int)
                }

    assigned = 0
    unassigned = 0

    with open(sigpairs, 'r') as file:
        header = file.readline()
        for line in file:
            X = line.rstrip().split('\t')
            gpair = X[0].split('-')
            g1, g2 = gpair[0], gpair[1]
            # prediction. Recombining = above, Non-Recombining = below
            prd = X[3]
            # get phylogroups
            if g1 in md and g2 in md:
                assigned += 1
                pg1, pg2 = md[g1], md[g2]
                pgpair = '-'.join(sorted([pg1, pg2]))

            else:
                unassigned += 1
                continue
            # assign to data category
            if prd == 'Recombining' and pg1 == pg2:
                data01['spg_above'] += 1
                data02['spg_above'][pgpair] += 1
            elif prd == 'Recombining' and pg1 != pg2:
                data01['dpg_above'] += 1
                data02['dpg_above'][pgpair] += 1
            elif prd == 'Non-recombining' and pg1 == pg2:
                data01['spg_below'] += 1
                data02['spg_below'][pgpair] += 1
            elif prd == 'Non-recombining' and pg1 != pg2:
                data01['dpg_below'] += 1
                data02['dpg_below'][pgpair] += 1

    print(f'\n\nGenome pairs assigned to phylogroup: {assigned}')
    print(f'Genome pairs unassigned to phylogroup: {unassigned}')

    return data01, data02


def build_bar_plot_01(data01, data03, outpre, norm):
    # simple bar plot for data01

    labels = ['A', 'B', 'C', 'D']
    counts = [
                data01['spg_above'], data01['dpg_above'],
                data01['spg_below'], data01['dpg_below']
                ]
    if norm == 'True':
        tots = data03['total']
        counts = [i/tots*100 for i in counts]

    fig, ax = plt.subplots(figsize=(7,5))

    ax.bar(
        labels, counts, color='#bdbdbd', edgecolor='black',
        hatch="/", alpha=0.5, #fill=False
        )

    # set plot title
    ax.set_title('Significant genome pairs')
    # change axis labels
    ax.set_xlabel('')
    ylab = 'Genome pairs (%)' if norm else 'Genome pairs'
    ax.set_ylabel(ylab, fontsize=12)
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

    txt = (
        'A - Same phylogroup above 97.5% c.i.\n'
        'B - Different phylogroup above 97.5% c.i.\n'
        'C - Same phylogroup below 2.5% c.i.\n'
        'D - Different phylogroup below 2.5% c.i.'
        )
    ax.annotate(
                txt, # text to annotate
                xy=(0.5, 0), # reference point for placement
                xycoords=('figure fraction', 'figure fraction'),
                xytext=(-100, 5), # distance from reference point
                textcoords='offset points',
                size=12, ha='left', va='bottom',
                annotation_clip=False
                )

    # adjust layout, save, and close
    #fig.set_tight_layout(True)
    plt.subplots_adjust(left=0.12, right=0.98, top=0.95, bottom=0.23)
    fig.savefig(f'{outpre}_F100_sigpairs-01.pdf')
    plt.close() 

    return True


def build_bar_plots_x(ctgy, pgs, data03, outfile, norm):
    # simple bar plots of the ctgy data

    ctgy_switch = {
                'spg_above': 'Same phylogroup above 97.5% c.i.',
                'spg_below': 'Same phylogroup below 2.5% c.i.',
                'dpg_above': 'Different phylogroup above 97.5% c.i.',
                'dpg_below': 'Different phylogroup below 2.5% c.i.',
                'gpg_counts': 'Phylogroup pairing genome counts'
                }
    title = ctgy_switch[ctgy]

    labels, counts = [], []
    for k,v in pgs.items():
        labels.append(k)
        counts.append(v)

    if norm and ctgy == 'gpg_counts':
        tots = sum(counts)
        counts = [i/tots*100 for i in counts]
    elif norm == 'True':
        counts = [j/data03[i]*100 for i,j in zip(labels, counts)]

    fig, ax = plt.subplots(figsize=(7,5))

    ax.bar(
        labels, counts, color='#bdbdbd', edgecolor='black',
        hatch="/", alpha=0.5, #fill=False
        )

    # set plot title
    ax.set_title(f'Breakdown: {title}')
    # change axis labels
    ax.set_xlabel(f'Phylogroup(s)', fontsize=12)
    ylab = "Genomes" if ctgy == 'gpg_counts' else "Genome pairs"
    if norm == 'True': ylab = ylab + " (%)"
    ax.set_ylabel(ylab, fontsize=12)
    # set the axis parameters / style
    ax.minorticks_on()
    ax.tick_params(axis='both', labelsize=12)
    ax.tick_params(axis='x', labelrotation=90)
    ax.tick_params(axis='x', which='minor', bottom=False)
    # set grid style
    ax.yaxis.grid(
        which="minor", color='#f0f0f0', linestyle='--', linewidth=.75
        )
    ax.yaxis.grid(
        which="major", color='#d9d9d9', linestyle='--', linewidth=1
        )
    ax.set_axisbelow(True)

    # adjust layout, save, and close
    fig.set_tight_layout(True)
    #plt.subplots_adjust(left=0.12, right=0.98, top=0.95, bottom=0.23)
    fig.savefig(outfile)
    plt.close() 


    return True


def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-md', '--metadata_file',
        help='Please specify the metadata file!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-sp', '--sig_pairs_file',
        help='Please specify *_sig-pairs.tsv file!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-o', '--output_file_prefix',
        help='Please specify a prefix for the output files!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-norm', '--normalize_by_group',
        help='(OPTIONAL): Set -norm True to divide each bar by group total!',
        metavar='',
        type=str,
        required=False
        )
    args=vars(parser.parse_args())

    # Do what you came here to do:
    print('\n\nRunning Script...')

    # define parameters
    metadata = args['metadata_file']
    sigpairs = args['sig_pairs_file']
    outpre = args['output_file_prefix']
    norm = args['normalize_by_group']

    # parse the data and write updated concatenated files
    print(f'\n\nParsing metadata ...')
    md, gpg_counts, data03 = parse_metadata_file(metadata)

    # parse sig pairs file
    print(f'\n\nParsing sig-pairs ...')
    data01, data02 = parse_sig_pairs_file(sigpairs, md)

    # plotting data
    print(f'\n\nPlotting figure 01 - high level data ...')
    _ = build_bar_plot_01(data01, data03, outpre, norm)

    print(f'\n\nPlotting detailed figures ...')
    for i, (ctgy, pgs) in enumerate(data02.items()):
        switch = {
                    'spg_above': 'A', 'dpg_above': 'B',
                    'spg_below': 'C', 'dpg_below': 'D'
                    }
        print(f'\t\tPlotting {switch[ctgy]}, figure 0{i+2} ...')
        outfile = f'{outpre}_F100_sigpairs-0{i+2}-{switch[ctgy]}.pdf'
        _ = build_bar_plots_x(ctgy, pgs, data03, outfile, norm)

    # plot total phylogroup pair genome counts
    print(f'\n\nPlotting figure 0{i+3} - phylogroup pair genome counts ...')
    # sort by value
    gpg_counts = {
                k: v for k, v in sorted(
                                    gpg_counts.items(),
                                    reverse=True,
                                    key=lambda item: item[1])
                                }
    outfile = f'{outpre}_F100_sigpairs-0{i+3}-GroupCounts.pdf'
    _ = build_bar_plots_x('gpg_counts', gpg_counts, data03, outfile, norm)

    print(f'\n\nComplete success space cadet!! Finished without errors.\n\n')

if __name__ == "__main__":
    main()

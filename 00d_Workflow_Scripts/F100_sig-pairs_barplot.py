 #!/usr/bin/env python

'''Barplot by metadata group.

Unfortunately the GAM model upper and lower confidence intervals reach
above 1 (y-axis) and 100 (x-axis) right around the genomovar group level
of ≥99.5% ANI and so we don't get significant pairs at this level.

But, for metadata groups below this level, this script can count the
number of genome pairs outside the upper and lower c.i.s within and
between groups. I'm using phylogroups from the 1st column of my
metadata files for this first draft.

Input:

- metadata file with genome name and phyllogroup column
- sig-pairs file from 02d_f100_scatter_pyGAM.py

Outputs: bar plots as vectorized pdf

- {outpre}_F100_sigpairs-01.pdf high level 4 category overview figure
- {outpre}_F100_sigpairs-{ctgy}-02.pdf category A by metadata groupings
- {outpre}_F100_sigpairs-{ctgy}-03.pdf category B by metadata groupings
- {outpre}_F100_sigpairs-{ctgy}-04.pdf category C by metadata groupings
- {outpre}_F100_sigpairs-{ctgy}-05.pdf category D by metadata groupings

Categories:
    A - Same phylogroup above 97.5% c.i.
    B - Different phylogroup above 97.5% c.i.
    C - Same phylogroup below 2.5% c.i.
    D - Different phylogroup below 2.5% c.i.

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
import pandas as pd; import numpy as np
import matplotlib.pyplot as plt


def parse_metadata_file(metadata):
    # reads metadata into md dict of {genome name: phylogroup}
    # return md dict
    md = {}

    with open(metadata, 'r') as file:
        for line in file:
            X = line.rstrip().split('\t')
            genome = X[0]
            phylogroup = X[1]
            md[genome] = phylogroup

    return md


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


def build_bar_plot_01(data01, outpre):
    # simple bar plot for data01

    labels = ['A', 'B', 'C', 'D']
    counts = [
                data01['spg_above'], data01['dpg_above'],
                data01['spg_below'], data01['dpg_below']
                ]

    fig, ax = plt.subplots(figsize=(7,5))

    ax.bar(
        labels, counts, color='#bdbdbd', edgecolor='black',
        hatch="/", alpha=0.5, #fill=False
        )

    # set plot title
    ax.set_title('Significant genome pairs')
    # change axis labels
    ax.set_xlabel('')
    ax.set_ylabel("Genome pairs", fontsize=12)
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


def build_bar_plots_x(ctgy, pgs, outfile):
    # simple bar plots of the ctgy data

    ctgy_switch = {
                'spg_above': 'Same phylogroup above 97.5% c.i.',
                'spg_below': 'Same phylogroup below 2.5% c.i.',
                'dpg_above': 'Different phylogroup above 97.5% c.i.',
                'dpg_below': 'Different phylogroup below 2.5% c.i.'
                }
    title = ctgy_switch[ctgy]

    labels, counts = [], []
    for k,v in pgs.items():
        labels.append(k)
        counts.append(v)

    fig, ax = plt.subplots(figsize=(7,5))

    ax.bar(
        labels, counts, color='#bdbdbd', edgecolor='black',
        hatch="/", alpha=0.5, #fill=False
        )

    # set plot title
    ax.set_title(f'Breakdown: {title}')
    # change axis labels
    ax.set_xlabel(f'Phylogroup(s)', fontsize=12)
    ax.set_ylabel("Genome pairs", fontsize=12)
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
    args=vars(parser.parse_args())

    # Do what you came here to do:
    print('\n\nRunning Script...')

    # define parameters
    metadata = args['metadata_file']
    sigpairs = args['sig_pairs_file']
    outpre = args['output_file_prefix']

    # parse the data and write updated concatenated files
    print(f'\n\nParsing metadata ...')
    md = parse_metadata_file(metadata)
    # parse sig pairs file
    print(f'\n\nParsing sig-pairs ...')
    data01, data02 = parse_sig_pairs_file(sigpairs, md)
    # plotting data
    print(f'\n\nPlotting figure 01 - high level data ...')
    _ = build_bar_plot_01(data01, outpre)
    print(f'\n\nPlotting detailed figures ...')
    for i, (ctgy, pgs) in enumerate(data02.items()):
        switch = {
                    'spg_above': 'A', 'dpg_above': 'B',
                    'spg_below': 'C', 'dpg_below': 'D'
                    }
        print(f'\t\tPlotting {switch[ctgy]}, figure 0{i+1} ...')
        outfile = f'{outpre}_F100_sigpairs-0{i+2}-{switch[ctgy]}.pdf'
        _ = build_bar_plots_x(ctgy, pgs, outfile)

    print(f'\n\nComplete success space cadet!! Finished without errors.\n\n')

if __name__ == "__main__":
    main()

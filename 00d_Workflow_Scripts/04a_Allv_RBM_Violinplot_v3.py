#!/usr/bin/env python

'''Violin plots from All vs. All RBM data

Creates violin plots of the Cumulative Recombining gene fraction within groups.

A) compare each genome to other genomes in the same genomovar
B) compare each genome to other genomes within each separate genomovar in the
genomes phylogroup
C) compare each genome to other genomes within each separate genomovar outside
the genomes phylogroup
D) compare each genome to other genomes of a different species
E) compare each genome to all genomes in its phylogroup excluding genomes from
its genomovar
F) compare each genome to all genomes within its species excluding genomes from
its phylogroup

Data are presented in hybrid violin plots where the top and bottom whiskers
show the minimum and maximum values, the middle whisker shows the median value,
the black boxes show the interquartile range, and the shaded light blue regions
show the density of values along the y-axis.

Each datapoint is 1 reference genome against all other genomes in each category
Each reference genome is represented in each category except genomes that
singleton genomovars are excluded from 1 and 2.

If you imagine a one vs many genomes rarefaction plot, each y-axis value in the
violin plot is the value at the end of the rarefaction plot, and the x-axis is
the specific groups of genome used for each rarefaction plot.
 
For A and B, first I collect the pairwise data for each reference genome to all
other genomes in its own phylogroup. Then I split these pairwise data for each
reference genome into the same genomovar or different genomovar (but not all
different genomovars in the species, only different genomovars in the same
phylogroup because the first step).
 
A) each reference genome is compared one to many genomes with genomes in
        its own genomovar (I think this one is clear)
 
B) each reference genome is compared one to many genomes with genomes in
        the other genomovars of its own phylogroup, but only one genomovar at
        a time. This excludes genomes from its own genomovar. (e.g. each 
        reference genome in genomovar A is compared to all genomes in 
        genomovar B, then genomovar C, then genomovar D for all genomovars in 
        the phylogroup with more than 1 genome).
 
* A and B exclude singleton genomovars.
 
###
 
For E and C, first I collect the pairwise data for each reference genome to all
other genomes that are NOT in the same genomovar. Then I split these pairwise 
data for each reference genome into the same phylogroup or different phylogroup
 
E) each reference genome is compared one to many genomes with genomes in
        its own phylogroup. This excludes genomes in the same genomovar but it 
        includes the genomes from singleton genomovars. Whereas B shows each 
        genomovar separately and excludes singletons, E compares each reference
        genome to all genomovars and singletons in the same phylogroup added up
        at once.
 
C) each reference genome is compared one to many genomes with genomes 
        outside of its own phylogroup, but only one phylogroup at a time. This
        includes the genomes from singleton genomovars.
 
###
 
For F and D, first I collect the pairwise data for each reference genome to all
other genomes that are NOT in the same phylogroup. Then I split these pairwise
data for each reference genome into the same species or different species.
 
F) each reference genome is compared one to many genomes with genomes
        outside of its own phylogroup, but to genomes of all other phylogroups
        at the same time. This includes the genomes from singleton genomovars.
        Whereas C shows each phylogroup separately, F compares each reference
        genome to genomes from all phylogroups in the species added up at once.
 
D) each reference genome is compared one to many genomes with genomes
        outside of the species – all genomes outside the species at the same
        time.


* Excludes self matches.
* Skips genomes not in metadata file.

All vs. All RBM only includes RBM genes of a genome and not all genes.
Gene list file has all genes for each genome and is used get the total
gene count for each genome used as the denominator for the recombining
gene fraction.

Genes of column 1 in the RBM file are always in order.
Appending them in a list with preserve the gene order.

Input files:

    - All vs. All RBM file
    - Full gene list - complete_gene_list.txt
    - Metadata tsv file with columns Genome, Genomovar, Phylogroup, Species
      * NOTE: Column order of Metadata tsv file is CRITICAL

    * to create complete_gene_list.txt:
      'grep '>' 03a_all_CDS.fnn > 08_complete_gene_list.txt'

Output files:

    - tsv data files
    - Vectorized PDF hybrid violin plot

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

import argparse, sys, string, re
from collections import defaultdict
from copy import copy
import numpy as np
import pandas as pd; import numpy as np
import matplotlib.pyplot as plt


def parse_metadata_file(metadata):
    # reads metadata into md dict of {genome name: phylogroup}
    # return md dict
    md = {}

    # track population structure
    ps = {'gv': {}, 'pg': {}, 'sp': {}}

    with open(metadata, 'r') as file:
        header = file.readline()
        for line in file:
            X = line.rstrip().split('\t')
            genome = X[0]
            genome = re.sub(r"[.-]", "_", genome)
            genomovar = X[1]
            genomovar = re.sub(r"[.-]", "_", genomovar)
            phylogroup = X[2]
            phylogroup = re.sub(r"[.-]", "_", phylogroup)
            species = X[3]
            species = re.sub(r"[.-]", "_", species)
            ps['gv'][genomovar] = ''
            ps['pg'][phylogroup] = ''
            ps['sp'][species] = ''
            # store metadata for each genome
            md[genome] = [genomovar, phylogroup, species]

    # get population structure label
    if len(ps['gv']) == 1:
        print(
            '\n\nERROR: This script is intended for use with datasets '
            'consisting of at least TWO or more genomovars.\nERROR: '
            'For single genomovar analysis please consider using this '
            'script: 03i_RBM-Clade_Rarefaction.py\n\n'
            )
        sys.exit(1)
    # if pg == 1 then sp will also always == 1
    # can't have two species in the same phylogroup.
    elif len(ps['pg']) == 1:
        pop_struc = 'gv'
    elif len(ps['sp']) == 1:
        pop_struc = 'gv_pg'
    # if the code passes the elif statements above we have all three.
    else:
        pop_struc = 'gv_pg_sp'

    return md, pop_struc


def parse_gene_list(genelist):
    # the RBM file only has RBM genes not all genes
    # need total gene count for denominator of recombining gene fraction
    # create gene position list of 0's to add RBM gene positions to
    
    contigs = defaultdict(list)
    pnct = {} # {contig_name: [0,0]} gene array of zeros
    ctgct =  defaultdict(list)# {genome_name: [contig_name_array]}

    with open(genelist, 'r') as file:
        header = file.readline()
        for line in file:
            X = line[1:].split(' ')[0]
            X = re.sub(r"[.-]", "_", X)
            x = X.split('_')
            genome = '_'.join(x[:-2])
            contig = '_'.join(x[:-1])
            gene = int(x[-1])
            contigs[contig].append(gene)
            ctgct[genome].append(contig)

    for contig, genes in contigs.items():
        gene_array = [0 for i in range(max(genes))]
        pnct[contig] = gene_array
        #g1 = '_'.join(contig.split('_')[:-1])

    for k, v in ctgct.items():
        ctgct[k] = set(v)

    return pnct, ctgct


def parse_gpairs(rbm, pnct, ctgct, recx):

    # create entry for each gpair and store binary gene list gnlst
    # 1 - RBM is above rec threshold, 0 - not.
    # recx is the sequence identity threshold that can be set by the user
    # to determine when an RBM should be considered as recombinant or not
    data1a = {}
    data1b = defaultdict(list)
    # returns data1b {gpair: gene list of [1,0]}

    with open(rbm, 'r') as file:

        for i, line in enumerate(file, start=1):
            if i % 1000000 == 0:
                sys.stdout.write(f"\r\t\tLines processed: {i}")
                sys.stdout.flush()
            X = line.rstrip().split()
            g1, g2 = X[0], X[1]
            g1, g2 = re.sub(r"[.-]", "_", g1), re.sub(r"[.-]", "_", g2)
            x = g1.split('_')
            g1 = '_'.join(x[:-2])
            contig = '_'.join(x[:-1])
            gene = int(x[-1]) - 1 # python zero indexed
            g2 = '_'.join(g2.split('_')[:-2])
            # skip self matches
            if g1 == g2: continue
            gpair = f'{g1}-{g2}'
            pid = float(X[2])
            # no need to change anything unless RBM is recombining
            #if pid != 100: continue
            score = 1 if pid >= recx else 0
            # add gpair to dict
            if gpair not in data1a:
                data1a[gpair] = {}
            # add contig to gpair dict
            if contig not in data1a[gpair]:
                data1a[gpair][contig] = copy(pnct[contig])
            # update gene value in data1a[pgair][contig]
            data1a[gpair][contig][gene] = score
    
    # iterate through data1a and merge contigs
    for gpair, contigs in data1a.items():
        g1 = gpair.split('-')[0]
        # get full contig length zero arrays for all contigs in genome
        all_contigs = {}
        for i in ctgct[g1]:
            c = copy(pnct[i])
            all_contigs[i] = c
        for contig, genes in contigs.items():
            all_contigs[contig] = copy(genes)
        # join all contigs and add to data1b
        joined_contigs = [i for c in all_contigs.values() for i in c]
        data1b[gpair] = copy(joined_contigs)

    return data1b


def sort_gpairs(data1, md):

    # create dict for each genome name by genome pair assignment
    # store gnlst so each nested dict contains a list of gnlst 
    # for each genome in each category.
    gvd = defaultdict(list) # genomovar
    pgd = defaultdict(list) # phylogroup
    spd = defaultdict(list) # species

    for gpair, gnlst in data1.items():
        gs = copy(gpair.split('-'))
        g1, g2 = gs[0], gs[1]

        # get metadata if exists
        if g1 in md and g2 in md:
            md1, md2 = copy(md[g1]), copy(md[g2]) # get metadata
            gv1, gv2 = md1[0], md2[0] # get genomovar assignments
            pg1, pg2 = md1[1], md2[1] # get phylogroup assignments
            sp1, sp2 = md1[2], md2[2] # get species

        else:
            print(f'\t\tSkipping gpair: One or both not in metadata: {g1}, {g2}')
            continue

        # alphabetical sort group names to combine g1-g2 and g2-g1 combos
        gvpair = '-'.join(sorted([gv1, gv2]))
        pgpair = '-'.join(sorted([pg1, pg2]))
        sppair = '-'.join(sorted([sp1, sp2]))

        # compare genomovars in the same phylogroup.
        if pg1 == pg2:
            gvd[f'{g1}::{gvpair}'].append(copy(gnlst))
        # compare same phylogroup outside of genomovar
        if gv1 != gv2:
            pgd[f'{g1}::{pgpair}'].append(copy(gnlst))
        # compare same species outside of phylogroup
        if pg1 != pg2:
            spd[f'{g1}::{sppair}'].append(copy(gnlst))

        # gvd contains binary gene lists for all genome pairs in the same
        # phylogroup that are in a a genomovar with more than 1 genome.

        # pgd contains scores for all genomes that are not in the same
        # genomovar

        # spd contains scores for all genomes that are not in the same
        # phylogroup

    return gvd, pgd, spd


def score_gpairs(list_dict, outfile):

    score_dict = defaultdict(list)

    for label, scores in list_dict.items():
        # label is one of:
        # f'{g1}::{gv1}-{gv2}'
        # f'{g1}::{pg1}-{pg2}'
        # f'{g1}::{sp1}-{sp2}'

        X = label.split('::')
        genome = X[0]
        grouping = X[1]
        dt = np.array([np.array(copy(score)) for score in scores])
        sm = dt.sum(axis=0)
        ty = [0 if i == 0 else 1 for i in sm]
        recombining_genes = sum(ty)
        total_genes = len(ty)
        recombining_fraction = recombining_genes / total_genes

        score_dict[grouping].append(recombining_fraction)

        outfile.write(f'{genome}\t{grouping}\t{recombining_fraction}\n')

    return score_dict


def sort_scores(score_dict, gpmin):

    same_list, diff_list = [], []
    same_dict, diff_dict = {}, {}

    for group, scores in score_dict.items():
        # skip group if fewer than gpmin pairs.
        if len(scores) < gpmin: continue
        # get group pairing
        X = group.split('-')
        g1, g2 = X[0], X[1]

        if g1 == g2:
            same_list.extend(scores)
            same_dict[group] = scores
        elif g1 != g2:
            diff_list.extend(scores)
            diff_dict[group] = scores
    
    return same_list, diff_list, same_dict, diff_dict


def score_genomes_by_category(gvd, pgd, spd, gpmin, outpre, pop_struc):

    # data3 stores recombining fraction values for each gv, pg, sp grouping
    header = "Genome\tPairing\tRecombining Genome Fraction\n"

    # gvd contains scores for all genome in the same phylogroup that are
    # in a a genomovar with more than 1 genome.
    gvdout = open(f'{outpre}_genomovar_data.tsv', 'w')
    gvdout.write(header)
    gvd = score_gpairs(gvd, gvdout) # genomovar
    gvdout.close()

    # pgd contains scores for all genomes that are not in the same
    # genomovar
    pgdout = open(f'{outpre}_phylogroup_data.tsv', 'w')
    pgdout.write(header)
    pgd = score_gpairs(pgd, pgdout) # phylogroup
    pgdout.close()

    # spd contains scores for all genomes that are not in the same
    # phylogroup
    spdout = open(f'{outpre}_species_data.tsv', 'w')
    spdout.write(header)
    spd = score_gpairs(spd, spdout) # species
    spdout.close()

    # create 6 category dict with lists
    # each value in the "sort_scores" lists is the fraction of recombining genes
    # in 1 genome, to the rest of genomes/genes in that category.
    
    sgvl, dgvl, sgvd, dgvd = sort_scores(gvd, gpmin)
    spgl, dpgl, spgd, dpgd = sort_scores(pgd, gpmin)
    sspl, dspl, sspd, dspd = sort_scores(spd, gpmin)

    if pop_struc == "gv_pg_sp":
        summary = {
                'A': sgvl, # same genomvar
                'B': dgvl, # different genomvar
                'C': dpgl, # different phylogroup
                'D': dspl, # different species
                'E': spgl, # same phylogroup
                'F': sspl, # same species
                }

    elif pop_struc == "gv_pg":
        summary = {
                'A': sgvl, # same genomvar
                'B': dgvl, # different genomvar
                'C': dpgl, # different phylogroup
                'D': spgl, # same phylogroup
                'E': sspl, # same species
                }

    elif pop_struc == "gv":
        summary = {
                'A': sgvl, # same genomvar
                'B': dgvl, # different genomvar
                'C': spgl, # same phylogroup
                }

    # create output list of data
    # note spd for species level is not split into same and different.
    data3 = [sgvd, dgvd, spgd, dpgd, spd, summary]

    return data3


def parse_rbm_file(rbm, md, pnct, ctgct, gpmin, outpre, pop_struc, recx):
    # parses rbm file, matches metadata
    # returns data dict of {}

    # store each gene pairs list of gene scores
    # gene score: 1 = Recombining; 0 = not Recombining
    data1 = parse_gpairs(rbm, pnct, ctgct, recx)

    # sort genome pairs into md categories and store gnlst
    # by g1 name. Each nested dict contains list of gnlst
    gvd, pgd, spd = sort_gpairs(data1, md)

    # for each genome in each category make a 2d numpy array of gnlst
    # sum numpy array by column. Any element ≥ 1 recombines within
    # that category and any 0 is a gene that does not.
    # sum the 1d numpy column sum and divide by length gives fraction
    # of genes in the genome that are recombining with other genes in 
    # other genomes in that category. Store that result by category
    # to plot in the boxplots
    data3 = score_genomes_by_category(gvd, pgd, spd, gpmin, outpre, pop_struc)

    return data3


def get_lw_xlabsize(ncat):

    if ncat <= 4:
        lw, xlabsize = 10, 12
    elif ncat <=8:
        lw, xlabsize = 8, 12
    elif ncat <=12:
        lw, xlabsize = 6, 10
    elif ncat <=16:
        lw, xlabsize = 4, 8
    elif ncat <= 24:
        lw, xlabsize = 2, 6
    elif ncat > 24:
        lw, xlabsize = 1, 4

    return lw, xlabsize


def build_violin_plot(data, title, outfile):
    # simple bar plot for data

    labels = list(data.keys())
    values = list(data.values())

    x, q1s, mds, q3s = [], [], [], []

    for i, d in enumerate(values):
        x.append(i+1)
        if len(d) >= 3:
            q1, md, q3 = np.percentile(d, [25, 50, 75])
        elif len(d) > 0:
            print(f'\n\n\t\tCAUTION: {labels[i]} has < 3 data points.\n\n')
            q1, md, q3 = np.mean(d), np.mean(d), np.mean(d)
        else:
            print(f'\n\n\t\tCAUTION: {labels[i]} has 0 data points.\n\n')
            q1, md, q3 = 0, 0, 0
            values[i].append(0)

        q1s.append(q1)
        mds.append(md)
        q3s.append(q3)

    if sum(mds) == 0:
        print((f'\n\n\t\tERROR: Not enough genomes per group to compare.\n\n'))

    # add inter quartile range boxplot bar
    # scale line width with number of x-axis categories
    ncat = len(labels) # number of x-axis categories
    lw, xlabsize = get_lw_xlabsize(ncat)

    fig, ax = plt.subplots(figsize=(7,5))

    parts = ax.violinplot(values, showmeans=False, showmedians=True)
    # add inter quartile range boxplot bar    
    ax.vlines(x, q1s, q3s, color='k', linestyle='-', lw=lw)
    

    # set plot title
    ax.set_title(f'Cumulative recombining RBM gene fractions by {title}')
    # change axis labels
    ax.set_xlabel('')
    ax.set_ylabel("Cumulative recombining gene fraction", fontsize=12)

    # set style for the axes
    ax.set_xticks(np.arange(1, len(labels) + 1), labels=labels)
    #ax.set_xlim(0.25, len(labels) + 0.75)

    # set the axis parameters / style
    ax.minorticks_on()
    ax.tick_params(axis='both', labelsize=xlabsize)
    ax.tick_params(axis='x', which='minor', bottom=False)
    # set grid style
    ax.yaxis.grid(
        which="minor", color='#f0f0f0', linestyle='--', linewidth=.75
        )
    ax.yaxis.grid(
        which="major", color='#d9d9d9', linestyle='--', linewidth=1
        )
    ax.set_axisbelow(True)

    if title == 'Summary':
        ax.tick_params(axis='x', labelrotation=0)
        plt.subplots_adjust(left=0.12, right=0.98, top=0.95, bottom=0.24)

        if ncat == 6:
            txt1 = (
                'A - same genomovar\n'
                'B - different genomovar\n'
                'C - different phylogroup\n'
                'D - different species'
                )

            txt2 = (
                'E - same phylogroup, different genomovars\n'
                'F - same species, different phylogroups'
                )

            ax.vlines(x=4.5, ymin=0, ymax=1, ls='--', lw=2, color='k')

        elif ncat == 5:
            txt1 = (
                'A - same genomovar\n'
                'B - different genomovar\n'
                'C - different phylogroup'
                )

            txt2 = (
                'D - same phylogroup, different genomovars\n'
                'E - same species, different phylogroups'
                )

            ax.vlines(x=3.5, ymin=0, ymax=1, ls='--', lw=2, color='k')

        elif ncat == 3:
            txt1 = (
                'A - same genomovar\n'
                'B - different genomovar'
                )

            txt2 = (
                'C - same phylogroup, different genomovars'
                )

            ax.vlines(x=2.5, ymin=0, ymax=1, ls='--', lw=2, color='k')

        ax.annotate(
                    txt1, # text to annotate
                    xy=(0.25, 0), # reference point for placement
                    xycoords=('figure fraction', 'figure fraction'),
                    xytext=(-100, 5), # distance from reference point
                    textcoords='offset points',
                    size=12, ha='left', va='bottom',
                    annotation_clip=False
                    )
        ax.annotate(
                    txt2, # text to annotate
                    xy=(0.65, 0), # reference point for placement
                    xycoords=('figure fraction', 'figure fraction'),
                    xytext=(-100, 5), # distance from reference point
                    textcoords='offset points',
                    size=12, ha='left', va='bottom',
                    annotation_clip=False
                    )
    else:
        ax.tick_params(axis='x', labelrotation=90)
        fig.set_tight_layout(True)
    
    # savefig and close
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
        '-gl', '--full_gene_list',
        help='Please specify the pancat file!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-rbm', '--allv_rbm_file',
        help='Please specify All vs. All RBM file!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-gpmin', '--minimum_genome_pairs',
        help='(OPTIONAL) Minimum genome pairs to include group (default: 2)!',
        metavar='',
        type=int,
        required=False,
        default=2
        )
    parser.add_argument(
        '-recx', '--Recombinat_Sequence_Identity_Threshold',
        help='(OPTIONAL) Set threshold to determine recombinant RBMs (default: 99.8)!',
        metavar='',
        type=int,
        required=False,
        default=2
        )
    parser.add_argument(
        '-o', '--output_file_prefix',
        help='Please specify an output file prefix!',
        metavar='',
        type=str,
        required=True
        )
    args=vars(parser.parse_args())

    ####################################################################

    # Do what you came here to do:
    print('\n\nRunning Script...')

    # define parameters
    metadata = args['metadata_file']
    genelist = args['full_gene_list']
    rbm = args['allv_rbm_file']
    gpmin = args['minimum_genome_pairs']
    recx = args['Recombinat_Sequence_Identity_Threshold']
    outpre = args['output_file_prefix']

    ####################################################################

    # parse metadata
    print(f'\n\nParsing metadata ...')
    md, pop_struc = parse_metadata_file(metadata)

    ####################################################################

    # parse gene list file
    # pnct is a dict with an array of zeros for gene positions for each contig
    # ctgct is a dict with a list of contig names for each genome
    print(f'\n\nParsing gene list ...')
    pnct, ctgct = parse_gene_list(genelist)

    ####################################################################

    # parse rbm file
    print(f'\n\nParsing RBM file ...')
    data3 = parse_rbm_file(rbm, md, pnct, ctgct, gpmin, outpre, pop_struc, recx)

    ####################################################################
    ####################################################################

    # plotting data
    # iterate data3 dict and plot each group of data.
    print(f'\n\nPlotting data ...')
    switch = {
                0: 'Same Genomovar', 1: 'Different Genomovar',
                2: 'Same Phylogroup', 3: 'Different Phylogroup',
                4: 'Species', 5: 'Summary'
                }
    for i, d in enumerate(data3):
        # Skip 3: 'different phylogorup' or 4: 'species' if empty
        if len(d) == 0: continue

        # Set the outfile name for the label matching the data
        olab = '_'.join(switch[i].split(' '))
        
        # Set the plot title matching the data
        title = switch[i]

        # if more than 40 category pairs on x-axis create additional plots.
        if len(d) > 40:
            labels = list(d.keys())
            letters = list(string.ascii_lowercase)
            abc = []
            for i in range(1,10):
                abc = abc + [A * i for A in letters]
            j = 0
            while len(labels) > 40:
                ilabs = labels[:40]
                nd = {k: d[k] for k in ilabs}
                outfile = f'{outpre}-figure0{i+1}{abc[j]}-{olab}.pdf'
                nt = f'{title} - {abc[j]}'
                _ = build_violin_plot(nd, nt, outfile)
                labels = labels[40:]
                j += 1
            # plot remaining
            nd = {k: d[k] for k in labels}
            outfile = f'{outpre}-figure0{i+1}{abc[j]}-{olab}.pdf'
            nt = f'{title} - {abc[j]}'
            _ = build_violin_plot(nd, nt, outfile)

        else:
            outfile = f'{outpre}-figure0{i+1}-{olab}.pdf'
            _ = build_violin_plot(d, title, outfile)


    print(f'\n\nComplete success space cadet!! Finished without errors.\n\n')

if __name__ == "__main__":
    main()


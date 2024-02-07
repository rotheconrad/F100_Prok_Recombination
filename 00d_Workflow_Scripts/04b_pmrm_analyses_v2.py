 #!/usr/bin/env python

'''Purged mutations vs. recent recombination events from pairwise RBM data

RBM data = pairwise all vs. all reciprocal best match genes.

To estimate more precisely the relative contribution of recombination compared
to diversifying mutation, we developed an empirical approach based on the 
sequence identity patterns across the genome. Our approach identifies 
recombined genes that represent recent events, i.e., showing 99.8-100% 
nucleotide sequence identity, and subsequently calculates how much sequence 
divergence these events presumably removed based on the ANI of the genomes 
compared. For example, if the ANI of the two genomes compared is 97%, this 
would mean that the divergence of the recombined genes was about 3% (+/-0.2%), 
on average, before the recombination took place, and thus recombination should 
have removed (purged) a total of nucleotide differences that should roughly be 
equal to [3% X total length of recombined genes]. During the same evolutionary 
time, diversifying mutation can create, at maximum, a total of nucleotide 
differences that should roughly be equal to [0.02 X total genome length] (for 
focusing on recent evolution only that corresponds to 99.8-100% identity events
 or 0.02% maximum sequence divergence accumulated). Using this approach, we 
 empirically estimated the magnitude of recombination (ρ) vs. point mutation 
 (θ) to be about 5, on average, for genome pairs of different genomovars.

This is defined as:

purged mutations (pm) = (100 - rbm_ani) / 100 * rec_length
recent mutations (rm) = divergence_time / 100 * total length

RBM = reciprocal best match gene
Rbm_ani is the ANI calculated from the average of RBM sequence identities

RBM genes ≥ 99.8% are classified as recent recombination events
when the default divergence_time is set to 0.2% sequence identity.

* default divergence_time is 0.2 but can be adjusted by the user with -dt

rec_length is the length of RBM genes identified as recent recombination.
total_length is the total length of all RBM genes between the genome pair.

If the ratio > 1 recombination removes more mutations
If the ratio < 1 evolution produces more mutations

###

I. Pairwise all vs. all data are written to a tab separated value file:
    _pmpr_pairwise_data.tsv. this data can be plotted against ANI or F100
    using the script: plot_allv_pmrm_data.py

II. Cumulative one vs. many results are plotted for each genome against the
    following different groupings of "many" genomes.

    These cumulative data are written to tsv files and hybrid violin plots
    in publication ready vectorized PDFs are also created.

1) we store the RBM sequence identity for each genome pair
2) we compute the ANI from the RBM genes of each genome pair
2) we store the gene lengths for each genome pair preserving gene position
3) we store the total genome length for each genome (based on RBMs gene length)
5) we compute ρ/θ (pm/rm) for each genome pair
6) we compute ρ/θ  (pm/rm) for each reference genome vs many grouping
   for this we use the average ANI of the many group to the reference

### 

One vs. many - "many" groupings are as follows:

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

###

* A note on the required gene list input:
    All vs. All RBM only includes RBM genes of a genome and not all genes.
    Gene list file has all genes for each genome and is used get the total
    gene count for each genome used as the denominator for the identical
    gene fraction.

    Genes of column 1 in the RBM file are always in order.
    Appending them in a list with preserve the gene order.

####

Input files:

    - All vs. All RBM file
    - Full gene list - complete_gene_list.txt
    - Metadata tsv file with columns Genome, Genomovar, phylogroup, species

Output files:

    - Pairwise all vs. all tsv file
    - cumulative one vs. many tsv files
    - Vectorized PDFs - hybrid violin plots

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: Jan 2024
License :: GNU GPLv3
Copyright 2024 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse
from collections import defaultdict
from copy import copy
import numpy as np
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
            genomovar = X[2]
            species = X[3]
            md[genome] = [phylogroup, genomovar, species]

    return md


def parse_gene_list(genelist):
    # the RBM file only has RBM genes not all genes
    # need total gene count for denominator of identical gene fraction
    # create gene position list of 0's to add RBM gene positions to
    
    contigs = defaultdict(list)
    pnct = {} # {contig_name: [0,0]} gene array of zeros
    ctgct =  defaultdict(list)# {genome_name: [contig_name_array]}

    with open(genelist, 'r') as file:
        header = file.readline()
        for line in file:
            X = line[1:].split(' ')[0]
            x = X.split('_')
            genome = '_'.join(x[:-2])
            contig = '_'.join(x[:-1])
            gene = int(x[-1])
            contigs[contig].append(gene)
            ctgct[genome].append(contig)

    for contig, genes in contigs.items():
        gene_array = [0 for i in range(max(genes))]
        pnct[contig] = gene_array

    for k, v in ctgct.items():
        ctgct[k] = set(v)

    return pnct, ctgct


def parse_rbm_file(rbm, md, pnct, ctgct):
    # parses rbm file, matches metadata
    # returns data dict of {gene pair: }

    d1 = {} # stores sequence identity (pid) and gene length (glen)
    d2 = defaultdict(list) # store pid to compute RBM ANI for each gpair
    dGL = defaultdict(int) # store glen to compute rbm gene base genome length
    d3 = defaultdict(list) # concatenated contig [pid, glen] by gpair

    # read in the RBM file and store pid and glen by gene position
    with open(rbm, 'r') as file:
        for line in file:
            X = line.rstrip().split('\t')
            x = X[0].split('_')
            g1 = '_'.join(x[:-2])
            contig = '_'.join(x[:-1])
            gene = int(x[-1]) - 1 # python zero indexed
            g2 = '_'.join(X[1].split('_')[:-2])
            # skip self matches
            if g1 == g2: continue
            gpair = f'{g1}-{g2}'
            pid = float(X[2])
            glen = max(int(X[12]), int(X[13])) # use longest gene length
            
            # add gpair to dict
            if gpair not in d1:
                d1[gpair] = {}
            # add contig to gpair dict
            if contig not in d1[gpair]:
                d1[gpair][contig] = copy(pnct[contig])
            # update gene value in data1a[pgair][contig]
            d1[gpair][contig][gene] = [pid, glen]

            # store pid for rbm ani
            d2[gpair].append(pid)

            # store glen for genome length total
            dGL[gpair] += glen

    # iterate through data1a and merge contigs
    for gpair, contigs in d1.items():
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
        d3[gpair] = copy(joined_contigs)

    # compute RBM ANI for each genome pair
    # the F100 score includes tied rbms while dict d3 will not
    dANI = {}
    dF100 = {}
    for gpair, ANIs in d2.items():
        ANI = sum(ANIs) / len(ANIs)
        dANI[gpair] = ANI
        rbms = len(ANIs)
        r100s = len([i for i in ANIs if i >= 100])
        F100 = r100s / rbms
        dF100[gpair] = F100

    return d3, dGL, dANI, dF100


def compute_allv_pmrm(d3, dGL, dANI, dF100, dt, outpre):

    # compute all vs all ρ/θ (pairwise each genome pair)
    # pmrm is th pm /rm or ρ/θ ratio
    # rlen is the rec_length or length of genes with pid ≥ divergence_time
    # default divergence_time is 0.2 coinciding with 99.8% ANI
    # tlen is totel length of RBM genes shared by the genome pair.
    
    allv = {'gpair': [], 'ρ/θ': [], 'ani': [], 'f100': [], 'rlen': [], 'tlen': [], 'rec_ani': []}
    for gpair, values in d3.items():
        rbm_ani = dANI[gpair]
        f100 = dF100[gpair]
        total_length = dGL[gpair]
        rec_length = 0
        rec_rbms = []
        for i in values:
            if i != 0:
                pid, glen = i[0], i[1]
                if pid >= 100 - dt:
                    rec_length += glen
                    rec_rbms.append(pid)

        # pm / rm - this is the goal
        pm_rate = (100 - rbm_ani) / 100 # remove the %
        pm = pm_rate * rec_length

        #rm_rate = dt / 100 # divergence time divided by 100 to remove the %
        rec_ani = sum(rec_rbms) / len(rec_rbms)
        rm_rate = (100 - rec_ani) / 100 

        rm = rm_rate * total_length
        pmrm = pm / rm if rm > 0 else 0

        # add everything to the allv dict
        allv['gpair'].append(gpair)
        allv['ρ/θ'].append(pmrm)
        allv['ani'].append(rbm_ani)
        allv['f100'].append(f100)
        allv['rlen'].append(rec_length)
        allv['tlen'].append(total_length)
        allv['rec_ani'].append(rec_ani)

    df = pd.DataFrame(allv)
    outfile = f'{outpre}_pmpr_pairwise_data.tsv'
    df.to_csv(outfile, sep='\t', index=False)

    print(f'\n\nWriting pairwise all vs. all data to: {outfile}:')

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
        '-o', '--output_file_prefix',
        help='Please specify an output file prefix!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-dt', '--divergence_time',
        help='(OPTIONAL) ANI divergence (default = 0.2)!',
        metavar='',
        type=float,
        required=False,
        default=0.2
        )
    parser.add_argument(
        '-gpmin', '--minimum_genome_pairs',
        help='(OPTIONAL) Minimum genome pairs to include group (default: 3)!',
        metavar='',
        type=int,
        required=False,
        default=3
        )
    parser.add_argument(
        '-lg', '--log_yaxis',
        help='OPTIONAL: Set -lg True to plot ln(yaxis) (Default=None).',
        metavar='',
        type=str,
        default=None,
        required=False
        )
    args=vars(parser.parse_args())

    # Do what you came here to do:
    print('\n\nRunning Script...')

    # define parameters
    metadata = args['metadata_file']
    genelist = args['full_gene_list']
    rbm = args['allv_rbm_file']
    outpre = args['output_file_prefix']
    dt = args['divergence_time']
    gpmin = args['minimum_genome_pairs']
    lg = args['log_yaxis']

    # parse metadata
    print(f'\n\nParsing metadata ...')
    md = parse_metadata_file(metadata)
    
    # parse gene list file
    # pnct is a dict with an array of zeros for gene positions for each contig
    # ctgct is a dict with a list of contig names for each genome
    print(f'\n\nParsing gene list ...')
    pnct, ctgct = parse_gene_list(genelist)

    # parse rbm file
    print(f'\n\nParsing RBM file ...')
    # store each gene pairs list of gene scores
    # d3 = {gpair: [RBM pID, gene_length]}
    d3, dGL, dANI, dF100 = parse_rbm_file(rbm, md, pnct, ctgct)

    # pairwise computation
    print(f'\n\nComputing all vs. all pairwise pm/rm ...')
    _ = compute_allv_pmrm(d3, dGL, dANI, dF100, dt, outpre)

    print(f'\n\nComplete success space cadet!! Finished without errors.\n\n')

if __name__ == "__main__":
    main()

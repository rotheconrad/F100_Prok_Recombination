 #!/usr/bin/env python

'''Process and random sample genome pairs 03c_AAI_RBM_F100.py output

This script is used for a directory with output files for  many species.
It reads all files in the directory and stores results with the file name
(Escherichia_coli_F100.tsv). From there it randomly samples n genome
pairs from each species with more than n genome pairs. It samples all
genomes pairs from species with fewer than n genome pairs.

Performs x number of random experiments and writes out a file for each
experiment. Results can be plotted with 03d_AAI_RBM_scatter.py.

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

import argparse, os, random
from pathlib import Path
from collections import defaultdict

random.seed(42)

"""
# Set the default sans-serif font to Helvetica
matplotlib.rcParams['font.sans-serif'] = "Helvetica"
# Set default font to sans-serif
matplotlib.rcParams['font.family'] = "sans-serif"
"""


def random_sample(data, n):
    """ select n random samples from each key: list in dict """

    sample = []

    for species, gpairs in data.items():
        samp = random.sample(gpairs, n)
        sample.extend(samp)

    return sample


def sample_rbm_directory(infile, outpre, n, x):
    """ Reads the file does the work """

    #initialize dictionAAary to store data
    data = defaultdict(list)

    # get the list of files from the genome directory
    flist = [f for f in os.listdir(inDir) if os.path.isfile(f'{inDir}{f}')]

    # read through files
    for file in flist:
        with open(file, 'r') as ifile:
            species = file.split('.')[0]
            for line in ifile:
                newline = f'{line.rstrip()}\t{species}\n'
                data[species].append(newline)

    for i in range(x):
        sample = random_sample(data, n)
        with open(f'{outpre}_experiment_{n}.tsv', 'o') as outfile:
            for line in sample: outfile.write(line)

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
        '-o', '--output_file_prefix',
        help='Please specify the output file prefix?',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-n', '--select_n_genome_pairs',
        help='OPTIONAL: Select n genome pairs (Default=10)',
        metavar='',
        type=int,
        default=10,
        required=False
        )
    parser.add_argument(
        '-x', '--repeat_experiment_x_times',
        help='OPTIONAL: Repeat random sampling x times (Default=10)',
        metavar='',
        type=int,
        default=10,
        required=False
        )
    args=vars(parser.parse_args())

    # Do what you came here to do:
    print('\n\nRunning Script...\n')

    # define parameters
    infile = args['input_file']
    outpre = args['output_file_prefix']
    n = args['select_n_genome_pairs']
    x = args['repeat_experiment_x_times']

    # run it
    _ = sample_rbm_directory(infile, outpre, n, x)



    print(f'\n\nComplete success space cadet!! Finished without errors.\n\n')

if __name__ == "__main__":
    main()

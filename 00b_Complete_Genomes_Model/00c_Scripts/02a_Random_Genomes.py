#!/usr/bin/env python

'''Random select n genomes within ANI range

Input is a directory with genomes in fasta format.
    1) Randomly selects SEED genome
    2) Randomly selects additional genomes and compares to SEED
    3) Keeps genome if within ANI range until n genomes are found.

Output is a new directory with sampled genomes.

Dependencies:
    fastANI in users PATHs.

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

import argparse, sys, os, random, subprocess
from pathlib import Path

def if_gzip(outDir, name):
    ''' unzips files that are gzipped '''
    _ = subprocess.run(['gunzip', f'{outDir}{name}'])
    name = '.'.join(name.split('.')[:-1])
    return name

def setup_the_bomb(n, inDir, outDir):
    ''' Selects SEED genome, creates new directory and list of files '''

    # get the list of files from the genome directory
    g_list = [g for g in os.listdir(inDir) if os.path.isfile(f'{inDir}{g}')]

    # make a new directory for this sampling
    _ = Path(outDir).mkdir(parents=True, exist_ok=True)

    # Select the first genome of this sample, copy, unzip, and rename it
    g_prime = random.choice(g_list) 
    _ = g_list.remove(g_prime)
    _ = subprocess.run(['cp', f'{inDir}{g_prime}', outDir])

    # if the fasta is gzipped unzip
    if g_prime.split('.')[-1] == 'gz': g_prime = if_gzip(outDir, g_prime)

    return g_list, g_prime

def select_genomes(n, inDir, outDir, ani_low, ani_high):

    # Select a SEED genome and setup output
    print('Selecting Random Genome 1 to use as SEED.')
    print('##################################################\n\n')
    g_list, g_prime = setup_the_bomb(n, inDir, outDir)
    prime = f'{outDir}{g_prime}'
    genomes_avail = len(g_list)+1

    # Keep track of genomes found and genomes tested. Fail cases.
    i, j = 1, 2
    while (j < n+1):
        if len(g_list) == 0:
            print('\n\n############# FAIL ####################################')
            print(
                f'Sampled all {genomes_avail} genomes in {inDir} directory and '
                f'only found {j} of the {n} matches requested between {atl}% '
                f'and {atu}% ANI!\n\n'
                )
            sys.exit()

        if i == 100000:
            print('\n\n############# FAIL ####################################')
            print(
                f'Sampled 100,000 Genomes and and only found {j} of the {n} '
                f'matches requested between {atl}% and {atu}% ANI\n\n'
                )
            sys.exit()

        # update i
        i += 1

        # Select and test genomes from the list.
        print(f'Testing Random Genome {i} of {genomes_avail}')
        g_new = random.choice(g_list)
        _ = g_list.remove(g_new)

        print('g_list len:', len(g_list))
        _ = subprocess.run(['cp', f'{inDir}{g_new}', outDir])
        if g_new.split('.')[-1] == 'gz': g_new = if_gzip(outDir, g_new)

        new = f'{outDir}{g_new}'
        temp = f'{outDir}temp.ani'
        fastANI = f"fastANI -r {prime} -q {new} -o {temp}"
        out = subprocess.Popen(fastANI, shell=True, stdout=subprocess.PIPE)
        out = out.communicate()[0]
        if os.stat(temp).st_size > 0:
            with open(temp, 'r') as f:
                ani = float(f.readline().split('\t')[2])
                print(ani)

            if ani < ani_low:
                _ = subprocess.run(['rm', new])
                print(f'FAIL: {ani}% ANI is below the {ani_low}% threshold!')

            elif ani > ani_high:
                _ = subprocess.run(['rm', new])
                print(f'FAIL: {ani}% ANI is above the {ani_high}% threshold!')

            elif ani >= ani_low and ani < ani_high:
                print(
                    f'PASS: {ani}% ANI is within the {ani_low}% to {ani_high}% '
                    f'threshold range!\nFound {j} genomes of the {n} requested.'
                    )
                j += 1

            else:
                print('FREAK OUT!!')

        else: subprocess.run(['rm', new])

    _ = subprocess.run(['rm', temp])

    print('\n\n############# SUCCESS ####################################')
    print(f'Sampled {j-1} of {n} genomes within {ani_low}% and {ani_high}% ANI')


def main():
    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-n', '--n_genomes',
        help='How many genomes would you like to select?',
        metavar='',
        type=int,
        required=True
        )
    parser.add_argument(
        '-i', '--input_directory',
        help='Please specify the genomes directory! (fasta or gz fasta)',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-o', '--output_directory',
        help='Please specify the output directory! (fasta or gz fasta)',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-low', '--ani_threshold_low',
        help='OPTIONAL: Lower ANI threshold (default=95.0)',
        metavar='',
        type=float,
        required=False,
        default=95.0
        )
    parser.add_argument(
        '-high', '--ani_threshold_high',
        help='OPTIONAL: Upper ANI threshold (default=100.0)',
        metavar='',
        type=float,
        required=False,
        default=100.0
        )
    args=vars(parser.parse_args())
    # Run this scripts main function
    print('\n\nRunning Script...\n\n')

    # define parameters
    n = args['n_genomes']
    inDir = args['input_directory']
    outDir = args['output_directory']
    ani_low = args['ani_threshold_low']
    ani_high = args['ani_threshold_high']

    # check directory name to include /
    if inDir[-1] != '/': inDir = inDir + '/'
    if outDir[-1] != '/': outDir = outDir + '/'

    # do the thing
    _ = select_genomes(n, inDir, outDir, ani_low, ani_high)

    print(f'\n\nComplete success space cadet!! Finished without errors.\n\n')

if __name__ == "__main__":
    main()

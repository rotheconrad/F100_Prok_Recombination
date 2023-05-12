 #!/usr/bin/env python

'''Generates file combinations to use for running a for loop with the
03f_Recombinant_pair_analysis.py script.

The genomes directory and the gene CDS directory should contain the exact same
files and filenames except for the file extentions.

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: July 2022
License :: GNU GPLv3
Copyright 2022 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse, os, itertools

def get_file_lists(indir):

    file_list = []

    for item in os.listdir(indir):
        if not item.startswith('.') and os.path.isfile(os.path.join(indir, item)):
            file_list.append(f'{indir}/{item}')

    return sorted(file_list)

def generate_file_combinations(genomes, genes, outfile):

    paired = list(zip(genomes, genes))

    combos = list(itertools.combinations(paired, 2))

    with open(outfile, 'w') as out:
        for (g1, c1), (g2, c2) in combos:
            line = f'{g1},{c1},{g2},{c2}\n'
            out.write(line)

    return True


###############################################################################
###############################################################################

def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-g', '--genomes_dir',
        help='Please specify the genomes fasta directory!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-c', '--genes_CDS_directory',
        help='Please specify the genes CDS fasta directory!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-o', '--output_file_name',
        help='Please specify the output file name (use .csv)!',
        metavar='',
        type=str,
        required=True
        )
    args=vars(parser.parse_args())

    # Do what you came here to do:
    print('\n\nRunning Script...\n')

    # define parameters
    genomes_dir = args['genomes_dir']
    genes_dir = args['genes_CDS_directory']
    outfile = args['output_file_name']

    genomes = get_file_lists(genomes_dir)
    genes = get_file_lists(genes_dir)
    _ = generate_file_combinations(genomes, genes, outfile)


    print(f'\n\nComplete success space cadet!! Finished without errors.\n\n')
    
if __name__ == "__main__":
    main()

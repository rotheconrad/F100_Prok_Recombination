 #!/usr/bin/env python

'''Match filename to metadata for REHAB Ecoli genomes

Simple quick script I wrote to match my filenames with the supplemental
metadata file.

Writes three output files:

    - {outpre}_GroupFiles.txt
            Contains file path pairs for *.fna and *.fnn needed to create the
            group input files for 03g_Recombinant_group_analysis.py
    - {outpre}_Metadata.txt
            Contains each genomes file name and associated meta data (tsv).
    - {outpre}_Metacolors.txt
            Contains each unique metadata value. Add a color for each.

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: June 2023
License :: GNU GPLv3
Copyright 2023 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse

def get_metadata(metadata, filenames, outpre):

    mdata = {}
    mvalues = {}

    with open(metadata, 'r') as file:
        header = file.readline()
        for line in file:
            X = line.rstrip().split('\t')
            accession = X[4]
            site = X[10]
            timepoint = X[7]
            phylogroup = X[19]
            if X[8][0] == 'W':
                niche = X[9]
            else:
                niche = X[8]
            mdata[accession] = [site, niche, timepoint, phylogroup]
            mvalues[site] = ''
            mvalues[niche] = ''
            mvalues[timepoint] = ''
            mvalues[phylogroup] = ''

    out1 = open(f'{outpre}_GroupFiles.txt', 'w')
    out2 = open(f'{outpre}_Metadata.txt', 'w')
    out2.write('Accession\tSite\tNiche\tTime Point\tPhylogroup\n')
    with open(filenames, 'r') as file:
        for line in file:
            X = line.rstrip().split('\t')
            accession = '_'.join(X[0].split('_')[:2])
            if accession in mdata:
                name = X[1].split('.')[0]
                genome = X[2]
                genes = X[3]
                out1.write(f'{genome}\t{genes}\n')
                out2.write(name + '\t' + '\t'.join(mdata[accession]) + '\n')

    out3 = open(f'{outpre}_MetaColors.txt', 'w')
    for mval, _ in mvalues.items():
        out3.write(f'{mval},\n')

    out1.close(), out2.close(), out3.close()

    return True

def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-f', '--filenames_file',
        help='Please specify the files names file!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-m', '--metadata_file',
        help='Please specify the meta data file!',
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

    # define input parameters
    filenames = args['filenames_file']
    metadata = args['metadata_file']
    outpre = args['output_file_prefix']

    print(f'\n\nRunning script ...')
    
    _ = get_metadata(metadata, filenames, outpre)

    print(f'\n\nComplete success space cadet!! Finished without errors.\n\n')
    
if __name__ == "__main__":
    main()


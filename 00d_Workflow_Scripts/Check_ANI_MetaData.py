 #!/usr/bin/env python

'''Check ANI and Meta Data

In the initial analysis we noticed several clusters of 100% identical genomes
in the REHAB E. coli collection. However, the assigned phylogroups were not
the same. 100% ANI genomes should have the same phylogroup.

This script takes the meta data file and the ANI file and returns a tsv of the
100% ANI genome pairs with their metadata to assist in investigation.

-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: August 2023
License :: GNU GPLv3
Copyright 2023 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse

def parse_metadata(metadata):

    data = {}

    with open(metadata, 'r') as file:
        header = file.readline().rstrip()
        for line in file:
            X = line.rstrip().split('\t')
            Accession = X[0]
            Site = X[1]
            Niche = X[2]
            TimePoint = X[3]
            Phylogroup = X[4]

            data[Accession] = [Site, Niche, TimePoint, Phylogroup]

    return data


def parse_ani(anifile, data, outfile):

    tracker = {}

    with open(anifile, 'r') as file:

        for line in file:
            X = line.rstrip().split('\t')
            gA = X[0].split('/')[-1].split('.')[0]
            gB = X[1].split('/')[-1].split('.')[0]

            if gA == gB: continue # remove self match

            if gA in data and gB in data:

                gS = '-'.join(sorted([gA,gB]))
                rANI = round(float(X[2]), 2)
                ANI = float(X[2])
                tfrags = int(X[4])

                if gS in tracker:
                    old_tfrags = tracker[gS][4]
                    if tfrags > old_tfrags:
                        tracker[gS] = [gA, gB, rANI, ANI, tfrags]

                else:
                    tracker[gS] = [gA, gB, rANI, ANI, tfrags]

    with open(outfile, 'w') as out:

        metaA = 'Phylogroup A\tPhylogroup B\tNiche A\tNiche B'
        metaB = 'Site A\tSite B\tTimePoint A\tTimePoint B'
        out.write(f'Genome A\tGenome B\trANI\tANI\t{metaA}\t{metaB}\n')

        for gS, vals in tracker.items():
            
            gA, gB, rANI, ANI = vals[0], vals[1], vals[2], vals[3]
            
            mA = data[gA]
            sA, nA, tA, pA = mA[0], mA[1], mA[2], mA[3]
            mB = data[gB]
            sB, nB, tB, pB = mB[0], mB[1], mB[2], mB[3]

            o = (
                f'{gA}\t{gB}\t{rANI}\t{ANI}\t{pA}\t{pB}\t'
                f'{nA}\t{nB}\t{sA}\t{sB}\t{tA}\t{tB}\n'
                )

            out.write(o)

    return True


def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-aa', '--ANI_file',
        help='Please specify the ANI file!',
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
        '-o', '--output_file',
        help='Please specify the output file name!',
        metavar='',
        type=str,
        required=True
        )
    args=vars(parser.parse_args())

    # define input parameters
    anifile = args['ANI_file']
    metadata = args['metadata_file']
    outfile = args['output_file']

    print(f'\n\nRunning script ...')
    
    data = parse_metadata(metadata)

    _ = parse_ani(anifile, data, outfile)

    print(f'\n\nComplete success space cadet!! Finished without errors.\n\n')
    
if __name__ == "__main__":
    main()


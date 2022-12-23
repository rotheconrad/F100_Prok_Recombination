 #!/usr/bin/env python

'''Parse ani.rb output to tsv file

Output columns:
genome 1, genome 2, one-way ANI 1, one-way ANI 2, Two way ANI,
frags 1, frags 2, frags 3, total frags 1, total frags 2

Where frags 1-3 correspond to the fragments used for the 3 ANIs and
total frags correspond to genome 1 and genome 2.


-------------------------------------------
Author :: Roth Conrad
Email :: rotheconrad@gatech.edu
GitHub :: https://github.com/rotheconrad
Date Created :: Dec 2022
License :: GNU GPLv3
Copyright 2022 Roth Conrad
All rights reserved
-------------------------------------------
'''

import argparse

def parse_anirb(infile, outfile):

    # parses the txt output from ani.rb to tsv

    with open(infile, 'r') as f, open(outfile, 'w') as o:
        # first time through lists will write the header
        # they get replaced with data for each genome pair thereafter
        gnames = ["Genome 1", "Genome 2"]
        anis = ["One-Way 1", "One-Way 2", "Two-Way"]
        frags = ["Frags 1", "Frags 2", "Frags 3"]
        tfrags = ["Total Frags 1", "Total Frags 2"]
        
        for line in f:
            # "Temporal" starts each new comparison. Write out the previous
            # genome pair data each time we start a new one and reset lineout
            if line.startswith('Temporal'):
                # join lists to write out the data
                lineout = gnames + anis + frags + tfrags
                o.write('\t'.join(lineout) + '\n')
                # reset lists for next genome pair
                gnames, anis, frags, tfrags = [], [], [], []
            elif line.startswith('  Reading'):
                # get the genome/long read file name
                gname = line.split(': ')[1].split('.')[0]
                gnames.append(gname)
            elif line.startswith('    Created'):
                tfrag = line.split(' ')[5]
                tfrags.append(tfrag)
            elif line.startswith('!'):
                # get the ANI and the fragments used
                X = line.rstrip().split('%')
                ani = X[0].split(': ')[1]
                frag = X[2].split(', ')[1].split(' ')[1]
                anis.append(ani)
                frags.append(frag)

    return True


def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-i', '--input_anirb_file',
        help='Please specify the input tabular BlastN file!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-o', '--output_file_name',
        help='Please specify the output file name (use .pdf)!',
        metavar='',
        type=str,
        required=True
        )
    args=vars(parser.parse_args())

    # Do what you came here to do:
    print('\n\nRunning Script...\n')

    # define parameters
    infile= args['input_anirb_file']
    outfile = args['output_file_name']

    _ = parse_anirb(infile, outfile)

    print(f'\n\nComplete success space cadet!! Finished without errors.\n\n')

if __name__ == "__main__":
    main()
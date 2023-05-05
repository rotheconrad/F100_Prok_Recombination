#!/usr/local/pacerepov1/python/2.7/bin/python

 #!/usr/bin/env python

'''Computes reciprocal best match between two gene fastas

Requires blast-plus in the system PATH.

Runs two-way blast search, selects shared top matches, filters out short
alignments (alignment length / shorter sequence length ≥ 50%), and 
retains ties for reciprocal best matches RBMs.

inputs are two gene fasta files as nucleotides.
output is a tabular blast file with the RBMs.

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

import argparse, subprocess
from pathlib import Path


def best_hits(data, query, bitscore, line):
    ''' compares bitscores of previous entries and keeps the best 
        tied bitscores are retained. '''

    # if the query is in the dict, check bitscore and replace or append
    if query in data:
        old_bitscore = float(data[query][0].split('\t')[11])

        if bitscore > old_bitscore:
            data[query] = [line]

        elif bitscore == old_bitscore:
            data[query].append(line)

    # otherwise add query to the dict
    else:
        data[query] = [line]

    return data


def parse_tab_blast(blast):
    ''' parses a tabular blast file into dictionary '''

    # this will store the best bitscore matches
    # the dictionary is {query: [line]} storing lines in a list
    # to retain multiple lines if they have tied bitscores
    data = {}

    with open(blast, 'r') as file:
        for line in file:
            X = line.rstrip().split()
            query = X[0]
            subject = X[1]
            alen = int(X[3]) # alignment length
            bitscore = float(X[11])
            qlen = int(X[12]) # query sequence length
            slen = int(X[13]) # subject sequence lenght
            short = min(qlen, slen) # get the short sequence lenght
            pmatch = alen / short # compute alignment ratio

            if pmatch >= 0.5:
                data = best_hits(data, query, bitscore, line)

    return data


def find_RBMs(query, subject1, b2_results, line, b1_RBMs, b2, out, P):

    ''' loops over b2_results to find RBMs '''

    # if no ties proceed with finding best match
    if len(b2_results) == 1:
        # subject2 should match query for an RBM
        subject2 = b2_results[0].split('\t')[1]
        # if it is a match
        if query == subject2:
            # write out the RBM line
            out.write(line)
            # Save queries that find RBMs to remove from dict
            b1_RBMs[query] = ''
            # remove the query from b2 dict
            b2.pop(subject1)
            # print the tie if P True
            if P: print(f'{line}{b2_results[0]}')
    # if there are b2 ties, sort through the ties to find an RBM
    elif len(b2_results) > 1:
        # create a list of subject2's to match with query
        rlist = []
        for match in b2_results:
            subject2 = match.split('\t')[1]
            rlist.append(subject2)
        # if query is in the list its an RBM
        if query in rlist:
            # write out the RBM line
            out.write(line)
            # Save queries that find RBMs to remove from dict
            b1_RBMs[query] = ''
            # get index position of b2 match
            x = rlist.index(query)
            # print tied matches to screen
            print(f'{line}{b2_results[x]}')
            # remove the match from b2 and update b2
            del b2_results[x]
            b2[subject1] = b2_results + ['tie\tremoved']
    else:
        print("There shouldn't be an else number 1")

    return b2, b1_RBMs


def do_the_RBM_thing(b1, b2, out_file):
    ''' finds RBMs, filters ≥ 50% alignment length of short sequence,
    retains ties for RBMs. '''

    # b1 queries with RBMs in b2. to remove from b1
    b1_RBMs = {}

    # open the output file
    out = open(out_file, 'w')

    ####################################################################
    print('\nThe following tied RBMs have been retained in the output:\n')
    # iterate through b1 and find RBMs
    for query, lines in b1.items():
        # iterate through lines
        # lines is a list of tab blast lines with tied bitscores
        # or lines is a list with a single entry if no ties
        # track if lines > 1 so we can print tied matches
        P = True if len(lines) > 1 else False
        for line in lines:
            # get the subject for the current query
            subject1 = line.split('\t')[1]
            # pull the reverse results from b2
            if subject1 in b2:
                b2_results = b2[subject1]
                b2, b1_RBMs = find_RBMs(
                                        query, subject1, b2_results,
                                        line, b1_RBMs, b2, out, P
                                        )
    ####################################################################

    # remove b1 RBMs from b1
    for r in b1_RBMs.keys():
        b1.pop(r)

    print('\n\n################## Finished b1 output ######################')
    print('\n### b1 lines remaining\n')
    # keeping track of how many tied RBMs
    count1, count2 = [], []

    # check out b1
    for query, lines in b1.items():
        # keep track of ties
        if len(lines) > 1: count1.append(len(lines))
        # print remaining blast matches
        for line in lines:
            print(line.rstrip())

    print('\n\n### b2 lines remaining\n')
    # check out b2
    for query, lines in b2.items():
        # filter out entries with the tie already removed
        if lines[0] == 'tie\tremoved':
            continue
        # keep track of ties
        if len(lines) > 1: count2.append(len(lines))
        # print remaining blast matches
        for line in lines:
            print(line.rstrip())

    print(f'\n\nNon-RBM Ties remaining in blast table 1: {count1}')
    print(f'Non-RBM Ties remaining in blast table 2: {count2}')

    # close output file
    out.close()

    return True


def main():

    # Configure Argument Parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-g1', '--gene_fasta_one',
        help='Please specify the first gene fasta file!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-g2', '--gene_fasta_two',
        help='Please specify the second gene fasta file!',
        metavar='',
        type=str,
        required=True
        )
    parser.add_argument(
        '-o', '--output_file',
        help='Please specify the output file!',
        metavar='',
        type=str,
        required=True
        )
    args=vars(parser.parse_args())

    # Do what you came here to do:
    print('\n\nRunning Script...\n')

    # define input params
    g1 = args['gene_fasta_one']
    g2 = args['gene_fasta_two']
    out = args['output_file']

    n1 = g1.split('/')[-1].split('.')[0]
    n2 = g2.split('/')[-1].split('.')[0]
    rtmp = f'{n1}-{n2}'
    # create temp directory
    #Path("rtmp").mkdir(parents=True, exist_ok=True)
    _ = subprocess.run(f'mkdir {rtmp}', shell=True)

    # first blast
    print(f'\nBuidling blast database for {g1} ...')
    db1_str = f'makeblastdb -in {g1} -dbtype nucl -out {rtmp}/db1'
    db1 = subprocess.run(db1_str, shell=True)
    print(f'\nRunning blast with {g2} as query against {g1} ...')
    b1_str = (
        f"blastn -max_target_seqs 10 -db {rtmp}/db1 -query {g2} -out {rtmp}/b1.blast "
        f"-subject_besthit -outfmt '6 qseqid sseqid pident length mismatch "
        f"gapopen qstart qend sstart send evalue bitscore qlen slen'"
        )
    blast1 = subprocess.run(b1_str, shell=True)
    print(f'\nParsing 1st blast result to dictionary ...')
    b1 = parse_tab_blast(f'{rtmp}/b1.blast')

    # second blast
    print(f'\nBuidling blast database for {g2} ...')
    db2_str = f'makeblastdb -in {g2} -dbtype nucl -out {rtmp}/db2'
    db2 = subprocess.run(db2_str, shell=True)
    print(f'\nRunning blast with {g1} as query against {g2} ...')
    b2_str = (
        f"blastn -max_target_seqs 10 -db {rtmp}/db2 -query {g1} -out {rtmp}/b2.blast "
        f"-subject_besthit -outfmt '6 qseqid sseqid pident length mismatch "
        f"gapopen qstart qend sstart send evalue bitscore qlen slen'"
        )
    blast2 = subprocess.run(b2_str, shell=True)
    print(f'\nParsing 2nd blast result to dictionary ...')
    b2 = parse_tab_blast(f'{rtmp}/b2.blast')

    # read in the fasta sequences
    print(f'\nComputing reciprocal best matches ...')
    _ = do_the_RBM_thing(b1, b2, out)

    # remove the tmp directory
    _ = subprocess.run(f'rm -r {rtmp}/', shell=True)

    
    print(f'\n\nComplete success space cadet!! Finished without errors.\n\n')

    
if __name__ == "__main__":
    main()


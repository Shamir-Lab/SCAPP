#!/usr/bin/env python

###
# Provide a command line script to classify sequences in a fasta file
###

from scapp_utils import readfq
from plasclass import plasclass

import argparse

def parse_user_input():

    parser = argparse.ArgumentParser(
        description=
        'classify_fasta classifies the sequences in a fasta file as plasmid origin or not',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('-f','--fasta',
     help='fasta file of the sequences to be classified',
     required=True, type=str
    )
    parser.add_argument('-o','--outfile',
     help='output file prefix',
     required=False, type=str
    )
    parser.add_argument('-p','--num_processes',
     help='Number of processes to use',
     required=False, type=int, default=8
    )

    return parser.parse_args()

def classify(infile, outfile, num_procs):
    ''' Run the classification
    '''
    c = plasclass.plasclass(num_procs)
    seq_names = []
    seqs = []
    i = 0
    fp = open(infile)
    with open(outfile,'w') as o:
        for name, seq, _ in readfq(fp):
            seq_names.append(name)
            seqs.append(seq)
            i += 1
            if i % 50000 == 0:
                probs = c.classify(seqs)
                for j,p in enumerate(probs):
                    o.write(seq_names[j] + '\t' + str(p) + '\n')
                seq_names = []
                seqs = []


        # last bunch of sequences:
        probs = c.classify(seqs)
        for j,p in enumerate(probs):
            o.write(seq_names[j] + '\t' + str(p) + '\n')

    fp.close()


def main(args):
    ''' Main function that classifies the sequences
    '''
    infile = args.fasta
    if args.outfile: outfile = args.outfile
    else: outfile = infile + '.probs.out'
    n_procs = args.num_processes

    classify(infile, outfile, n_procs)

    print("Finished classifying")
    print("Class scores written in: {}".format(outfile))

if __name__=='__main__':
    args = parse_user_input()
    main(args)

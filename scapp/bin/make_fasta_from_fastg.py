#!/usr/bin/env python
import re, argparse, os
from recyclelib.utils import readfq

def parse_user_input():
    parser = argparse.ArgumentParser(
        description=
        'make_fasta_from_fastg converts fastg assembly graph to fasta format'
        )
    parser.add_argument('-g','--graph',
     help='(spades 3.50+) FASTG file to process [recommended: before_rr.fastg]',
     required=True, type=str
     )
    parser.add_argument('-o','--output',
     help='output file name for FASTA of cycles',
     required=False, type=str
     )

    return parser.parse_args()

def parse_lines(fastg, ofile_name):
    ofile = open(ofile_name, 'w')
    fp = open(fastg, 'r')
    seen = set() ##
    for name,seq,qual in readfq(fp):
        name = re.sub('[:,]'," ", name[:-1]).split(" ")[0]
        if name[-1] == "'": name = name[:-1]
        if name in seen: continue
        else: seen.add(name)
        line = ">"+name+"\n"+seq+"\n"
        ofile.write(line)
    fp.close()
    ofile.close()

if __name__=='__main__':
    args = parse_user_input()
    fastg = args.graph

    # output 1 - fasta of sequences
    if args.output:
        fasta_ofile = args.output
    else:
        (root,ext) = os.path.splitext(fastg.name)
        fasta_ofile = root + ext.replace(".fastg", ".nodes.fasta")

    parse_lines(fastg, fasta_ofile)

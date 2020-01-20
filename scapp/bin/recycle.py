#!/usr/bin/env python

import argparse, os
import copy
import logging
import multiprocessing as mp
import pysam


from scapp_utils import *

import PARAMS

def parse_user_input():
    parser = argparse.ArgumentParser(
        description=
        'SCAPP extracts likely plasmids (and other circular DNA elements) from de novo assembly graphs'
        )
    parser.add_argument('-g','--graph',
     help='(spades 3.50+) assembly graph FASTG file to process; recommended for spades 3.5: before_rr.fastg, for spades 3.6+:assembly_graph.fastg',
     required=True, type=str
     )
    parser.add_argument('-k','--max_k',
        help='integer reflecting maximum k value used by the assembler',
        required=True, type=int, default=55
        )
    parser.add_argument('-b','--bam',
        help='BAM file resulting from aligning reads to contigs file, filtering for best matches',
        required=True, type=str
        )
    parser.add_argument('-l', '--length',
     help='minimum length required for reporting [default: 1000]',
     required=False, type=int, default=1000
     )
    parser.add_argument('-m', '--max_CV',
     help='coefficient of variation used for pre-selection [default: 0.5, higher--> less restrictive]',
      required=False, default=1./2, type=float
      )
    parser.add_argument('-i','--iso',
        help='True or False value reflecting whether data sequenced was an isolated strain',
        required=False, type=bool, default=False
        )
    parser.add_argument('-o','--output_dir',
        help='Output directory',
        required=True, type=str
        )
    parser.add_argument('-p','--num_processes',
        help='Number of processes to use',
        required=False, type=int, default=1
        )
    parser.add_argument('-s','--scores',
            help='Contig plasmid scores file',
            required=False, type=str
        )
    parser.add_argument('-gh','--gene_hits',
            help='Contig plasmid gene hits file',
            required=False, type=str
        )

    # internal thresholds to be passed down to utils functions
    parser.add_argument('-clft','--classification_thresh',
        help='threshold for classifying potential plasmid [0.5]',
        required=False, type=float
    )
    parser.add_argument('-gm','--gene_match_thresh',
        help='threshold for % identity and fraction of length to match plasmid genes [0.75]',
        required=False, type=float
    )
    parser.add_argument('-sls','--selfloop_score_thresh',
        help='threshold for self-loop plasmid score [0.9]',
        required=False, type=float
    )
    parser.add_argument('-slm','--selfloop_mate_thresh',
        help='threshold for self-loop off loop mates [0.1]',
        required=False, type=float
    )
    parser.add_argument('-cst','--chromosome_score_thresh',
        help='threshold for high confidence chromosome node score [0.2]',
        required=False, type=float
    )
    parser.add_argument('-clt','--chromosome_len_thresh',
        help='threshold for high confidence chromosome node length [10000]',
        required=False, type=float
    )
    parser.add_argument('-pst','--plasmid_score_thresh',
        help='threshold for high confidence plasmid node score [0.9]',
        required=False, type=float
    )
    parser.add_argument('-plt','--plasmid_len_thresh',
        help='threshold for high confidence plasmid node length [10000]',
        required=False, type=float
    )
    parser.add_argument('-cd','--good_cyc_dominated_thresh',
        help='threshold for # of mate-pairs off the cycle in dominated node [0.5]',
        required=False, type=float
    )

    return parser.parse_args()

def run_scapp(fastg, outdir, bampath, num_procs, max_k, \
                    genes_file, use_genes, scores_file, use_scores, \
                    max_CV, min_length, ISO=False):
    ''' Run SCAPP'''

    logger = logging.getLogger("scapp_logger")

    basename, _ = os.path.splitext(os.path.basename(fastg))
    fasta_ofile = os.path.join(outdir, basename+".cycs.fasta")
    cycs_ofile = os.path.join(outdir, basename+".cycs.paths.txt")
    loop_ofile = os.path.join(outdir,basename+".self_loops.fasta")
    f_cycs_fasta = open(fasta_ofile, 'w') # output 1 - fasta of sequences
    f_cyc_paths = open(cycs_ofile, 'w') # output 2 - file containing path name (corr. to fasta),
    f_long_self_loops = open(loop_ofile,'w') # output 3 - file of self-loop fasta sequences
    bamfile = pysam.AlignmentFile(bampath)

    # graph processing begins
    G = get_fastg_digraph(fastg)

    logger.info("Removing %d isolate nodes" % len(list(nx.isolates(G))))
    G.remove_nodes_from(list(nx.isolates(G)))

    cov_vals = [get_cov_from_spades_name(n) for n in G.nodes()]
    MED_COV = np.median(cov_vals)
    STD_COV = np.std(cov_vals)
    # set thresholds for max. CV, min
    # path coverage for allowing cross mappings
    if ISO:
        thresh = np.percentile(cov_vals, 95)
    else:
        thresh = np.percentile(cov_vals, 75)

    logger.info("Coverage threshold:%4f" % thresh)
    print(MED_COV, STD_COV, thresh)
    path_count = 0
    SEQS = get_fastg_seqs_dict(fastg,G)


    # add a score to every node, remove long nodes that are most probably chrom.
    if use_scores:
        get_node_scores(scores_file,G)

    # keep track of the nodes that have plasmid genes on them
    if use_genes:
        get_gene_nodes(genes_file,G)

    # gets set of long simple loops, removes short
    # simple loops from graph
    long_self_loops = get_long_self_loops(G, min_length, SEQS, bamfile, use_scores, use_genes, max_k)

    for nd in long_self_loops:
        name = get_spades_type_name(path_count, nd,
        SEQS, max_k, G, get_cov_from_spades_name(nd[0]))
        path_count += 1

        seq = get_seq_from_path(nd, SEQS, max_k_val=max_k)
        print(nd)
        print(" ")
        if len(seq)>=min_length:
            f_cycs_fasta.write(">" + name + "\n" + seq + "\n")
            f_long_self_loops.write(">" + name + "\n" + seq + "\n")
            f_cyc_paths.write(name + "\n" +str(nd[0])+ "\n" +
             str(get_num_from_spades_name(nd[0])) + "\n")

    comps = (G.subgraph(c).copy() for c in nx.strongly_connected_components(G))

    if use_scores:
        # Remove nodes that are most likely chromosomal
        # This may
        smaller_comps = set()
        for comp in comps:
            remove_hi_confidence_chromosome(comp)
            smaller_comps.update((comp.subgraph(c).copy() for c in nx.strongly_connected_components(comp)))
        comps = smaller_comps

    ###################################
    # iterate through SCCs looking for cycles
    ###################################

    #multiprocessing to find shortest paths
    pool = mp.Pool(num_procs)

    print("================== Added paths ====================")

    VISITED_NODES = set([]) # used to avoid problems due to RC components
    redundant = False

    for c in sorted(comps,key=len):

	    # check if any nodes in comp in visited nodes
        # if so continue
        for node in c.nodes():
             if node in VISITED_NODES:
                 redundant = True
                 break
        if redundant:
             redundant = False
             continue # have seen the RC version of component
        COMP = nx.DiGraph()
        COMP = c.to_directed()

        rc_nodes = [rc_node(n) for n in COMP.nodes()]
        VISITED_NODES.update(COMP.nodes())
        VISITED_NODES.update(rc_nodes)

        path_set = process_component(COMP, G, max_k, min_length, max_CV, SEQS, thresh, bamfile, pool, use_scores, use_genes, num_procs)

        for p in path_set:
            name = get_spades_type_name(path_count, p[0], SEQS, max_k, G, p[1])
            seq = get_seq_from_path(p[0], SEQS, max_k_val=max_k)
            print(p[0])
            print(" ")
            if len(seq)>=min_length:
                f_cycs_fasta.write(">" + name + "\n" + seq + "\n")
                f_cyc_paths.write(name + "\n" +str(p[0])+ "\n" +
                 str([get_num_from_spades_name(n) for n in p[0]]) + "\n")
            path_count += 1

    pool.close()
    pool.join()
    f_cycs_fasta.close()
    f_cyc_paths.close()
    f_long_self_loops.close()


def main():

    ####### entry point  ##############

    args = parse_user_input()
    num_procs = args.num_processes
    fastg = args.graph
    max_k = args.max_k
    files_dir = os.path.dirname(fp.name)
    if args.scores:
        use_scores = True
    else: use_scores = False
    if args.gene_hits:
        use_genes = True
    else: use_genes = False

    bampath = args.bam
    ISO = args.iso

    # default threshold vaariables
    PARAMS.load_params_json()
    if args.max_CV:
        PARAMS.MAX_CV = args.max_CV
    if args.length:
        PARAMS.MIN_LENGTH = args.length
    if args.classification_thresh:
        PARAMS.CLASSIFICATION_THRESH = args.classification_thresh
    if args.gene_match_thresh:
        PARAMS.GENE_MATCH_THRESH = args.gene_match_thresh
    if args.selfloop_score_thresh:
        PARAMS.SELF_LOOP_SCORE_THRESH = args.selfloop_score_thresh
    if args.selfloop_mate_thresh:
        PARAMS.SELF_LOOP_MATE_THRESH = args.selfloop_mate_thresh
    if args.chromosome_score_thresh:
        PARAMS.CHROMOSOME_SCORE_THRESH = args.chromosome_score_thresh
    if args.chromosome_len_thresh:
        PARAMS.CHROMOSOME_LEN_THRESH = args.CHROMOSOME_LEN_THRESH
    if args.plasmid_score_thresh:
        PARAMS.PLASMID_SCORE_THRESH = args.plasmid_score_thresh
    if args.plasmid_len_thresh:
        PARAMS.PLASMID_LEN_THRESH
    if args.good_cyc_dominated_thresh:
        PARAMS.GOOD_CYC_DOMINATED_THRESH = args.good_cyc_dominated_thresh

    # Set up logging and write config and options to the log file
    logfile = os.path.join(args.output_dir,"scapp.log")
    logging.basicConfig(filemode='w', filename=logfile, level=logging.INFO, format='%(asctime)s: %(message)s', datefmt='%d/%m/%Y %H:%M')
    logger = logging.getLogger("scapp_logger")

    run_scapp(fastg, args.output_dir, bampath, num_procs, max_k, \
                    args.gene_hits, use_genes, args.scores, use_scores, \
                    max_CV, min_length, ISO)


if __name__ == '__main__':
    main()

from __future__ import division
import numpy as np
import math
import networkx as nx
import re, pysam
import logging

import multiprocessing as mp
from multiprocessing import Manager

import time

import PARAMS

complements = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
logger = logging.getLogger("scapp_logger")

def readfq(fp): # this is a generator function
    """ # lh3's fast fastX reader:
        https://github.com/lh3/readfq/blob/master/readfq.py
    """
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, ''.join(seqs), None # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, seq, None # yield a fasta record instead
                break


def get_node_scores(scores_file,G):
    """ Write the plasmid scores into each node in the graph
    """
    scores = {}
    with open(scores_file) as f:
        for line in f:
            split = line.strip().split()
            scores[split[0]] = float(split[1])
    for nd in G.nodes():
        G.add_node(nd, score=scores[nd])

def get_gene_nodes(genes_file,G):
    """ Annotate each node in the graph whether it has plasmid gene in it
    """
    gene_nodes = set()
    with open(genes_file) as f:
        for line in f:
            gene_nodes.add(line.strip())
    for nd in G.nodes():
        if nd in gene_nodes:
            G.add_node(nd, gene=True)
        else:
            G.add_node(nd,gene=False)

def rc_seq(dna):
    rev = reversed(dna)
    return "".join([complements[i] for i in rev])

def get_num_from_spades_name(name):
    name_parts = name.split("_")
    contig_length = name_parts[1]
    return int(contig_length)

def get_length_from_spades_name(name):
    name_parts = name.split("_")
    contig_length = name_parts[3]
    return int(contig_length)

def get_cov_from_spades_name(name):
    name_parts = name.split("_")
    cov = name_parts[5]
    if cov[-1]=="'": cov=cov[:-1]
    return float(cov)

def get_fastg_digraph(fastg_name):
    """ scans through fastg headers as an adjacency list
        builds and returns a nx directed graph using adjacencies
        note: no connections are created between each node and its
        rc node - we need to take care to maintain these
    """
    lines = []
    fp = open(fastg_name, 'r')
    for name,seq,qual in readfq(fp):
        name = re.sub('[:,]'," ", name[:-1])
        lines.append(name)
    G = nx.DiGraph()
    return nx.parse_adjlist(lines, create_using=G)

def get_fastg_seqs_dict(fastg_name, G):
    """ returns a dictionary of sequences in graph
        where node names are keys and sequence strings
        are values; useful for saving memory when G
        is a subgraph (e.g., a component)
    """
    fp = open(fastg_name, 'r')
    seqs = {}
    for name,seq,qual in readfq(fp):
        name_parts = re.sub('[:,]'," ", name[:-1]).split()
        node = name_parts[0]
        seqs[node] = seq
    return seqs

def rc_node(node):
    """ gets reverse complement
        spades node label
    """
    if node[-1] == "'": return node[:-1]
    else: return node + "'"


def get_cov_from_spades_name_and_graph(name,G):
    if name not in G:
        return 0
    if 'cov' in G.nodes[name]:
        return G.nodes[name]['cov']
    else:
        return get_cov_from_spades_name(name)

def update_node_coverage(G, node, new_cov):
    """ changes coverage value stored in 'cov'
        field on both F and R version of a node
        if new_cov is 0, both node versions are removed
    """
    if node not in G.nodes(): # nothing to be done, perhaps already removed
        return
    if new_cov == 0:
        G.remove_node(node)
        if rc_node(node) in G.nodes():
            G.remove_node(rc_node(node))
    else:
        G.add_node(node, cov=new_cov)
        if rc_node(node) in G.nodes():
          G.add_node(rc_node(node), cov=new_cov)

def get_spades_base_mass(G, name):
    length = get_length_from_spades_name(name)
    coverage = get_cov_from_spades_name_and_graph(name,G)
    if coverage <= 0.0: coverage = 1.0/float(length) # ensure no division by zero, consider more principled way to do this
    return length * coverage

def get_seq_from_path(path, seqs, max_k_val=55, cycle=True):
    """ retrieves sequence from a path;
        instead of specifying cycles by having the first and
        last node be equal, the user must decide if a cycle is
        the intent to avoid redundant k-mers at the ends
    """
    start = seqs[path[0]]
    if len(path)==1:
        if cycle:
            return start[max_k_val:]
        else:
            return start
    else:
        seq = ''
        for p in path:
            seq += seqs[p][max_k_val:]
        if cycle: return seq
        else: return start[:max_k_val] + seq

def get_wgtd_path_coverage_CV(path, G, seqs, max_k_val=55):
    if len(path)< 2: return 0
    mean, std = get_path_mean_std(path, G, seqs, max_k_val,discount=True)
    if mean<=0: return 0
    return std/mean

def get_node_cnts_hist(path):
    d = {}
    for p in path:
        # always count based on positive node
        pos_name = p if (p[-1]!="'") else p[:-1]
        d[pos_name] = d.get(pos_name,0) + 1
    return d


def get_discounted_node_cov(node,path,G):
    """ Return the coverage of the node, discounted by the coverage of neighbouring
        nodes not in the path
    """
    pred_covs = [(get_cov_from_spades_name_and_graph(p,G),p) for p in G.predecessors(node)]
    succ_covs = [(get_cov_from_spades_name_and_graph(s,G),s) for s in G.successors(node)]

    non_path_cov = sum([p[0] for p in pred_covs if p[1] not in path]) + sum([s[0] for s in succ_covs if s[1] not in path])
    in_path_cov = sum([p[0] for p in pred_covs if p[1] in path]) + sum([s[0] for s in succ_covs if s[1] in path])
    node_cov = get_cov_from_spades_name_and_graph(node,G)
    node_cov *= in_path_cov/(non_path_cov + in_path_cov)
    ###################### A possible alternative would be to discount in this way for both the in- and out-neighbours
    ###################### and then average the two discounted
    return node_cov


def get_path_covs(path,G,discount=False):
    if discount:
        covs = [get_discounted_node_cov(n,path,G) for n in path]
        # discount weight of nodes that path passes through multiple times
        cnts = get_node_cnts_hist(path)
        for i in range(len(path)):
            p = path[i]
            pos_name = p if (p[-1]!="'") else p[:-1]
            if cnts[pos_name] > 1:
                covs[i] /= cnts[pos_name]
    else:
        covs = [get_cov_from_spades_name_and_graph(n,G) for n in path]

    return covs

def get_path_mean_std(path, G, seqs, max_k_val=55,discount=True):
    covs = np.array(get_path_covs(path,G,discount))
    wgts = np.array([(get_length_from_spades_name(n)-max_k_val) for n in path])
    tot_len = len(get_seq_from_path(path, seqs, max_k_val, cycle=True))
    if tot_len<=0: return (0,0)
    wgts = np.multiply(wgts, 1./tot_len)
    mean = np.average(covs, weights = wgts)
    std = np.sqrt(np.dot(wgts,(covs-mean)**2))
    return (mean,std)

def update_path_coverage_vals(path, G, seqs, max_k_val=55,):
    mean, _ = get_path_mean_std(path, G, seqs, max_k_val) ## NOTE: CAN WE STILL GUARANTEE CONVERGENCE WHEN DISCOUNTING COVERAGE ??!
    covs = get_path_covs(path,G)
    new_covs = covs - mean
    logger.info("Path: %s Mean: %s Covs: %s" % (str(path),str(mean),str(covs) ))
    for i in range(len(path)):
        if new_covs[i] > 0:
            update_node_coverage(G,path[i],new_covs[i])
        else:
            update_node_coverage(G,path[i],0)
    return new_covs

def update_path_with_covs(path, G, covs):
    for i in range(len(path)):
        if covs[i] > 0:
            update_node_coverage(G,path[i],covs[i])
        else:
            update_node_coverage(G,path[i],0)

def get_total_path_mass(path,G):
    return sum([get_length_from_spades_name(p) * \
        get_cov_from_spades_name_and_graph(p,G) for p in path])


def get_long_self_loops(G, min_length, seqs, bamfile, use_scores=True, use_genes=True, max_k_val=55, score_thresh=0.9, mate_thresh = 0.1):
    """ returns set of self loop nodes paths that are longer
        than min length and satisfy mate pair requirements;
        removes those and short self loops from G
    """
    potential_plasmids = set([])
    to_remove = []

    for nd in list(nx.nodes_with_selfloops(G)):
        if (rc_node(nd),) in potential_plasmids: continue
        nd_path = (nd,)
        path_len = len(get_seq_from_path(nd_path, seqs, max_k_val))

        # check whether it is isolated or connected to other nodes:
        isolated_loop = False
        if G.in_degree(nd) == 1 and G.out_degree(nd)== 1:
            isolated_loop = True
        if isolated_loop:
            if path_len < min_length:
                to_remove.append(nd)
                continue

            # take nodes that have plasmid genes or very high plasmid scores
            if use_scores and use_genes:
                logger.info("SLS: %f" % PARAMS.SELF_LOOP_SCORE_THRESH)
                if G.nodes[nd]['score'] > PARAMS.SELF_LOOP_SCORE_THRESH or G.nodes[nd]['gene']==True:
                    potential_plasmids.add(nd_path)
                    logger.info("Added path: %s - high scoring long self-loop" % nd)
                    to_remove.append(nd)
                    continue

            off_node_mate_count, on_node_mate_count = count_selfloop_mates(nd,bamfile,G)
            if float(off_node_mate_count) > PARAMS.SELF_LOOP_MATE_THRESH*float(on_node_mate_count):
                logger.info('Self loop %s has %2f percent off-node mate-pairs. Removing' % (nd,PARAMS.SELF_LOOP_MATE_THRESH))
                to_remove.append(nd)
            else:
                potential_plasmids.add(nd_path)
                logger.info("Added path: %s  - long self loop" % nd)
                to_remove.append(nd)
        else: # non-isolated loop
            if path_len < min_length: continue

            off_node_mate_count, on_node_mate_count = count_selfloop_mates(nd,bamfile,G)
            if float(off_node_mate_count) > PARAMS.SELF_LOOP_MATE_THRESH*float(on_node_mate_count):  # TODO: could be different than for isolated loop
                                                                                    # Maybe - func of node length (and read length, insert size???)
                logger.info('Self loop %s has %2f percent off-node mate-pairs.' % (nd,PARAMS.SELF_LOOP_MATE_THRESH))
            else:
                potential_plasmids.add(nd_path)
                logger.info("Added path: %s  - long self loop" % nd)
                to_remove.append(nd)

    for nd in to_remove:
        update_node_coverage(G, nd, 0)
    logger.info("Removing %d self-loop nodes" % len(to_remove))
    return potential_plasmids

def remove_hi_confidence_chromosome(G):
    """ Remove the long nodes that are predicted to likely be chromosomal
    """
    to_remove = []
    for nd in G.nodes():
        if get_length_from_spades_name(nd) > PARAMS.CHROMOSOME_LEN_THRESH and \
            G.nodes[nd]['score'] < PARAMS.CHROMOSOME_SCORE_THRESH:
            to_remove.append(nd)
            to_remove.append(rc_node(nd))
    G.remove_nodes_from(to_remove)
    logger.info("Removed %d long, likely chromosomal nodes" % len(set(to_remove)))

def get_hi_conf_plasmids(G):
    """ Return a list of nodes that are likely plasmids
    """

    hi_conf_plasmids = [nd for nd in G.nodes() if (get_length_from_spades_name(nd) > PARAMS.PLASMID_LEN_THRESH and \
                        G.nodes[nd]['score'] > PARAMS.PLASMID_SCORE_THRESH)]

    logger.info("Found %d long, likely plasmid nodes" % len(hi_conf_plasmids))
    return hi_conf_plasmids

def get_plasmid_gene_nodes(G):
    """ Return list of nodes annotated as having a plasmid gene
    """
    plasmid_gene_nodes = [nd for nd in G.nodes() if G.nodes[nd]['gene']==True]
    logger.info("Found %d nodes with plasmid genes" % len(plasmid_gene_nodes))
    return plasmid_gene_nodes

def get_unoriented_sorted_str(path):
    """ creates unique, orientation-oblivious string representation of path,
        used to make sure node covered whenever rc of node is;
        lets us avoid issue of rc of node having different weight than node
    """
    all_rc_path = []
    for p in path:
        if p[-1] != "'": p = p+"'"
        all_rc_path.append(p)
    return "".join(sorted(all_rc_path))

def get_shortest(args_array):
    """ Worker function for getting shortest path to each node in parallel
    """

    node, G, paths_list = args_array
    shortest_score = float("inf")
    path = None
    for pred in G.predecessors(node):
         try:
             path_len,shortest_path = nx.bidirectional_dijkstra(G, node, pred, weight='cost')
             if path_len < shortest_score:
                 path = shortest_path
                 shortest_score = path_len
         except nx.exception.NetworkXNoPath:
             continue
    if path is not None: paths_list.append(path)
    # done

def enum_high_mass_shortest_paths(G, pool, use_scores=False, use_genes=False, seen_paths=None):
    """ given component subgraph, returns list of paths that
        - is non-redundant (includes) no repeats of same cycle
        - includes shortest paths starting at each node n (assigning
        node weights to be 1/(length * coverage)) to each of
        its predecessors, and returning to n
    """
    if seen_paths == None:
        seen_paths = []
    unq_sorted_paths = set([])
    # in case orientation obliv. sorted path strings passed in
    for p in seen_paths:
        unq_sorted_paths.add(p)
    paths = []

    nodes = []
    nodes = list(G.nodes()) # creates a copy

    logger.info("Getting edge weights")

    # use add_edge to assign edge weights to be 1/mass of starting node
    # TODO: only calculate these if they haven't been/need to be updated
    for e in G.edges():
        if use_genes and G.nodes[e[1]]['gene'] == True:
            G.add_edge(e[0], e[1], cost = 0.0)
        elif use_scores==True:
            G.add_edge(e[0], e[1], cost = (1.-(G.nodes[e[1]]['score']))/get_spades_base_mass(G, e[1]))
        else:
            G.add_edge(e[0], e[1], cost = (1./get_spades_base_mass(G, e[1])))

    logger.info("Getting shortest paths")
    paths_list = []
    if pool._processes > 1 and pool._processes <= 2*len(nodes): # otherwise, run single threaded
        paths_list=Manager().list()
        pool.map(get_shortest, [[node, G, paths_list] for node in nodes])
    else:
        for node in nodes:
            get_shortest([node,G,paths_list])

    for path in paths_list:
        # below: create copy of path with each node as rc version
        # use as unique representation of a path and rc of its whole
        unoriented_sorted_path_str = get_unoriented_sorted_str(path)

        # here we avoid considering cyclic rotations of identical paths
        # by sorting their string representations
        # and comparing against the set already stored
        if unoriented_sorted_path_str not in unq_sorted_paths:
            unq_sorted_paths.add(unoriented_sorted_path_str)
            paths.append(tuple(path))

    return paths


def get_high_mass_shortest_path(node,G,use_scores,use_genes):
    """ Return the shortest circular path back to node
    """
    # TODO: potentially add check for unique paths so that don't check same cycle
    # twice if there are two potential plasmid nodes in it

    for e in G.edges():
        if use_genes and G.nodes[e[1]]['gene'] == True:
            G.add_edge(e[0], e[1], cost = 0.0)
        elif use_scores == True:
            G.add_edge(e[0], e[1], cost = (1.-(G.nodes[e[1]]['score']))/get_spades_base_mass(G, e[1]))
        else:
            G.add_edge(e[0], e[1], cost = (1./get_spades_base_mass(G, e[1])))

    shortest_score = float("inf")
    path = None
    for pred in G.predecessors(node):
        try:
            path_len = nx.shortest_path_length(G, source=node, target=pred, weight='cost')
            if path_len < shortest_score:
                path = tuple(nx.shortest_path(G, source=node,target=pred, weight='cost'))
                shortest_score = path_len
        except nx.exception.NetworkXNoPath:
            continue

    return path

def get_non_repeat_nodes(G, path):
    """ returns a list of all non-repeat (in degree and out-degree
        == 1) nodes in a path; if there are no such nodes,
        returns an empty list
        NB: G input should be whole graph, not specific SCC, to avoid
        disregarding isolated nodes
    """
    sing_nodes = []
    for nd in path:
        if G.out_degree(nd)==1 and G.in_degree(nd)==1:
            sing_nodes.append(nd)
    return sing_nodes


def get_spades_type_name(count, path, seqs, max_k_val, G, cov=None):
    path_len = len(get_seq_from_path(path,seqs,max_k_val))
    if cov==None:
        cov = get_total_path_mass(path,G)/float(path_len)
    info = ["RNODE", str(count+1), "length", str(path_len),
     "cov", '%.5f' % (cov)]
    return "_".join(info)


def count_selfloop_mates(node,bamfile,G):
    """ Counts the number of off-node and on-node mate pairs of
        a self-loop node
    """
    off_node_count = 0
    on_node_count = 0
    if node[-1] == "'": node = node[:-1]
    try:
        for hit in bamfile.fetch(node):
            nref = bamfile.getrname(hit.next_reference_id)
            if nref != node:
                off_node_count += 1
            else: on_node_count += 1

    except ValueError:
        pass

    return off_node_count, on_node_count

def get_contigs_of_mates(node, bamfile, G):
    """ retrieves set of nodes mapped to by read pairs
        having one mate on node; discards isolated nodes
        because they tend to reflect irrelevant alignments
    """
    mate_tigs = set([])
    if node[-1] == "'": node=node[:-1]
    try:
        for hit in bamfile.fetch(node):
            nref = bamfile.getrname(hit.next_reference_id)
            if nref != node:
                mate_tigs.add(nref)

    except ValueError:
        pass

    to_remove = set([])
    for nd in mate_tigs:
        if (G.in_degree(nd)==0 and G.out_degree(nd)==0) or \
        (not G.has_node(nd)):
            to_remove.add(nd)
        # see if nd reachable by node or vice-versa
        # try both flipping to rc and switching source and target
        elif G.has_node(rc_node(node)) and not any([nx.has_path(G, node, nd), nx.has_path(G, rc_node(node),nd),
                nx.has_path(G, nd, node), nx.has_path(G, nd, rc_node(node))]):
            to_remove.add(nd)
        elif not any([nx.has_path(G,node,nd),nx.has_path(G,nd,node)]):
            to_remove.add(nd)
    mate_tigs -= to_remove

    return mate_tigs

####### Updated - exclude path with node that mostly has mates outside of path
def is_good_cyc(path, G, bamfile):
    """ check all non-repeat nodes only have mates
        mapping to contigs in the cycle
    """

    sing_nodes = set()
    for node in path:
      if node[-1] == "'": node = node[:-1]
      sing_nodes.add(node)
    if len(sing_nodes)==0: return True

    non_path_dominated_nodes = 0

    for nd in sing_nodes:
        mate_tigs = get_contigs_of_mates(nd, bamfile, G)
        # NOTE: ^ this only gets mates that are reachable from nd in G
        logger.info("\tNode: %s" % nd)
        logger.info("\t\tMates: %s" % ", ".join(mate_tigs))

        # need to check against F and R versions of path nodes
        path_rc = [rc_node(x) for x in path]
        num_mates_in_path = sum([1 for x in mate_tigs if (x in path or x in path_rc)])
        num_mates_not_in_path = len(mate_tigs)-num_mates_in_path
        if num_mates_in_path < num_mates_not_in_path:
       ########### if len(mate_tigs)>1 and num_mates_in_path < num_mates_not_in_path:
            non_path_dominated_nodes += 1
    if float(non_path_dominated_nodes)/float(len(sing_nodes)) > PARAMS.GOOD_CYC_DOMINATED_THRESH:
        logger.info("Too many nodes with majority of mates not on path")
        return False
    else: return True

#########################
def process_component(COMP, G, max_k, min_length, max_CV, SEQS, thresh, bamfile, pool, use_scores=False, use_genes=False, num_procs=1):
    """ run recycler for a single component of the graph
        use multiprocessing to process components in parallel
    """

        ###############MOVED FROM OUTER CODE ON WHOLE G
    if use_scores: remove_hi_confidence_chromosome(COMP) ##################################

    # initialize shortest path set considered
    path_count = 0
    seen_unoriented_paths = set([])
    paths_set = set([]) #the set of paths found


    # first look for paths starting from the nodes annotated with plasmid genes
    if use_genes:
        plasmid_gene_nodes = get_plasmid_gene_nodes(COMP)
        potential_plasmid_mass_tuples = [(get_spades_base_mass(COMP,nd),nd) for nd in plasmid_gene_nodes]
        potential_plasmid_mass_tuples.sort(key = lambda n: n[0])
        while potential_plasmid_mass_tuples: # could be removing other nodes from the list
            top_node = potential_plasmid_mass_tuples.pop() # highest mass node
            top_node_name = top_node[1]
            path = get_high_mass_shortest_path(top_node_name,COMP,use_scores,use_genes) #######
            if path is None: continue
            # check coverage variation
            path_CV = get_wgtd_path_coverage_CV(path,G,SEQS,max_k_val=max_k)
            logger.info("Plasmid gene path: %s, CV: %4f" % (str(path),path_CV))
            if path_CV <= max_CV and is_good_cyc(path,G,bamfile):
                logger.info("Added plasmid gene path %s" % (str(path)))

                # prevent checking nodes that have been removed
                i = 0
                while i < len(potential_plasmid_mass_tuples):
                    if potential_plasmid_mass_tuples[i][1] in path or \
                        rc_node(potential_plasmid_mass_tuples[i][1]) in path:
                        potential_plasmid_mass_tuples.pop(i)
                    else: i += 1

                seen_unoriented_paths.add(get_unoriented_sorted_str(path))
                before_cov, _ = get_path_mean_std(path, G, SEQS, max_k)
                covs = update_path_coverage_vals(path, G, SEQS, max_k)
                update_path_with_covs(path, COMP, covs)
                path_count += 1
                paths_set.add((path,before_cov))

            else:
                logger.info("Did not add plasmid gene path: %s" % (str(path)))

        # then look for circular paths that start from hi confidence plasmid nodes
    if use_scores:
        potential_plasmid_nodes = get_hi_conf_plasmids(COMP)
        potential_plasmid_mass_tuples = [(get_spades_base_mass(COMP,nd),nd) for nd in potential_plasmid_nodes]
        potential_plasmid_mass_tuples.sort(key = lambda n: n[0])
        while potential_plasmid_mass_tuples: # could be removing other nodes from the list
            top_node = potential_plasmid_mass_tuples.pop() # highest mass node
            top_node_name = top_node[1]
            path = get_high_mass_shortest_path(top_node_name,COMP,use_scores,use_genes)
            if path is None: continue
            # check coverage variation
            path_CV = get_wgtd_path_coverage_CV(path,G,SEQS,max_k_val=max_k)
            logger.info("Hi conf path: %s, CV: %4f" % (str(path),path_CV))

            if path_CV <= max_CV and is_good_cyc(path,G,bamfile):
                logger.info("Added hi conf path %s" % (str(path)))

                # prevent checking nodes that have been removed
                i = 0
                while i < len(potential_plasmid_mass_tuples):
                    if potential_plasmid_mass_tuples[i][1] in path or \
                        rc_node(potential_plasmid_mass_tuples[i][1]) in path:
                        potential_plasmid_mass_tuples.pop(i)
                    else: i += 1

                seen_unoriented_paths.add(get_unoriented_sorted_str(path))
                before_cov, _ = get_path_mean_std(path, G, SEQS, max_k)
                #before_cov, _ = get_path_mean_std(path, COMP, SEQS, max_k)
                covs = update_path_coverage_vals(path, G, SEQS, max_k)##########################
                update_path_with_covs(path, COMP, covs) ####################################
                path_count += 1
                paths_set.add((path,before_cov))

            else:
                logger.info("Did not add hi-conf path: %s" % (str(path)))

        # 3rd step. Run Recycler algorithm that looks for circular high mass shortest
        # paths and accept them as plasmid predictions if the coverages and mate pairs
        # match the required thresholds
#######################################################################################
#######################################################################################


    paths = enum_high_mass_shortest_paths(COMP, pool, use_scores,use_genes,seen_unoriented_paths)
    last_path_count = 0
    last_node_count = 0

        # continue as long as you either removed a low mass path
        # from the component or added a new path to final paths
    while(path_count!=last_path_count or\
        len(COMP.nodes())!=last_node_count):

        last_node_count = len(COMP.nodes())
        last_path_count = path_count

        # make tuples of (CV, path)
        path_tuples = []
        for p in paths:
            if len(get_seq_from_path(p, SEQS, max_k_val=max_k)) < min_length:
                seen_unoriented_paths.add(get_unoriented_sorted_str(p))
                logger.info("Num seen paths: %d" % (len(seen_unoriented_paths)))
                continue
            path_tuples.append((get_wgtd_path_coverage_CV(p,G,SEQS,max_k_val=max_k), p))

        logger.info("Num path tuples: %d" % (len(path_tuples)))
        if(len(path_tuples)==0): break

        # sort in ascending CV order
        path_tuples.sort(key=lambda path: path[0])

        for pt in path_tuples:
            curr_path = pt[1]
            curr_path_CV = pt[0]
            logger.info("Path: %s" % (",".join(curr_path)))
            if get_unoriented_sorted_str(curr_path) not in seen_unoriented_paths:

                ## only report if low CV and matches mate pair info
                if (curr_path_CV <= (max_CV) and \
                    is_good_cyc(curr_path,G,bamfile)):

                    logger.info("Added path %s" % ", ".join(curr_path))
                    logger.info("\tCV: %4f" % curr_path_CV)
                    seen_unoriented_paths.add(get_unoriented_sorted_str(curr_path))
                    #before_cov, _ = get_path_mean_std(curr_path, COMP, SEQS, max_k)
                    before_cov, _ = get_path_mean_std(curr_path, G, SEQS, max_k)
                    covs = update_path_coverage_vals(curr_path, G, SEQS, max_k)
                    update_path_with_covs(curr_path, COMP, covs)
                    path_count += 1
                    paths_set.add((curr_path,before_cov))
                    break

                else:
                    logger.info("Did not add path: %s" % (", ".join(curr_path)))
                    logger.info("\tCV: %4f" % curr_path_CV)
                    if curr_path_CV > max_CV:
                        break # sorted by CV
                    else: # not good mate pairs
                        seen_unoriented_paths.add(get_unoriented_sorted_str(curr_path))

        # recalculate paths on the component
        print(str(len(COMP.nodes())) + " nodes remain in component")
        logger.info("Remaining nodes: %d" % (len(COMP.nodes())))
        paths = enum_high_mass_shortest_paths(COMP, pool, use_scores,use_genes,seen_unoriented_paths)

    #end while
    return paths_set

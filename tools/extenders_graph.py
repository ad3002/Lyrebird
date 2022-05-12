#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 21.05.2014
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com

from trseeker.tools.raw_reads_tools import raw_reads_continue_kmer_right
from trseeker.tools.raw_reads_tools import raw_reads_continue_kmer_left
from trseeker.tools.sequence_tools import get_revcomp
from collections import defaultdict
from PyBioSnippets.Lyrebird.tools.graphs import Graph


def expand_kmer(qkmer, seen_fragments, all_kmers):
    """
    """
    for i, kmer in enumerate(all_kmers):
        if qkmer in kmer:
            seen_fragments.add(kmer)
    return seen_fragments


def get_right_extension(kmer, kmer2tf, settings):
    """ Extender for repeat graph.
    Returns:
        new_kmer, ptf, tf, total_freq, R
    """
    if kmer2tf:
        kmer_a = kmer[1:]+"A"
        kmer_c = kmer[1:]+"C"
        kmer_t = kmer[1:]+"T"
        kmer_g = kmer[1:]+"G"
        R = [
            (kmer2tf[kmer_a]+kmer2tf[get_revcomp(kmer_a)], kmer_a),
            (kmer2tf[kmer_c]+kmer2tf[get_revcomp(kmer_c)], kmer_c),
            (kmer2tf[kmer_t]+kmer2tf[get_revcomp(kmer_t)], kmer_t),
            (kmer2tf[kmer_g]+kmer2tf[get_revcomp(kmer_g)], kmer_g),
        ]
    else:
        R = raw_reads_continue_kmer_right(kmer, settings["jf_db"])
        R.reverse()
    R.sort()
    total_freq = float(sum([x[0] for x in R]))
    if settings["verbose"]:
        if settings["verbose"] > 3:
            for tf, new_kmer in R:
                print "\t", new_kmer, round(tf/total_freq, 4), tf
    if not total_freq:
        yield None, None, None, None, None
    for tf, new_kmer in R:
        ptf = round(tf/total_freq, 4)
        if ptf < settings["pcutoff"]:
            continue
        if tf < settings["coverage"]:
            continue
        if settings["verbose"] == 3:
            print "\t", new_kmer, round(tf/total_freq, 4), tf
        yield new_kmer, ptf, tf, total_freq, R
    yield None, None, None, None, None


def get_left_extension(kmer, kmer2tf, settings):
    """ Extender for repeat graph.
    Returns:
        new_kmer, ptf, tf, total_freq, R
    """
    if kmer2tf:
        kmer_a = "A" + kmer[:-1]
        kmer_c = "C" + kmer[:-1]
        kmer_t = "T" + kmer[:-1]
        kmer_g = "G" + kmer[:-1]
        R = [
            (kmer2tf[kmer_a]+kmer2tf[get_revcomp(kmer_a)], kmer_a),
            (kmer2tf[kmer_c]+kmer2tf[get_revcomp(kmer_c)], kmer_c),
            (kmer2tf[kmer_t]+kmer2tf[get_revcomp(kmer_t)], kmer_t),
            (kmer2tf[kmer_g]+kmer2tf[get_revcomp(kmer_g)], kmer_g),
        ]
    else:
        R = raw_reads_continue_kmer_left(kmer, settings["jf_db"])
        R.reverse()
        
    R.sort()
    total_freq = float(sum([x[0] for x in R]))
    if settings["verbose"]:
        if settings["verbose"] > 3:
            for tf, new_kmer in R:
                print "\t", new_kmer, round(tf/total_freq, 4), tf
    if not total_freq:
        yield None, None, None, None, None
    for tf, new_kmer in R:
        ptf = round(tf/total_freq, 4)
        if ptf < settings["pcutoff"]:
            continue
        if tf < settings["coverage"]:
            continue
        if settings["verbose"] == 3:
            print "\t", new_kmer, round(tf/total_freq, 4), tf
        yield new_kmer, ptf, tf, total_freq, R
    yield None, None, None, None, None


def get_repeat_graph(kmer, settings, kmer2tf):
    """
    """
    next_kmer = kmer
    seen = set()
    queue = [kmer]
    G = Graph()
    step = 0
    while queue:

        print step, len(queue), len(seen)

        seed_kmer = queue.pop(0)
        if seed_kmer in seen:
            continue
        seen.add(seed_kmer)
        for new_kmer, ptf, tf, total_freq, R in get_right_extension(seed_kmer, kmer2tf, settings):
            if not new_kmer:
                if settings["verbose"] > 3:
                    print "\tSkipped by not kmer to extend"
                break
            if settings["verbose"] == 3:
                print new_kmer, ptf, tf, total_freq, R 
            fraction_of_started = 1
            if kmer2tf:
                fraction_of_started = float(tf)/kmer2tf[kmer]
                if fraction_of_started < settings["ptf_cutoff"]:
                    if settings["verbose"]:
                        print "\tSkipped by fraction_of_started:", fraction_of_started
                    continue
            if not new_kmer in seen and not new_kmer in queue:
                queue.append(new_kmer)
            G.add_child(seed_kmer, new_kmer, ptf, tf, fraction_of_started)
        for new_kmer, ptf, tf, total_freq, R in get_left_extension(seed_kmer, kmer2tf, settings):
            if not new_kmer:
                if settings["verbose"]:
                    print "\tSkipped by not kmer to extend"
                break
            if settings["verbose"] == 3:
                print new_kmer, ptf, tf, total_freq, R 
            fraction_of_started = 1
            if kmer2tf:
                fraction_of_started = float(tf)/kmer2tf[kmer]
                if fraction_of_started < settings["ptf_cutoff"]:
                    if settings["verbose"]:
                        print "\tSkipped by fraction_of_started:", fraction_of_started
                    continue
            if not new_kmer in seen and not new_kmer in queue:
                queue.append(new_kmer)
            G.add_child(new_kmer, seed_kmer, ptf, tf, fraction_of_started) 
        G.check()
        print "queue:", len(queue)
        step += 1
    return G



def prolong(kmer, all_kmers, lyrebird, inner_name, inner_family, prev_known):
    
    if prev_known is None:
        prev_known = {}

    j = 0
    forward_candidates = set()
    rev_candidates = set()
    candidates = set()
    revkmer = get_revcomp(kmer)
    rev_candidates.add(revkmer)
    kmer2hd = {}
    for j in xrange(0, 23-17+1):
        forward_candidates = expand_kmer(kmer[j:17+j], forward_candidates, all_kmers)
    for c in forward_candidates:
        hd = hamming(kmer, c)
        if hd < 3:
            kmer2hd[c] = hd
            candidates.add(c)
    for j in xrange(0, 23-17+1):
        rev_candidates = expand_kmer(revkmer[j:17+j], rev_candidates, all_kmers)
    for c in rev_candidates:
        hd = hamming(revkmer, c)
        if hd < 3:
            kmer2hd[c] = hd
            candidates.add(c)
    print "\tEXPANDED %s for %s new kmers" % (kmer, len(candidates))
    added = 0
    candidates = [x for x in candidates if not x in prev_known]
    for new_kmer in candidates:
        new_kmer_rev = get_revcomp(new_kmer)

        prev_known[new_kmer] = (inner_name, inner_family)
        prev_known[new_kmer_rev] = (inner_name, inner_family)
        hd = kmer2hd[new_kmer]
        if new_kmer > new_kmer_rev:
            new_kmer = new_kmer_rev
        data = {
                'kid': -1,
                'kmer': new_kmer,
                'inner_name': inner_name,
                'inner_family': inner_family,
                'ref_kmer': kmer,
                'ed': hd,
                'hd': hd,
                'dust': get_dust_score(new_kmer),
                'priority': 100,
            }
        lyrebird.add_force(data, rewrite_name=True, skip_if_known=False, db=None)
        added += 1
    print "\tADDED %s for %s" % (added, kmer)
    return prev_known


if __name__ == '__main__':
    
    kmer = "AGAGTTAAACAGAGGCAAACAGA"
    jellyfish_db = "/mnt/gonduras/akomissarov/genome_white_shark/jellyfish/contigs.23.jf"
    status = "Zero"
    k = len(kmer)

    status, sequence, length, path = extend_kmer_right(kmer, jellyfish_db, tf_cutoff=0, fraq_cutoff=0.0)
    print extend_kmer_left(sequence, status, jellyfish_db, k, tf_cutoff=0, fraq_cutoff=0.0)
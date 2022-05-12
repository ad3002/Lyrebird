#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 21.05.2014
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com

from trseeker.tools.ngrams_tools import print_next_cutoff
from trseeker.tools.ngrams_tools import print_prev_cutoff


# def extend_with_consensus(kmer, kmer2tf, tf_cutoff=0, fraq_cutoff=0.0, verbose=False, lyrebird=None):

#     tf = kmer2tf[kmer]
#     consensus = kmer
#     consensus_cov = [tf]
#     instances = [kmer]

#     while True:

#         next_nucleotides = []
        
#         counter = [0,0,0,0,0]

#         for subseq in instances:
#             next_kmer = subseq[-k:]
#             R, n, nucleotides, max_hits = print_next_cutoff(next_kmer, kmer2tf, cutoff=tf_cutoff)
#             next_nucleotides.append(nucleotides)

#             for nucl, tf in R.items():
#                 if nucl == 'a':
#                     counter[0] += tf
#                 elif nucl == 'c':
#                     counter[0] += tf
#                 elif nucl == 't':
#                     counter[0] += tf
#                 else nucl == 'g':
#                     counter[0] += tf


def extend_kmer_right(kmer, kmer2tf, tf_cutoff=0, fraq_cutoff=0.0, verbose=False, lyrebird=None):
    """ Extend kmer right according to jellyfish db.

    R = [(freq, nucleotide), ...] - sorted by freq

    Iteratively extended until tf_cutoff or fraq_cutoff will be reached or TRs will be found.

    Returns: (status, sequence, length or (lengths, next_kmer), path)

    Path list:
        (norm_tf, tf, next_kmer, nucleotide)
    Status:
        Zero - reached tf_cutoff or fraq_cutoff
        TRs - TRs found
        Loop - Loop found

    """
    step = 0
    k = len(kmer)
    seen = set()
    seen.add(kmer)
    next_kmer = kmer
    sequence = kmer

    path = []

    while True:
        # get most frequent next kmer
        R, n, nucleotides, max_hits = print_next_cutoff(next_kmer, kmer2tf)
        max_hits.sort()
        sum_tf = float(sum(R.values()))
        # print "->", next_kmer, max_hits
        if sum_tf > 0:
            norm_tfs = [(round(x[0]/sum_tf, 3), x[0], next_kmer[1:]+x[1], x[1]) for x in max_hits]
        else:
            (norm_tf, tf, next_kmer, nucleotide) = (0, 0, None, None) 
            return "Zero_cov", sequence, (len(sequence)-k+1, next_kmer), path, None   
        (norm_tf, tf, next_kmer, nucleotide) = norm_tfs[-1]

        path.append(("right", step, norm_tfs))

        if lyrebird:
            known = lyrebird.query(next_kmer)
        else:
            known = None

        if known:
            print(step, known)
            # raw_input("right known?")
            return "Known", sequence, (len(sequence)-k+1, next_kmer), path, known

        step += 1
        if tf <= tf_cutoff:
            # return without last zero letter
            if verbose:
                print("Kmer with tf < %s found from right side" % tf_cutoff, sequence, len(sequence)-k+1)
            return "Zero_cov", sequence, len(sequence), path, known

        if norm_tf <= fraq_cutoff:
            # return without last zero letter
            if verbose:
                print("Kmer with norm_tf < %s found from right side" % norm_tf, sequence, len(sequence)-k+1)
            return "Zero_fraq", sequence, len(sequence), path, known

        if next_kmer in seen:
            if next_kmer == kmer:
                # return only monomer
                if verbose:
                    print("Found TRs", sequence[:-k+1], len(sequence)-k+1)
                return "TRs", sequence[:-k+1], len(sequence)-k+1, path, known
            else:
                # return sequence with last letter
                sequence += next_kmer[-1]
                if verbose:
                    print("Found loop inside path", sequence, len(sequence))
                return "Loop", sequence, (len(sequence)-k+1, next_kmer), path, known

        sequence += next_kmer[-1]

        seen.add(next_kmer)
        print("Extended right on step %s" % step, "\r", end=" ")
        

def extend_kmer_left(sequence, status, kmer2tf, k, tf_cutoff=0, fraq_cutoff=0.0, verbose=False, lyrebird=None):
    """ Extend kmer right according to jellyfish db.

    R = [(freq, nucleotide), ...] - sorted by freq

    Iteratively extended until tf_cutoff or fraq_cutoff will be reached.

    Returns: (status, sequence, length or (lengths, next_kmer), path)

    Path list:
        (norm_tf, tf, next_kmer, nucleotide)
    Status:
        Zero:prev_status - reached tf_cutoff or fraq_cutoff
        Loop:prev_status - Loop found

    """
    step = 0
    seen = []
    for i in range(0, len(sequence)-k+1):
        seen.append(sequence[i:i+k])
    next_kmer = seen[0]
    path = []

    while True:
        # get most frequence next kmer
        L, n, nucleotides, max_hits = print_prev_cutoff(next_kmer, kmer2tf)
        max_hits.sort()
        sum_tf = float(sum(L.values()))
        if sum_tf > 0:
            norm_tfs = [(round(x[0]/sum_tf, 3), x[0], x[1]+next_kmer[:-1], x[1]) for x in max_hits]
            (norm_tf, tf, next_kmer, nucleotide) = norm_tfs[-1]
        else:
            (norm_tf, tf, next_kmer, nucleotide) = (0, 0, None, None)

            return "Zero_cov", sequence, len(sequence), path, None

        path.append(("left", step, norm_tfs))

        if lyrebird:
            known = lyrebird.query(next_kmer)
        else:
            known = None

        if known:
            print(step, known)
            # raw_input("left known?")
            return("Known", sequence, (len(sequence)-k+1, next_kmer), path, known)

        step += 1
        if tf <= tf_cutoff:
            # return without last zero letter
            if verbose:
                print("Kmer with tf < %s found from left side" % tf_cutoff, sequence, len(sequence)-k+1)
            return "Zero_cov", sequence, len(sequence), path, known
        if norm_tf <= fraq_cutoff:
            # return without last zero letter
            if verbose:
                print("Kmer with norm_tf < %s found from left side" % norm_tf, sequence, len(sequence)-k+1)
            return "Zero_freq", sequence, len(sequence), path, known
        if next_kmer in seen:
            # return seqeunce with last left letter
            sequence = next_kmer[0] + sequence
            if verbose:
                print("Found loop inside path", sequence, len(sequence))
            return "Loop", sequence, (len(sequence)-k+1, next_kmer), path, known

        sequence = next_kmer[0] + sequence
    
        seen.append(next_kmer)

        print("Extended left on step %s" % step, "\r", end=" ")


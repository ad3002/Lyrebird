#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 21.05.2014
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com

from trseeker.seqio.tab_file import sc_iter_simple_tab_file
from collections import defaultdict
from trseeker.tools.sequence_tools import get_revcomp
from trseeker.seqio.fasta_file import sc_iter_fasta


def load_kmers(kmer_file, n_cutoff=None, tf_cutoff=None):
    """ Load kmers from kmer file. You can provide two cutoff:
    1) by number of kmers
    2) by minimal tf
    Returns list of kmers and dictionary kmer to tf.
    """
    all_kmers = []
    kmer2tf = defaultdict(int)
    for i, (kmer, tf) in enumerate(sc_iter_simple_tab_file(kmer_file)):
        if not kmer:
            continue
        all_kmers.append(kmer)
        kmer2tf[kmer] = int(tf)
        kmer2tf[get_revcomp(kmer)] = int(tf)
        if n_cutoff and i > n_cutoff:
            break
        if tf_cutoff and tf < tf_cutoff:
            break
    return all_kmers, kmer2tf


def load_fasta_file(fasta_file, first_n=None, size=None):
    """ Load fasta file as $-separated string.
    """
    sequences = []
    for i, seq_obj in enumerate(sc_iter_fasta(fasta_file)):
        if size:
            sequences.append(seq_obj.sequence.upper()[3000000:size])
        else:
            sequences.append(seq_obj.sequence.upper())
        if first_n is not None:
            if i == first_n-1:
                break
    return "$".join(sequences)

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
from PyExp import Timer
from trseeker.seqio.fasta_file import sc_iter_fasta
from PyBioSnippets.Lyrebird.tools.lyrebird_io import load_kmers, load_fasta_file
from collections import defaultdict, Counter
from trseeker.tools.sequence_tools import get_revcomp


def build_rs_index(sequence, k):
    n = len(sequence)
    index = {}
    kmer2tf = defaultdict(int)
    for i in xrange(0,n-k+1):
        if i % 1000000 == 0:
            print i, round(100.*i/n, 2)
        kmer = sequence[i:i+k]
        if "N" in kmer:
            continue
        if "$" in kmer:
            continue
        index.setdefault(kmer, [])
        index[kmer].append(i)
        kmer2tf[kmer] += 1
    return index, kmer2tf



if __name__ == '__main__':

    fasta_file = "/home/akomissarov/ecoli.fa"
    k = 23
    freq_limit = 1

    with Timer("Load genome"):
        sequence = load_fasta_file(fasta_file)

    with Timer("Build kmer index"):
        index, kmer2tf = build_rs_index(sequence, k)

    with Timer("Sort kmer index"):
        kmer_with_tf = [(freq,kmer) for kmer,freq in kmer2tf.items() if freq > freq_limit]
        kmer_with_tf.sort(reverse=True)

    for tf, kmer in kmer_with_tf:

        print tf, kmer

        pos = index[kmer]

        # assert len(pos) > 1
        for x in pos:
            print sequence[x-50:x+50], x

        raw_input("Next?")


    
#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 07.09.2010
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com 
"""
Functions related to assembly statistics
"""


def get_n50(lengths, ref=None):
    """ Get (N50 contig length, N50, shortest contig, longest contig) statistics for list of contigs lengths.
    @param lengths: a list of contig lengths
    @return (N50, L50, shortest contig, longest contig)
    """
    if not lengths:
        return 0, 0, 0, 0
    lengths.sort(reverse=True)
    if ref:
        total = ref
    else:
        total = sum(lengths)
    n50 = 0
    l50 = 0
    shortest_seq = min(lengths)
    longest_seq = max(lengths)
    for x in lengths:
        l50 += 1
        n50 += x
        if n50 >= total/2:
            return x, l50, shortest_seq, longest_seq


def get_ng50(lengths, genome_size, verbose=False):
    """ Get (NG50 contig length, NG50, shortest contig, longest contig) statistics for list of contigs lengths.
    @param lengths: a list of contig lengths
    @return (N50, L50, shortest contig, longest contig, NG50, LG50)
    """
    if not lengths:
        return 0, 0, 0, 0
    lengths.sort(reverse=True)
    total = sum(lengths)
    n50 = 0
    l50 = 0
    shortest_seq = min(lengths)
    longest_seq = max(lengths)
    ng50 = None
    for x in lengths:
        l50 += 1
        n50 += x
        if n50 >= genome_size/2:
            if ng50 is None:
                if verbose:
                    print("NG50: %s (%s%%)" % (x, round(100.*x/genome_size, 2)))
                    if len(lengths) > 1:
                        print("LG50: %s (%s%%)" % (l50, round(100.*(1 - float(l50-1)/(len(lengths)-1)), 2)))
                    else:
                        print("LG50: %s (%s%%)" % (l50, 100.00))
                ng50 = x
                lg50 = l50
        if n50 >= total/2:
            if verbose:
                print("N50: %s (%s%%)" % (x, round(100.*x/genome_size, 2)))
                if len(lengths) > 1:
                    print("L50: %s (%s%%)" % (l50, round(100.*(1 - float(l50-1)/(len(lengths)-1)), 2)))
                else:
                    print("LG50: %s (%s%%)" % (l50, 100.00))
                print("Contigs: ", len(lengths))
                print("Shortest:", shortest_seq)
                print("Longest: %s (%s%%)" % (longest_seq, round(100.*longest_seq/genome_size, 2)))
                print("Total length: %s (x%s)" % (sum(lengths), round(1.*sum(lengths)/genome_size, 2)))
            return x, l50, shortest_seq, longest_seq, ng50, lg50

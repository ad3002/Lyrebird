#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (C) 2007-2009 Aleksey Komissarov ( ad3002@gmail.com )
# All rights reserved.
#
# This software is licensed as described in the file COPYING, which
# you should have received as part of this distribution.
"""
Short functions for working with collections of raw reads.
"""
from trseeker.seqio.sra_file import fastq_reader
from trseeker.tools.jellyfish_tools import query_kmers, query_kmers_new
from collections import Counter, defaultdict
from trseeker.tools.edit_distance import get_ed_similarity


def raw_reads_load_fastq_reads(fastq_file, limit=1000000):
    ''' Return list of raw reads from fastq file.

    @param fastq_file: fastq file with reads.
    @param limit: maximum number of reads
    '''
    result = []
    for i, read in enumerate(fastq_reader(fastq_file)):
        print(i, "\r", end=" ")
        result.append(read.sequence)
        if i > limit:
            break
    print
    return result


def get_reads_by_kmer(kmer, all_reads):
    ''' Return reads with given kmer.

    @param kmer: kmer
    @param all_reads: list of reads
    '''
    result = []
    for read in all_reads:
        if kmer in read:
            result.append(read)
    return result


def raw_reads_get_flanks(kmer, reads):
    ''' Get kmer flanks with raw reads.

    @param kmer: given kmer
    @param reads: list of reads with kmer
    @return: [(left_flank, right_flank),...], unsplitted_reads
    '''
    k = len(kmer)
    result = []
    errors = []
    for read in reads:
        try:
            left, right = read.split(kmer)
        except:
            errors.append(read)
        result.append((left, right))
    return result, errors


def raw_reads_get_next_kmers(right_flanks, k, i=0):
    ''' Return next kmers for right flanks

    @param right_flanks: list of flanks
    @param k: kmer length
    @param i: start position
    '''
    result = []
    for flank in right_flanks:
        if len(flank)-i < k:
            continue
        result.append(flank[i:k])
    return result


def raw_reads_get_shifted_kmers(right_flanks, kmer, i):
    ''' Return list of shifted kmers by i bases from right_flanks.

    @param right_flanks: list of flanks
    @param kmer: given kmer
    @param i: shift value
    '''
    result = []
    k = len(kmer)
    errors = 0
    if i < k:
        prefix = kmer[i:]
        j = 0
        l = k - len(prefix)
    else:
        prefix = None
        j = i - k
        l = k
    for flank in right_flanks:
        if prefix:
            data = prefix+flank[j:j+l]
        else:
            data = flank[j:j+l]
        if len(data) == k:
            result.append(data)
        else:
            errors += 1
    return result, errors


def raw_reads_get_variants(kmer, case="upper"):
    ''' Return left and rigth shift variants for given kmer.

    @param kmer: given kmer
    @param case: case upper (default) or lower
    '''
    if case == "upper":
    	alphabet = ["A", "C", "G", "T"]
    else:
    	alphabet = ["a", "c", "g", "t"]
    kmer = kmer.strip()
    left_kmer = kmer[1:]
    right_kmer = kmer[:-1]
    left_data = []
    right_data = []
    for letter in alphabet:
        right_data.append(left_kmer+letter)
        left_data.append(letter+right_kmer)
    return left_data, right_data


def raw_reads_continue_kmer_right(kmer, jellyfish_db, cutoff=0):
    ''' Return next kmers for given kmer according to jellyfish database.

    @param kmer: given kmer
    @param jellyfish_db: jellyfish database
    '''
    left_data, right_data = raw_reads_get_variants(kmer)
    if not isinstance(jellyfish_db, defaultdict):
        R = query_kmers(jellyfish_db, right_data, both_strands=True, verbose=False)
        R = [(int(v),k) for k,v in R.items()]
    else:
        R = [(int(jellyfish_db[_kmer]), _kmer) for _kmer in right_data]
    if cutoff:
        for i, x in enumerate(R):
            if x[0] < cutoff:
                R[i] = (0, x[1])
    R.sort()    
    return R


def raw_reads_continue_kmer_left(kmer, jellyfish_db, cutoff=0):
    ''' Return previous kmers for given kmer according to jellyfish database.

    @param kmer: given kmer
    @param jellyfish_db: jellyfish database
    '''
    left_data, right_data = raw_reads_get_variants(kmer)
    if not isinstance(jellyfish_db, defaultdict):
        L = query_kmers(jellyfish_db, left_data, both_strands=True, verbose=False)
        L = [(int(v),k) for k,v in L.items()]
    else:
        L = [(int(jellyfish_db[_kmer]), _kmer) for _kmer in left_data]
    L.sort()
    if cutoff:
        for i, x in enumerate(L):
            if x[0] < cutoff:
                L[i] = (0, x[1])
    return L


def raw_reads_get_variability(right_flanks, k, used_kmers, ed_cutoff=80):
    ''' Return variability of shifted kmers.

    @param right_flanks: list of flanks
    @param k: k
    @param used_kmers: set of previously used kmers
    @param ed_cutoff: cutoff of edit distance
    @return: (fraction_of_ok, ref_kmer, used_kmers, new_kmers)
    '''
    next_kmers = raw_reads_get_next_kmers(right_flanks, k)
    c = Counter(next_kmers)
    ref_kmer = c.most_common()[0][0]
    new_kmers = set()

    ok = 0
    error = 0
    for x in c:
        ed = get_ed_similarity(ref_kmer, x)
        if ed > ed_cutoff:
            ok += c[x]
            used_kmers.add(x)
            used_kmers.add(get_revcomp(x))
            new_kmers.add(x)
        else:
            error += c[x]
    return float(ok)/(error + ok), ref_kmer, used_kmers, new_kmers


def group_kmers_by_hamming(kmers, d=1):
    ''' Return group of kmer joined by hamming distance.

    @param kmers: list of kmers
    @param d: hamming distance
    '''
    pass
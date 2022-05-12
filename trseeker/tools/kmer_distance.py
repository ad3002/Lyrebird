#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (C) 2007-2013 Aleksey Komissarov ( ad3002@gmail.com )
# All rights reserved.
#
# This software is licensed as described in the file COPYING, which
# you should have received as part of this distribution.
'''
For input list of sequence compute distance between them by similar kmers.

Parameters:
k - kmer length
'''
import os
from trseeker.seqio.tab_file import sc_iter_tab_file
from trseeker.models.trf_model import TRModel
from trseeker.tools.ngrams_tools import *
from trseeker.settings import NGRAM_N, NGRAM_LENGTH
from collections import defaultdict
from trseeker.seqio.tr_file import get_all_trf_objs

class KmerBasedDistance(object):

    SKIP_N = True
    verbose = False

    def __init__(self, data, k=NGRAM_LENGTH, verbose=False):
        self.data = data
        self.kmers = []
        self.kmer2p = {}
        self.pointer = -1
        self.k = k
        self.tf = defaultdict(int)
        self.df = defaultdict(int)
        self.p2docids = defaultdict(list)
        self.p2doctfs = defaultdict(list)
        self.p2docnfs = defaultdict(list)
        self.N = len(self.data)
        self.D = defaultdict(float)
        if self.N>100 or verbose:
            self.verbose = True
       

    def count_kmers(self, docid):
        ''' Update tf and df data with k-mers from given sequence.
        '''
        # import cProfile, pstats, io
        # pr = cProfile.Profile()
        # pr.enable()

        sequence = self.data[docid]
        _local_tf = defaultdict(int)
        self.local_nf = {}
        self.local_tf = defaultdict(int)
        local_seen = {}
        n = float(len(sequence))
        for i in range(0, len(sequence) - self.k + 1):
            kmer = sequence[i:i + self.k]
            if self.SKIP_N and 'n' in kmer:
                continue
            _local_tf[kmer] += 1
        for kmer in _local_tf:
            
            if not kmer in self.kmers:
                self.pointer += 1
                self.kmer2p[kmer] = self.pointer
                self.kmers.append(kmer)
            p = self.kmer2p[kmer]
            if not p in local_seen:
                self.df[p] += 1
                local_seen[p] = True
            self.tf[p] += _local_tf[kmer]
            self.local_tf[p] = _local_tf[kmer]
        for p in self.local_tf:
            self.local_nf[p] = self.local_tf[p]/n

        # pr.disable()
        # pr.print_stats()
        

    def count_tf_df_for_data(self):
        '''
        '''
        print "Count tf/df for data"
        for docid in xrange(self.N):

            import cProfile, pstats, io
            pr = cProfile.Profile(timeunit=1)
            pr.enable()
            

            
            
            #self.count_kmers(docid)

            sequence = self.data[docid]

            if self.verbose:
                print "Process tf/df: ", docid, self.N, len(sequence)
            
            _local_tf = defaultdict(int)
            local_nf = {}
            local_tf = defaultdict(int)
            local_seen = {}
            n = float(len(sequence))
            for i in xrange(0, len(sequence) - self.k + 1):
                kmer = sequence[i:i + self.k]
                if self.SKIP_N and 'n' in kmer:
                    continue
                _local_tf[kmer] += 1
            for kmer in _local_tf:
                
                if not kmer in self.kmers:
                    self.pointer += 1
                    self.kmer2p[kmer] = self.pointer
                    self.kmers.append(kmer)
                p = self.kmer2p[kmer]
                if not p in local_seen:
                    self.df[p] += 1
                    local_seen[p] = True
                self.tf[p] += _local_tf[kmer]
                local_tf[p] = _local_tf[kmer]
            for p in local_tf:
                local_nf[p] = local_tf[p]/n

            for p in local_tf:
                self.p2docids[p].append(docid)
                self.p2doctfs[p].append(local_tf[p])
                self.p2docnfs[p].append(local_nf[p])
            pr.disable()
            pr.print_stats()
        
        print

    def compute_distances(self):
        '''
        '''
        print "Compute distances"
        kmer_n = len(self.df)
        for m, p in enumerate(self.df):
            print m, kmer_n, "\r",
            if self.df[p] == 1:
                continue
            docids = self.p2docids[p]
            for i, docid_a in enumerate(docids):
                for j, docid_b in enumerate(docids[i+1:]):
                    if docid_a < docid_b:
                        key = "%s\t%s" % (docid_a, docid_b)
                    else:
                        key = "%s\t%s" % (docid_b, docid_a)
                    self.D[key] += min(
                                self.p2docnfs[p][i],
                                self.p2docnfs[p][j]
                    )
        print

    def save_index_to_file(self, file_name):
        '''
        '''
        print "Save index data"
        with open(file_name, "w") as fh:
            for p, kmer in enumerate(self.kmers):
                d = (kmer, 
                         get_revcomp(kmer),
                         self.tf[p],
                         self.df[p],
                         ",".join(map(str, self.p2docids[p])),
                         ",".join(map(str, self.p2doctfs[p])),
                         ",".join(map(str, self.p2docnfs[p])),
                    )
                d = "\t".join(map(str, d)) 
                fh.write("%s\n" % d)

    def save_distances_to_file(self, file_name):
        '''
        '''
        print "Save distance data"
        with open(file_name, "w") as fh:
            for key, value in self.D.items():
                fh.write("%s\t%s\n" % (key, value))

def compute_kmer_profiles_for_trs(trf_large_file, output_folder, k):
    '''
    '''
    if not os.path.isdir(output_folder):
        os.makedirs(output_folder)
    for i, trf_obj in enumerate(sc_iter_tab_file(trf_large_file, TRModel)):
        print "Compute for", i, "\r",
        file_name = os.path.join(output_folder, "%s.ngram" % str(trf_obj.trf_id))
        ngrams = get_ngrams_freq(trf_obj.trf_array, m=NGRAM_N, n=k)
        with open(file_name, "w") as fh:
            for (ngram, tf, nf) in ngrams:
                rngram = get_revcomp(ngram)
                if ngram < rngram:
                    data = "%s\t%s\t%s\t%s\n" % (ngram, rngram, tf, nf)
                else:
                    data = "%s\t%s\t%s\t%s\n" % (rngram, ngram, tf, nf)
                fh.write(data)
    print


def get_dust_score(sequence, k=4):
    ''' Return DUST score for given sequence and kmer length.
    '''
    d = defaultdict(int)
    for i in range(0, len(sequence)-k+1):
        d[sequence[i:i+k]] += 1
    score = 0.
    total = 0.
    for v in d.values():
        score += v*(v-1)/2.
        total += score
    return total/(len(sequence)-k+1)



def compile_ngrams(trf_large_file, ngram_index_file, k=NGRAM_LENGTH, cutoff=None, dust=False):
    """ Compile ngrams collection for given project."""

    data = []
    trf_index = []
    print "Read arrays..."
    trf_index = get_all_trf_objs(trf_large_file)
    trf_index = [trf_obj for trf_obj in trf_index if trf_obj.trf_array_length >= k]
    data = [trf_obj.trf_array for trf_obj in trf_index]
    print "Skipped %s arrays as short" % (len(trf_index) - len(data))
    print "Process kmers..."
    index_data = process_list_to_kmer_index(data, k, docids=True, cutoff=cutoff)
    _process_index_data_to_file(ngram_index_file, index_data, k, dust=dust, trf_index=trf_index)

def _process_index_data_to_file(ngram_index_file, index_data, k, dust=False, trf_index=None):
    print "Sort data..."
    result = []
    skipped_by_dust = 0
    for i, (key, revkey, tf, df, docids, freqs) in enumerate(index_data):
        if dust:
            if get_dust_score(key) > dust:
                skipped_by_dust += 1
                continue
        new_doc_ids = []
        for j in docids:
            new_doc_ids.append(trf_index[j].trf_id)
        # sort docsis by freqs
        items = []
        for pos in range(len(new_doc_ids)):
            all_items = trf_index[pos].trf_array_length - k + 1
            if all_items <= 0:
                continue
            items.append((
                    freqs[pos] * 1. / all_items,
                    new_doc_ids[pos]
                ))
        items.sort()
        new_doc_ids = [str(x[1]) for x in items]
        freqs = [str(x[0]) for x in items]

        data = [key, revkey, tf, df, ",".join(new_doc_ids), ",".join(freqs)]
        data = "%s\n" % "\t".join(map(str, data))
        result.append(data)
    if dust:
        print "Skipped by dust:", skipped_by_dust
    print "Save data to %s..." % ngram_index_file
    with open(ngram_index_file, "w") as fh:
        fh.writelines(result)

def compile_kmer_index_from_arrays(arrays, ngram_index_file, k=NGRAM_LENGTH):
    """ Compile ngrams collection for given project."""
    print "Process kmers..." 
    index_data = process_list_to_kmer_index(arrays, k, docids=True)
    print "Sort data..."
    result = []
    for i, (key, revkey, tf, df, docids, freqs) in enumerate(index_data):
        new_doc_ids = []
        for j in docids:
            new_doc_ids.append(j)
        data = [key, revkey, tf, df, ",".join(map(str, new_doc_ids)), ",".join(map(str, freqs))]
        data = "%s\n" % "\t".join(map(str, data))
        result.append(data)
    print "Save data..."
    with open(ngram_index_file, "w") as fh:
        fh.writelines(result)

def compile_kmer_index_from_kmer_profiles():
    pass

def compute_distances(index_objs, id2length, index_function=None):
    ''' Compute distance between kmer profiles.
    '''
    print "Compute distances..."
    n = len(index_objs)
    D = {}
    KEYS_CUTOFF = 1000000
    chunk = 0
    for k, kmer_index in enumerate(index_objs):
        ids = kmer_index.docs
        print "Current size", len(ids)
        for i, trid_a in enumerate(ids):
            print k, n, i, len(ids), "\r",
            for j, trid_b in enumerate(ids[i+1:]):
                if trid_a < trid_b:
                    key = "%s\t%s" % (trid_a, trid_b)
                else:
                    key = "%s\t%s" % (trid_b, trid_a)
                D.setdefault(key, 0.0)
                D[key] += min(
                            kmer_index.freqs[i]*1./id2length[trid_a],
                            kmer_index.freqs[j]*1./id2length[trid_b]
                )
        if len(D) > KEYS_CUTOFF:
            print "save D for chunk", chunk
            file_name = "chunk%s.dat" % chunk
            chunk += 1
            with open(file_name, "w") as fh:
                for key in D:
                    s = "%s\t%s\n" % (key, D[key])
                    fh.write(s)
            D = {}
        index_objs[k] = None
    return D

def compute_distances_for_index(index, id2length, index_function=None):
    ''' Compute distance between kmer profiles.
    '''
    print "Compute distances..."
    n = len(index)
    D = {}
    for k, kmer_index in enumerate(index):
        ids = kmer_index[4]
        for i, trid_a in enumerate(ids):
            print k, n, i, len(ids), "\r",
            for j, trid_b in enumerate(ids[i+1:]):
                if trid_a < trid_b:
                    key = "%s\t%s" % (trid_a, trid_b)
                else:
                    key = "%s\t%s" % (trid_b, trid_a)
                D.setdefault(key, 0.0)
                D[key] += min(
                            kmer_index[5][i]*1./id2length[trid_a],
                            kmer_index[5][j]*1./id2length[trid_b]
                )
        index[k] = None
    return D

def compute_distances_for_index_by_raw_kmers(index, id2length, index_function=None):
    ''' Compute distance between kmer profiles.
    '''
    print "Compute distances..."
    n = len(index)
    D = {}
    for k, kmer_index in enumerate(index):
        ids = kmer_index[4]
        for i, trid_a in enumerate(ids):
            print k, n, i, len(ids), "\r",
            for j, trid_b in enumerate(ids[i+1:]):
                if trid_a < trid_b:
                    key = "%s\t%s" % (trid_a, trid_b)
                else:
                    key = "%s\t%s" % (trid_b, trid_a)
                D.setdefault(key, 0.0)
                D[key] += min(
                            kmer_index[5][i],
                            kmer_index[5][j]
                )
        index[k] = None
    return D
#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 08.09.2011
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com 
'''
Working with ngram indexes.
'''

def save_ngram_index(ngram_index_file, hash2id, result_tf,
                     result_df, result_rtf, result_rdf,
                     seen_rev, hash2rev, ngram2kmerann):
    ''' Write ngram index. 
    
    Format: id, rev_id, hash, rev_hash, tf, df, etf, edf.

    etf and edf are tf and df taking into account reverse complement.'''
    
    with open(ngram_index_file, "w") as fh:
        for kmer, id in hash2id.items():
            if not kmer in seen_rev:
                rev_id = hash2id[hash2rev[kmer]]

                if type(ngram2kmerann) == dict:
                    kmer_ann = [str(x) for x in ngram2kmerann[kmer]]
                    kmer_ann = ",".join(kmer_ann)
                else:
                    kmer_ann = ""

                data = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
                                 id,
                                 rev_id,
                                 kmer,
                                 hash2rev[kmer],
                                 result_tf[id],
                                 result_df[id],
                                 result_rtf[id],
                                 result_rdf[id],
                                 kmer_ann)
                fh.write(data)

def save_ngram_pos_index(ngram_trids_file, id2trids, id2trid2tf):
    ''' Write ngrams position data. 
    Format: id, comma separated list of trids,
    comma separated list of tf in trid.
    '''
    with open(ngram_trids_file, "w") as fh:
        for id, trids in id2trids.items():
            trids = list(trids)
            tfs = []
            for trid in trids:
                tfs.append(id2trid2tf[id][trid])
            trids = ",".join([str(x) for x in trids])
            tfs = ",".join([str(x) for x in tfs])

            data = "%s\t%s\t%s\n" % (id, trids, tfs)
            fh.write(data)

def save_distance_data(dist_file, distances):
    ''' Save distance data: a b dist'''

    print "Sort distances"
    distances = [(key, d) for key, d in distances.items()]
    distances.sort(reverse=True, key=lambda x: x[1])

    print "Save result"
    with open(dist_file, "w") as fh:
        for key, value in distances:
            fh.write("%s\t%s\n" % (key, value))
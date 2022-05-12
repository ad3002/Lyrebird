#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 30.06.2011
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com 

from PyExp import AbstractModel
from trseeker.seqio.tab_file import sc_iter_tab_file, sc_iter_simple_tab_file
from trseeker.tools.sequence_tools import get_dust_score

class KmerIndexModel(AbstractModel):
    ''' Model for ngram to TR data.

    '''
    dumpable_attributes = [
                           "kmer",
                           "rkmer",
                           "tf",
                           "df",
                           "docs",
                           "freqs"
                           ]
    int_attributes = [
                      "tf",
                      "df",
                      ]

    float_attributes = []

    list_attributes = ["docs",
                       "freqs",
                       ]

    list_attributes_types = {"docs":int,
                             "freqs":float,
                             }

class KmerSliceModel(AbstractModel):
    '''
    '''
    dumpable_attributes = ["kmer",
                             "local_tf",
                             "df",
                             "f_trs",
                             "f_wgs",
                             "f_mwgs",
                             "f_trace",
                             "f_sra",
                             "f_ref1",
                             "f_ref2",
                             "f_caroli",
                             ]
    int_attributes = ["local_tf",
                             "df",
                        ]

    def __init__(self, data):
      super(KmerSliceModel, self).__init__()
      data = data.split(":")
      if len(data) < len(self.dumpable_attributes):
        for i in xrange(len(self.dumpable_attributes) - len(data)):
          data.append("-1")
      for i, attr in enumerate(self.dumpable_attributes):
        setattr(self, attr, data[i])
      self.local_tf = int(self.local_tf)
      self.df = int(self.df)

    def set_with_dict(self, dictionary):
      raise NotImplemented 

    def set_with_list(self, dictionary):
      raise NotImplemented

    def __str__(self):
      return ":".join([str(getattr(self, x)) for x in self.dumpable_attributes])


class SliceTreeModel(AbstractModel):
    ''' Model for kmers tree node.
    '''

    dumpable_attributes = [
                           "deep",
                           "size",
                           "blast_fams",
                           "maxdf",
                           "nmaxdf",
                           "pmaxdf",
                           "gc_var",
                           "units",
                           "trs",
                           "kmers",
                           ]

    int_attributes = [
                           "deep",
                           "size",
                           "maxdf",
                           "nmaxdf",
                      ]

    float_attributes = [
                        "pmaxdf",
                        "gc_var",
    ]

    list_attributes = ["blast_fams",
                       "units",
                       "trs",
                       "kmers",
                       ]

    list_attributes_types = {"blast_fams":str,
                             "units":int,
                             "trs":int,
                             "kmers":KmerSliceModel,
                             }

def sc_ngram_reader(file_name):
    ''' Read file with ngrams data.'''
    for obj in sc_iter_tab_file(file_name, KmerIndexModel):
        yield obj

def sc_ngram_trid_reader(file_name):
    ''' Read file with ngram to TRs+tfs data.'''
    for nid, sid_list, tf_list  in sc_iter_simple_tab_file(file_name):
        sids = sid_list.split(",")
        tfs = tf_list.split(",")

        result = []
        for i in xrange(0, len(sids)):
            result.append((int(sids[i]), float(tfs[i])))
        yield int(nid), result

def read_kmer_index(ngram_index_file, micro_kmers, cutoff=1, dust=0):
    ''' Read kmer index data as list
    '''
    data = []
    for index_obj in sc_iter_tab_file(ngram_index_file, KmerIndexModel):
        print index_obj.df, "\r",
        if micro_kmers:
            if index_obj.kmer in micro_kmers:
                if index_obj.tf < micro_kmers[index_obj.kmer]:
                    continue
        if dust:
          if get_dust_score(index_obj.kmer) < dust:
            continue
        if index_obj.df == cutoff:
            break
        data.append(index_obj)
    print 
    return data
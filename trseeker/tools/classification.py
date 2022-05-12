#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 07.09.2010
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com
"""

"""

from collections import defaultdict
from collections import Counter
import math
from trseeker.tools.ngrams_tools import process_list_to_kmer_index


class KmerClassificator(object):
    """

    """

    def __init__(self, kmer_dict, k, kmer2df, n, fam2n, fam2ntf):
        self.kmer2data = kmer_dict
        self.k = k
        self.kmer2df = kmer2df
        self.N = float(n)
        self.fam2N = fam2n
        self.fam2NTF = fam2ntf

    def add_gc(self, family2gc):
        self.fam2gc = family2gc

    def classify_by_array(self, array):
        """
        """
        n = len(array)
        fam2p = Counter()
        for i in xrange(n-self.k+1):
            kmer = array[i:i+self.k]
            families = self.kmer2data[kmer]
            for f in families:
                fam2p[f[0]] += f[1] * f[2]
        return fam2p.items()

    def cosine_distance(self, array):
        """
        """
        n = len(array)
        fam2numerators = defaultdict(list)
        fam2denominators = defaultdict(list)

        kmer2tf = defaultdict(int)
        idf = {}
        for i in xrange(n-self.k+1):
            kmer = array[i:i+self.k]
            kmer2tf[kmer] += 1
        for kmer in kmer2tf:
            # print self.N, 
            if not kmer in self.kmer2df:
                self.kmer2df[kmer] = 1
            idf[kmer] = math.log( (1+self.kmer2df[kmer]) / (self.N+1))
            kmer2tf[kmer] = (kmer2tf[kmer]/(n-self.k+1)) * idf[kmer]
        # print kmer2tf
        # print idf

        for kmer in kmer2tf:
            families = self.kmer2data[kmer]
            for (f, tf, df) in families:
                #print f, tf, df, idf[kmer], kmer2tf[kmer]
                w = tf * idf[kmer]
                fam2numerators[f].append(
                            w * kmer2tf[kmer]
                    )
                fam2denominators[f].append(
                            (kmer2tf[kmer]**2, w**2)
                    )
        D = []
        # print fam2denominators
        # print fam2numerators
        for f in fam2numerators:
            n1 = sum([x[0] for x in fam2denominators[f]])
            n2 = sum([x[1] for x in fam2denominators[f]])
            denom = math.sqrt(n1) * math.sqrt(n2)
            if denom:
                D.append(
                        (f, sum(fam2numerators[f])/denom)
                    )
        return D

    def cosine_distance2(self, array, array_gc, cutoff=0.001):
        """
        10000 *
            qtf/|qkmers| *
                [1+dtf] / [1+|dkmers|] *
                    [1+df] / [1+|docs|] *
                        1 - |arrayGC-familyGC|

        """
        n = len(array)
        fam2denominators = defaultdict(list)

        kmer2tf = defaultdict(float)
        idf = {}
        all_tf = float(n - self.k + 1)
        for i in xrange(n-self.k+1):
            kmer = array[i:i+self.k]
            kmer2tf[kmer] += 1.

        for kmer in kmer2tf:
            kmer2tf[kmer] = kmer2tf[kmer]/all_tf
            
        for kmer in kmer2tf:
            if not kmer in self.kmer2data:
                continue
            families = self.kmer2data[kmer]
            for (f, tf, df) in families:
                #print f, tf, df, idf[kmer], kmer2tf[kmer]
                fam2denominators[f].append(
                            10000*kmer2tf[kmer] * ((1+float(tf))/(1+self.fam2NTF[f])) * ((1+float(df))/(1+self.fam2N[f])) * (1-abs(self.fam2gc[f] - array_gc))
                    )
        D = []
        # print fam2denominators
        # print fam2numerators
        for f in fam2denominators:
            value = sum([x for x in fam2denominators[f]])
            if value > cutoff:
                D.append(
                        (f, value)
                )
        return D


def compute_classification_dictionary(fam2trs, k, docs):
    """
    kmer -> [(family, tf, df),]
    family -> (N, sum(|n|) - kN + N)
    """
    fam2index = {}
    fam2N = {}
    fam2L = {}
    result = defaultdict(list)
    # kmer2tf = defaultdict(int)
    kmer2df = defaultdict(int)

    fam2N = {}
    fam2NTF = {}
    for fid, trs_objs in fam2trs.items():
        arrays = [x.trf_array for x in trs_objs]
        N = len(arrays)
        NTF = 0
        L = sum([len(x) for x in arrays]) - k*N + N
        index = process_list_to_kmer_index(arrays, k, docids=False)
        for data in index:
            kmer = data[0]
            rkmer = data[1]
            tf = data[2]
            df = data[3]
            NTF += tf
            # kmer2tf[kmer] += tf
            kmer2df[kmer] += df
            # kmer2tf[rkmer] += tf
            kmer2df[rkmer] += df
            result[kmer].append((fid, tf, df))
        fam2NTF[fid] = NTF
        fam2N[fid] = N
    classificator = KmerClassificator(result, k, kmer2df, docs, fam2N, fam2NTF)
    return classificator
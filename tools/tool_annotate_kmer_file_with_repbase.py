#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 16.01.2019
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com 

from collections import defaultdict
import sys, os

humrep1 = "/mnt/voronezh/pipeline_ak/RepBase24.01.fasta/humrep.ref"
humrep2 = "/mnt/voronezh/pipeline_ak/RepBase24.01.fasta/appendix/humapp.ref"

kmer2repeat = defaultdict(list)
k = 23


def sc_iter_fasta_brute(file_name, inmem=False):
    """ Iter over fasta file."""
    
    header = None
    seq = []
    with open(file_name) as fh:
        if inmem:
            data = fh.readlines()
        else:
            data = fh
        for line in data:
            if line.startswith(">"):
                if seq:
                    sequence = "".join(seq)
                    yield (header, sequence)
                header = line
                seq = []
                continue
            seq.append(line)
        if seq or header:
            sequence = "".join(seq)
            yield (header, sequence)
            

for file_name in (humrep1, humrep2):
    for (header, sequence) in sc_iter_fasta_brute(file_name):
        header = header[1:].split()[0]
        sequence = sequence.upper()
        for i in xrange(len(sequence)-k+1):
            kmer = sequence[i:i+k]
            kmer2repeat[kmer].append((header,i))

print "Reading TRs..."
trf_all_file = "/mnt/voronezh/pipeline_ak/human_38.trf_all.trf"
TRs = []
with open(trf_all_file) as fh:
    for line in fh:
        d = line.strip().split("\t")
        TRs.append(d)

if __name__ == '__main__':

    
    kmer_file = sys.argv[1]
    output_file = sys.argv[2]

    with open(output_file, "w") as fw:
        with open(kmer_file) as fh:
            for i, line in enumerate(fh):
                if not line.strip():
                    continue
                kmer, rkmer, df, tf, hits, gcs = line.upper().strip().split("\t")

                # print i, df, tf, kmer2repeat[kmer], kmer2repeat[rkmer]

                d = "%s\t%s\t%s\n" % (line.strip(), kmer2repeat[kmer], kmer2repeat[rkmer])
                fw.write(d)

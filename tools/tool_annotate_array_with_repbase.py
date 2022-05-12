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



def get_revcomp(sequence):
    '''Return reverse complementary sequence.

    >>> complementary('AT CG')
    'CGAT'

    '''
    c = dict(zip('ATCGNatcgn~[]', 'TAGCNtagcn~]['))
    return ''.join(c.get(nucleotide, '') for nucleotide in reversed(sequence))


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
            
print "Reading repeats db..."
for file_name in (humrep1, humrep2):
    for (header, sequence) in sc_iter_fasta_brute(file_name):
        header = header[1:].split()[0]
        sequence = sequence.upper()
        for i in xrange(len(sequence)-k+1):
            kmer = sequence[i:i+k]
            kmer2repeat[kmer].append((header,i))
            rkmer = get_revcomp(kmer)
            kmer2repeat[rkmer].append((header,i))

print "Load centromeres"
file_name = "/home/genomerussia/resources/b38/bwa/full/20150713_location_of_centromeres_and_other_regions.txt"
incen = {}
with open(file_name) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        d = line.strip().split("\t")
        if not d:
            continue
        name, chrid, start, end, tfeature, _, _, _ = d
        incen["CHR%s" % chrid] = (int(start), int(end))


print incen

file_name = "/home/genomerussia/resources/b38/bwa/full/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai"
lengths = {}
with open(file_name) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        d = line.strip().split("\t")
        if not d:
            continue
        name, leng, _, _, _ = d
        lengths[name.upper()] = int(leng)



if __name__ == '__main__':

    
    trf_all_file = sys.argv[1]
    output_file = sys.argv[2]

    print "Reading TRs..."
    TRs = []
    with open(trf_all_file) as fh:
        for line in fh:
            d = line.upper().strip().split("\t")
            TRs.append(d)

    with open(output_file, "w") as fw:
        for rid, d in enumerate(TRs):

            d[6] = int(d[6])
            d[7] = int(d[7])
            
            trf_array = d[14]
            trf_monomer = d[13]
            match = defaultdict(int)
            ann_cov = 0
            for i in xrange(len(trf_array)-k+1):
                kmer = trf_array[i:i+k]
                known = False
                kmer2repeat[kmer].sort()
                for repeat,pos in kmer2repeat[kmer]:
                    match[repeat] += 1
                    known = True
                if known:
                    ann_cov += 1
            match = match.items()
            match.sort(key=lambda x: x[1], reverse=True)
            chrm = d[18].split()[0].upper()
            if "_" in chrm or "decoy" in chrm:
                d1 = "ChrUN"
            elif chrm in incen and abs(incen[chrm][0] - d[6]) < 5000000 or abs(incen[chrm][1] - d[6]) < 5000000:
                d1 = ("CEN", abs(incen[chrm][0] - d[6])/1000000, abs(incen[chrm][1] - d[6])/1000000)
            elif d[6] < 5000000:
                d1 = ("START", d[6]/1000000)
            elif chrm in lengths and lengths[chrm] - d[6] < 5000000:
                d1 = ("END", (lengths[chrm] - d[6])/1000000)
            else:
                d1 = None
                
            r = d[8], d[9], d[10], i, ann_cov, round(100.*ann_cov/(len(trf_array)-k+1), 2), match[:5], "%s: %s-%s" % (chrm, d[6], d[7]), d1
            print r
            fw.write("%s\n" % map(str,r))
            

            # d = "%s\t%s\t%s\n" % (line.strip(), kmer2repeat[kmer], kmer2repeat[rkmer])
            # fw.write(d)

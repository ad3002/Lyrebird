#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 09.06.2010
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com 
'''
Functions related to short read data (SRA).
'''
import os
from trseeker.tools.other_tools import sort_dictionary_by_value
from trseeker.tools.ngrams_tools import generate_ngrams

def sra_fastaq_reader(file_name):
    ''' Iterate over fastaq data.'''
    with open(file_name, "rb") as fh:
        i = 0
        for line in fh:
            i += 1
            if i == 1:
                title = line.split()[0]
            elif i == 2:
                seq = line.strip()
                yield title, seq.lower()
            elif i == 4:
                i = 0
            else:
                continue

def sra_fasta_reader(file_name):
    ''' Iterate over fasta SRA data.'''
    seq = None
    with open(file_name, "rb") as fh:
        title = None
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if title:
                    yield title, seq
                title = line
                seq = ""
                continue
            else:
                seq += line
        if title:
            yield title, line

def read_fastaq_freq_to_memory(file_name):
    ''' Read fastaq data with correspoding TF.'''
    SRA = {}
    for title, seq in sra_fastaq_reader(file_name):
        SRA.setdefault(seq, 0)
        SRA[seq] += 1
    return SRA

def fastaq_to_fasta(file_name, output_file):
    ''' Convert fastaq to fasta.'''
    n = 0
    with open(output_file, "w") as fh:
        for title, seq in sra_fastaq_reader(file_name):
            s = ">%s\n%s\n" % (title, seq)
            fh.write(s)
            n += 1
    print n

def write_fastaq_repeats(input_file, output_file, min_tf=1000):
    ''' Write file with fastaq repeated reads with exact match.'''
    SRA = read_fastaq_freq_to_memory(input_file)
    i = 0
    with open(output_file, "w") as fw:
        for k in SRA:
            if SRA[k] >= min_tf:
                i += 1
                fw.write(">%s-%s\n%s\n" % (i, SRA[k], k))

def seq_to_bin(seq):
    ''' Write binary representation of sequence.
    TODO: fixit
    '''
    LBT = {'a':'00',
     'c':'01',
     'g':'10',
     't':'11',
     'A':'00',
     'C':'01',
     'G':'10',
     'T':'11', }

    b = '0b1'
    for l in seq:
        b += LBT[l]
    return eval(b)

def bin_to_seq(bseq):
    ''' Write string representation of binary sequence data.
    TODO: fixit
    '''

    BLT = {'00' : 'a',
       '01':'c',
       '10':'g',
       '11':'t', }

    bseq = str(bin(bseq))[3:]
    result = ''
    for i in xrange(0, len(bseq), 2):
        result += BLT[bseq[i:i + 2]]
    return result

def write_reduced_fasta(input_file, output_reduced):
    ''' Read SRA fasta data without read repeats.
    Output format: read tf
    Return number of skipped reads.
    '''
    seen = {}
    n = 0
    print "Reading..."
    skipped = 0
    for title, seq in sra_fasta_reader(input_file):
        if 'N' in seq:
            skipped += 1
            continue

        seq = seq_to_bin(seq)
        seen.setdefault(seq, 0)
        seen[seq] += 1
        n += 1
    print "Sorting..."
    seen = sort_dictionary_by_value(seen, reverse=True)
    print "Formating..."
    seen = ["%s\t%s\n" % (k, v) for v, k in seen]
    print "Saving..."
    with open(output_reduced, "w") as fh:
        for line in seen:
            fh.write(line)
    print "Done."
    print "Skipped: ", skipped

    return skipped

def write_ngrams(input_file, output_ngram, NGRAM_N):
    ''' Write ngrams data from SRA fasta data.
    Output format: binary_ngram freq
    Return: Nngram, skipped read
    '''
    seen = {}
    n = 0
    print "Reading..."
    skipped = 0
    for title, seq in sra_fasta_reader(input_file):
        for i, ngram in generate_ngrams(seq, n=NGRAM_N):
            if 'N' in ngram:
                skipped += 1
                continue
            ngram = seq_to_bin(ngram)
            seen.setdefault(ngram, 0)
            seen[ngram] += 1
            n += 1
    print "Sorting..."
    seen = sort_dictionary_by_value(seen, reverse=True)
    print "Formating..."
    seen = ["%s\t%s\n" % (bin_to_seq(k), v) for v, k in seen]
    print "Saving..."
    with open(output_ngram, "w") as fh:
        for line in seen:
            fh.write(line)
    print "Done."
    print "ngrams, skipped: ", n, skipped
    return n, skipped
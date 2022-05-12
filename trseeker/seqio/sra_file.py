#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 05.06.2011
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com 
'''
SRA files.
TODO: implement this.
TODO: separate model from reader
'''

from trseeker.tools.sequence_tools import get_revcomp, get_gc
from trseeker.tools.sequence_tools import clear_sequence
import re
try:
    from Bio import pairwise2
except:
    print("Not Bio installed")
import random
from PyExp import WiseOpener

class FastqObj(object):

    def __init__(self, head, seq, strain, qual, phred33=False):
        self.head = head.strip()
        self.seq = seq.strip()
        self.qual = qual.strip()
        self.strain = strain.strip()
        self.id = head.strip().replace("#", " ").split()[0].replace("@","")
        # trimmeing input
        self.trimmed_seq = self.seq.strip()
        self.trimmed_qual = qual.strip()
        if phred33:
            self.qv = [ord(x)-33 for x in qual.strip()]
        # trimming results
        self.adaptor_positions = []
        self.adapter_contamination = None
        self.qual_junk = None
        self.status = True
        self.parts = []
        self.id = self.read_id_hiseq15
        


    @property
    def fastq(self):
        return "%s\n%s\n%s\n%s\n" % (
                    self.head,
                    self.seq.upper(),
                    self.strain,
                    self.qual,
            )


    def change_strain(self):
        if self.strain == "+":
            self.strain = "-"
        elif self.strain == "-":
            self.strain = "+"

    def split_read(self, position, step):
        ''' Return two new read_objs splitted by position.
        '''
        head1 = self.head+":s1"
        head2 = self.head+ ":s2"
        seq1 = self.seq[:position]
        seq2 = self.seq[position-step:]
        qual1 = self.qual[:position]
        qual2 = self.qual[position-step:]
        read1 = FastqObj(head1, seq1, self.strain, qual1)
        read2 = FastqObj(head2, seq2, self.strain, qual2)
        return read1, read2


    def check_and_fix_read_q_legnths(self):

        l1 = len(self.seq)
        l2 = len(self.qual)
        if l1 != l2:
            l = min(l1, l2)
            self.seq = self.seq[:l]
            self.qual = self.qual[:l]
            return 1
        return 0


    def make_reversed(self):
        '''
        '''
        self.seq = self.seq[::-1]
        self.qual = self.qual[::-1]

    def fastq_with_error(self, error):
        return "@%s__%s\n%s\n%s\n%s\n" % (
                    error,
                    self.head[1:],
                    self.seq,
                    self.strain,
                    self.qual,
            )

    def fastq_with_cov(self, cov):
        cov = ",".join(map(str, cov))
        return "%s%s\n%s\n%s\n" % (
                    self.head,
                    self.trimmed_seq,
                    cov,
                    self.trimmed_qual,
            )

    @property
    def trimmed_fastq(self):
        return "%s%s\n%s%s\n" % (
                    self.head,
                    self.trimmed_seq,
                    self.strain,
                    self.trimmed_qual,
            )

    @property
    def sequence(self):
        return self.seq.strip()

    @property
    def fasta(self):
        return ">%s\n%s\n" % (
                    self.head.strip(),
                    self.seq.strip(),
                    )

    @property
    def gc(self):
        return get_gc(self.seq)

    @property
    def length(self):
        return len(self.seq)

    @property
    def read_id(self):
        self.id = int(self.head.split()[0].split(".")[-1])
        return self.id

    @property
    def read_id_and_pair(self):
        try:
            self.id = int(self.head.split()[0].split(".")[-1])
        except:
            self.id = self.head.split()[0]
        self.pair = int(self.head.split()[1].split(":")[0])
        return self.id, self.pair

    @property
    def read_id_hiseq15(self):
        self.id = self.head.replace("#", " ").split()[0]
        return self.id

    def trim(self):
        '''
        '''
        raise NotImplemented
                    
    def trim_by_quality_cutoff(self, cutoff=30):
        '''
        '''
        for i, q in enumerate(self.qv):
            if q < cutoff:
                self.trimmed_seq = self.trimmed_seq[:i]
                self.trimmed_qual = self.trimmed_qual[:i]
                self.qual_junk = self.trimmed_seq[i:]
                break

    def trim_exact_adaptor(self, adaptor, verbose=False):
        '''
        '''
        raise NotImplemented
       
    def trim_inexact_adaptor(self, adaptor, verbose=False, num_errors=2):
        '''
        '''
        raise NotImplemented

    def iter_kmers(self, k=23):
        '''
        '''
        for i in xrange(len(self.seq)-k+1):
            yield i, self.seq[i:i+k]

class PERun(object):
    pass

def fastq_reader(fastq_file, phred33=False, head_pattern=None):
    with WiseOpener(fastq_file) as fh:
        while True:
            try:
                head = fh.readline()
                if head == '':
                    break
                if head_pattern:
                    if not re.match(head_pattern, head):
                        continue
                seq = fh.readline()
                strain = fh.readline()
                qual = fh.readline()
                fastq_obj = FastqObj(head, seq, strain, qual, phred33=phred33)
                yield fastq_obj
            except:
                break

def fastq_pe_reader(fastq_file1, fastq_file2, phred33=False, head_pattern=None):
    with WiseOpener(fastq_file1) as fh1:
        with WiseOpener(fastq_file2) as fh2:
            while True:
                try:
                    head1 = fh1.readline()
                    head2 = fh2.readline()
                    if head1 == '':
                        break
                    if head_pattern:
                        if not re.match(head_pattern, head1):
                            continue
                    seq1 = fh1.readline()
                    strain1 = fh1.readline()
                    qual1 = fh1.readline()
                    fastq_obj1 = FastqObj(head1, seq1, strain1, qual1, phred33=phred33)

                    seq2 = fh2.readline()
                    strain2 = fh2.readline()
                    qual2 = fh2.readline()
                    fastq_obj2 = FastqObj(head2, seq2, strain2, qual2, phred33=phred33)

                    yield fastq_obj1, fastq_obj2
                except:
                    break

def fastq_iter_seqs_pe(fastq_file1, fastq_file2):
    with WiseOpener(fastq_file1) as fh1:
        with WiseOpener(fastq_file2) as fh2:
            while True:
                try:
                    fh1.readline()
                    seq1 = fh1.readline()
                    fh1.readline()
                    fh1.readline()
                    
                    fh2.readline()
                    seq2 = fh2.readline()
                    fh2.readline()
                    fh2.readline()
                    
                    yield seq1, seq2
                except:
                    break

def read_qual_reader(read_file, qual_file, phred33=False):
    """ Read fastq_obj from two files:
        with reads
        with corresponding quals
    """
    with open(read_file) as fh:
        with open(qual_file) as qfh:
            rid = 0
            while True:
                try:
                    seq = fh.readline()
                    qual = qfh.readline()
                    fastq_obj = FastqObj("@%s" % rid, seq, "+", qual, phred33=phred33)
                    yield fastq_obj
                    rid += 1
                except Exception as e:
                    print(e)
                    break

def read_simple_qual_reader(read_file, qual_file, phred33=False):
    """ Read fastq_obj from two files:
        with reads
        with corresponding quals
    """
    with open(read_file) as fh:
        with open(qual_file) as qfh:
            rid = 0
            while True:
                try:
                    yield rid, fh.readline().strip(), qfh.readline().strip()
                    rid += 1
                except Exception as e:
                    print(e)
                    break

def error_aware_fastq_reader(fastq_file, phred33=False, header=None):
    errors = 0
    ok = 0
    data = []
    with open(fastq_file) as fh:
        while True:
            # check head
            head = fh.readline()
            if not head:
                return
            if header and not head.startswith(header):
                print("Error with head:", head.strip())
                errors += 1
                continue
            seq =  fh.readline()
            if not seq[0].upper() in ["A", "C", "T", "G", "N"]:
                print()
                print("Error with seq:", seq.strip())
                errors += 1
                continue
            strain  = fh.readline()
            if header and strain.startswith(header):
                    print()
                    print("Error with strain:", strain.strip())
                    errors += 1
                    continue
            qual = fh.readline()
            if header and qual.startswith(header):
                    print()
                    print("Error with qual:", qual.strip())
                    errors += 1
                    continue
            if not len(qual) == len(seq):
                    print()
                    print("Error with qual length:", qual.strip())
                    errors += 1
                    continue
            ok += 1
            fastq_obj = FastqObj(head, seq, strain, qual, phred33=phred33)
            # print errors, ok
            yield fastq_obj


def generate_fastq(seq, i, pair):
        head = "@ART_%s %s" % (i, pair)
        qual = "".join([random.choice(["A","B","C","D","E","F","G","I"]) for x in range(len(seq))])
        return "%s\n%s\n+\n%s\n" % (
                    head,
                    seq.upper(),
                    qual,
            )


def iter_sorted_pe_data(fastq1_file, fastq2_file):
    ''' Iterate over sorted PE fastq files.
    '''
    reader1 = error_aware_fastq_reader(fastq1_file)
    reader2 = error_aware_fastq_reader(fastq2_file)
    while True:
        read1 = None
        read2 = None
        if reader1 is None and reader2 is None:
            return
        if reader1:
            try:
                read1 = next(reader1)
            except StopIteration:
                read1 = None
                reader1 = None
        if reader2:
            try:
                read2 = next(reader2)
            except StopIteration:
                read2 = None
                reader2 = None
        if not read1:
            yield None, read2
            continue
        if not read2:
            yield read1, None
            continue
        if read1 and read2 and read1.id == read2.id:
            yield read1, read2
            continue
        if read1 and read2 and read1.id < read2.id:
            while True:
                yield read1, None
                if reader1:
                    try:
                        read1 = next(reader1)
                    except StopIteration:
                        read1 = None
                        reader1 = None
                        break
                else:
                    break   
                if read1.id == read2.id:
                    yield read1, read2
                    break
        if read1 and read2 and read1.id > read2.id:
            while True:
                yield None, read2
                if reader2:
                    try:
                        read2 = next(reader2)
                    except StopIteration:
                        read2 = None
                        reader2 = None
                        break
                else:
                    break
                if read1.id == read2.id:
                    yield read1, read2
                    break
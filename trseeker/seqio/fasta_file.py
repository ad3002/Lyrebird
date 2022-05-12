#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 05.06.2011
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com 
"""
Classes:

- FastaFileIO(AbstractBlockFileIO)
  
"""
import os
import pickle
from trseeker.seqio.block_file import AbstractBlockFileIO
from trseeker.models.sequence_model import SequenceModel


class FastaFileIO(AbstractBlockFileIO):
    """  Working with multi fasta files, 
    where each block starts with '>' token.
    
    Overrided public methods:
    
    - __init__(self)
    
    Inherited public properties:
    
    - data  - iterable data, each item is tuple (head, body)
    - N     - a number of items in data
    
    Inherited public methods:
    
    - read_from_file(self, input_file)
    - read_online(self, input_file) ~> item
    - get_block_sequence(self, head_start, next_head, fh)
    - get_blocks(self, token, fh)
    - gen_block_sequences(self, token, fh)
    - read_from_db(self, db_cursor)
    - write_to_file(self, output_file)
    - write_to_db(self, db_cursor)
    - read_as_iter(self, source)
    - iterate(self) ~> item of data 
    - do(self, cf, args) -> result
    - process(self, cf, args)
    - clear(self)
    - do_with_iter(self, cf, args) -> [result,]
    - process_with_iter(self, cf, args)
        
    """

    def __init__(self):
        """ Overrided. Hardcoded start token."""
        token = ">"
        super(FastaFileIO, self).__init__(token)


def sc_iter_reads(file_name):
    """
    """
    with open(file_name) as fh:
        for sequence in fh:
            seq_obj = SequenceModel()
            seq_obj.set_ncbi_sequence("unknown", sequence.strip())
            yield seq_obj            


def sc_iter_fasta(file_name, lower=False, protein=False, inmem=True, skip_clean=False):
    """ Iter over fasta file."""
    header = None
    seq = []
    with open(file_name) as fh:
        data = fh
        if inmem:
            data = fh.readlines()
        for line in data:
            if line.startswith(">"):
                if seq:
                    sequence = "".join(seq)
                    seq_obj = SequenceModel(lower=lower, protein=protein)
                    seq_obj.set_ncbi_sequence(header, sequence, skip_clean=skip_clean)            
                    yield seq_obj
                header = line
                seq = []
                continue
            seq.append(line)
        if seq or header:
            sequence = "".join(seq)
            seq_obj = SequenceModel(lower=lower, protein=protein)
            seq_obj.set_ncbi_sequence(header, sequence)            
            yield seq_obj

def sc_iter_fasta_brute(file_name, inmem=False, lower=False):
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
                if seq or header:
                    sequence = "".join(seq)
                    if lower:
                        sequence = sequence.lower()
                    yield header, sequence
                header = line.strip()
                seq = []
                continue
            seq.append(line.strip())
        if seq or header:
            sequence = "".join(seq)
            if lower:
                sequence = sequence.lower()
            yield header, sequence


def sc_iter_fasta_for_refseq(file_name):
    """ Iter over fasta file."""

    header = None
    seq = []
    with open(file_name) as fh:
        data = fh.readlines()
        for line in data:
            if line.startswith(">"):
                if seq:
                    sequence = "".join(seq)
                    yield (header, sequence)
                header = line.strip().split()[0][1:]
                seq = []
                continue
            seq.append(line.upper().strip())
        if seq or header:
            sequence = "".join(seq)
            yield (header, sequence)

def sc_iter_fasta_slow(file_name, lower=True, protein=False):
    """ Iter over fasta file."""

    reader = FastaFileIO()
    for (head, sequence, start, next) in reader.read_online(file_name):
        seq_obj = SequenceModel(lower=lower, protein=protein)
        seq_obj.set_ncbi_sequence(head, sequence)
        yield seq_obj

def fasta_reader(file_name):
    """ Synonym for  sc_iter_fasta.
    """
    return sc_iter_fasta(file_name)

def sc_iter_fasta_simple(file_name, lower=True):
    """ Iter over fasta file."""

    reader = FastaFileIO()
    for (head, sequence, start, next) in reader.read_online(file_name):
        seq_obj = SequenceModel(lower=lower)
        seq_obj.set_ncbi_sequence(head, sequence)
        yield seq_obj.seq_gi, seq_obj.sequence

def save_all_seq_with_exact_substring(fasta_file, substring, output_file):
    """ Save all fasta sequences with exact match with given substring in output_file.
    """
    with open(output_file, "w") as fh:
        for seq_obj in sc_iter_fasta(fasta_file):
            if substring in seq_obj.seq_sequence:
                fh.write(seq_obj.fasta)


def sort_fasta_file_by_length(file_name):
    """ Sort by length and save fasta file
    """
    objs = []
    for i, seq_obj in enumerate(sc_iter_fasta(fasta_file)):
        objs.append(seq_obj)
    objs.sort(key=lambda x: x.length)
    with open(fasta_file, "w") as fh:
        for obj in objs:
            fh.write(obj.fasta)


def sc_get_fasta_lengths(fasta_file):
    """ Compute or read from precomputed pickle file
    dict header to seq length for input fasta_file and save it with pickle.
    """
    length_file = fasta_file + ".length.p"
    if os.path.isfile(length_file):
        try:
            with open(length_file) as fh:
              head2lengths = pickle.load(fh)
            return head2lengths
        except Exception as e:
            print(e)
            print("Can't read file %s" % length_file)
    head2lengths = {}
    for seq_obj in sc_iter_fasta(fasta_file):
        head2lengths[seq_obj.header] = seq_obj.length
    with open(length_file, "w") as fh:
        pickle.dump(head2lengths, fh)
    return head2lengths


def sc_chr2seq_from_fasta(file_name):
    """ Iter over fasta file and return chr name to sequence."""
    chr2seq = {}
    header = None
    seq = []
    with open(file_name) as fh:
        data = fh.readlines()
        for line in data:
            if line.startswith(">"):
                if seq:
                    chr2seq[header] = "".join(seq)
                header = line.strip().split()[0][1:]
                seq = []
                continue
            seq.append(line.upper().strip())
        if seq or header:
            chr2seq[header] = "".join(seq)
    return chr2seq


def sc_add_prefix_to_fasta_header(fasta_file, output_file, prefix, clean=False):
    '''Add custom prefix to each header in multifasta.
    '''
    with open(output_name, "w") as fw:
        for header, seq in sc_iter_fasta_brute(fasta_file):
            if clean:
                name = header[1:].split()[0]
            header = "%s_%s" % (prefix, name)
            fw.write(">%s\n%s\n" % (header, seq)) 

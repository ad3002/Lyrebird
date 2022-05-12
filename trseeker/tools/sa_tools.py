#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 11.04.2010
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com 
'''
Functions related to suffix array construction.
TODO: checkit
'''
import os
import pickle

from trseeker.tools.sequence_tools import clear_sequence
from trseeker.seqio.fasta_file import sc_iter_fasta_simple

def fasta_to_sa_input(fasta_file, sa_input_file, index_file_name=None, file_name=None, start_id=None, increment=False):
    ''' Fasta to SA input. Text file of sequences delimeted by $ symbol.'''
    if not index_file_name:
        index_file_name = fasta_file + ".sa.index"

    if not increment:
        if os.path.isfile(sa_input_file):
            os.unlink(sa_input_file)
        if os.path.isfile(index_file_name):
            os.unlink(index_file_name)

    with open(index_file_name, "a") as index_fw:
        with open(sa_input_file, "a") as fa_fw:
            if start_id:
                i = start_id
            else:
                i = 0
            for seq_obj in sc_iter_fasta_simple(fasta_file):

                # write sa input
                if i > 0:
                    fa_fw.write("$")
                i += 1
                sequence = seq_obj.sequence
                head = seq_obj.seq_head
                fa_fw.write(sequence)

                # write index
                if file_name:
                    head_info = "%s:%s" % (file_name, head.strip())
                else:
                    head_info = head.strip()
                head_info = head_info.replace(">", "")
                index_fw.write("%s\t%s\t%s\n" % (i, head_info.replace("\t", " "), len(sequence)))
    return i

def sa_input_to_fasta(input_file, output_file):
    ''' SA input file to fasta file.'''
    if os.path.isfile(output_file):
        os.remove(output_file)
    with open(input_file) as fh:
        data = fh.read()
    data = data.strip().split("$")
    with open(output_file, "a") as fw:
        for i, item in enumerate(data):
            item = item.strip()
            if not item:
                continue
            fw.write(">%s\n%s\n" % (i, item))

def pickle_dictionary_for_docid_trid(sa_doc_index, doc_to_trf_file, trf_to_doc_file):
    ''' Function precompiles doc2id and id2doc pickled dictionary.
    '''
    doc_to_trf = {}
    trf_to_doc = {}
    with open(sa_doc_index) as fh:
        for line in fh:
            data = line.strip().split("\t")
            docid = int(data[0]) - 1
            id = data[1].split(">")[-1]
            try:
                trfid = int(id)
            except:
                trfid = docid + 1
            doc_to_trf[docid] = trfid
            trf_to_doc[trfid] = docid

    with open(doc_to_trf_file, "w") as dffh:
        with open(trf_to_doc_file, "w") as tdfh:
            pickle.dump(doc_to_trf, dffh)
            pickle.dump(trf_to_doc, tdfh)

def filter_sa_dataset(sa_file, output_file, min_tf, min_df):
    ''' Function filters and writes sa data with given min tf and df.
    '''
    dataset = []
    with open(sa_file) as fh:
        for line in fh:
            if line.startswith("type"):
                continue
            data = line.strip().split("\t")
            tf = int(data[4])
            td = int(data[5])
            if tf < min_tf and td < min_df:
                continue
            dataset.append((td, line))
        dataset.sort(reverse=True)
        with open(output_file, "a") as fw:
            for item in dataset:
                fw.write(item[1])

def iterate_sa_corpus(corpus):
    ''' Yields corpus texts. Corpus is $ delimited text file.
    '''
    left_pos = None
    right_pos = None
    N = len(corpus)
    while True:
        if left_pos is None:
            left_pos = corpus.find("$", 0)
            if left_pos < 0:
                yield corpus
                break
            yield corpus[0:left_pos]
        right_pos = corpus.find("$", left_pos + 1)
        if left_pos == N - 1:
            break
        if right_pos < 0:
            yield corpus[left_pos + 1:]
            break
        yield corpus[left_pos + 1:right_pos]
        left_pos = right_pos
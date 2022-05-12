#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 21.05.2014
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com

from PyBioSnippets.Lyrebird.tools.checkers import *
from PyBioSnippets.Lyrebird.tools.adders import *
from PyBioSnippets.Lyrebird.tools.classifiers import *
from PyBioSnippets.Lyrebird.tools.uploaders import *
from PyBioSnippets.Lyrebird.tools.connectors import *
import math
from collections import defaultdict

def learn_add_trs(i, sequence, kmer2tf, kmer, k, taxon, prefix, lyrebird, tr_family_to_letter, annotation_cache, fast=False):
    '''
    '''
    print("UNKNOWN  %s %s %s ?" % (i, kmer2tf[kmer], kmer))
    meta = None
    length = len(sequence)
    sequence, inner_name, inner_family = submit_trs_dialog(sequence, length, k, prefix, lyrebird, kmer2tf, tr_family_to_letter, annotation_cache, fast=fast)
    if inner_family:
        add_lyrebird_seq(lyrebird, sequence, taxon, inner_name, inner_family)
    else:
        choice = raw_input("Submit fragment (start:end)?") or None
        if choice:
            start, end = map(int, choice.strip().split(":"))
            sequence = sequence[start-1:end-1]
            submit_fragment_dialog(sequence, len(sequence), k, prefix, lyrebird)
        else:
            print("Repeat was skipped")
    return inner_name, inner_family, []


def learn_add_te(i, sequence, kmer2tf, kmer, k, taxon, prefix, lyrebird, annotation_cache, default_name=None, default_family=None, fast_te_adding=False):
    '''
    '''
    print("UNKNOWN  %s %s %s ?" % (i, kmer2tf[kmer], kmer))
    meta = None
    length = len(sequence)
    sequence, inner_name, inner_family = submit_te_dialog(sequence, length, k, prefix, lyrebird, kmer2tf, annotation_cache, default_name=default_name, default_family=default_family, fast_te_adding=fast_te_adding)
    if inner_family:
        add_lyrebird_seq(lyrebird, sequence, taxon, inner_name, inner_family)
    else:
        raise NotImplemented
        choice = raw_input("Submit fragment (start:end)?") or None
        if choice:
            start, end = map(int, choice.strip().split(":"))
            sequence = sequence[start-1:end-1]
            submit_fragment_dialog(sequence, len(sequence), k, prefix, lyrebird)
        else:
            print("Repeat was skipped")
    return inner_name, inner_family, []


def ask_priority(default=100):
    ''' Ask priority.
    '''
    decision = raw_input("Say priority (number, default 100)?") or default
    try:
        decision = int(decision)
        return decision
    except:
        print("Priority should be a number")
        return ask_priority()


def submit_trs_dialog(sequence, length, k, prefix, lyrebird, kmer2tf, tr_family_to_letter, annotation_cache, skip=False, fast=False):
    ''' TRs submitting pipeline. 
    '''
    print("TRs classification:")
    monomer, inner_name, inner_family = ask_trs_name(sequence, length, prefix, tr_family_to_letter, fast=fast)
    print("TRs name: %s, family name: %s" % (inner_name, inner_family))
    if skip:
        return sequence, inner_name, inner_family
    # extend sequences
    sequence = monomer_to_sequence(monomer, k)
    annotation, summary = print_annotation(length, sequence, k, lyrebird, kmer2tf, annotation_cache)
    
    decision = 'y'
    if not len(summary) == 1:
        decision = raw_input("Submit (y/n, default y)?") or 'y'
    if decision == 'y':
        priority = 100
        # priority = ask_priority()
        add_kmers(lyrebird, sequence, inner_name, inner_family, priority, k=k, prev_known=None, rewrite_name=False)
        return sequence, inner_name, inner_family
    return sequence,None,None


def submit_te_dialog(sequence, length, k, prefix, lyrebird, kmer2tf, annotation_cache, skip=False, default_name=None, default_family=None, fast_te_adding=False):
    ''' TRs submitting pipeline. 
    '''
    print("TE classification:")
    annotation, summary = print_annotation(length, sequence, k, lyrebird, kmer2tf, annotation_cache)
    if fast_te_adding:
        inner_name, inner_family = default_name, default_family
    else:
        inner_name, inner_family = ask_name(summary=summary, default_name=default_name, default_family=default_family)
    print("TE name: %s, family name: %s" % (inner_name, inner_family))
    if skip:
        return sequence, inner_name, inner_family
    # extend sequences
    decision = 'y'
    if not len(summary) == 1:
        decision = raw_input("Submit (y/n/f, default y, f - fragment)?") or 'y'
    if decision == 'y':
        priority = 100
        # priority = ask_priority()
        add_kmers(lyrebird, sequence, inner_name, inner_family, priority, k=k, prev_known=None, rewrite_name=False)
        return sequence, inner_name, inner_family
    elif decision == 'f':
        raise NotImplemented
    return sequence,None,None


def submit_by_known(sequence, k, known_obj, lyrebird):
    '''
    '''
    print("Found sequence:", sequence)
    inner_name = known_obj["inner_name"]
    inner_family = known_obj["inner_family"]
    
    q = raw_input("Submit by as %s:%s repeat?" % (inner_name, inner_family)) or "y"
    if q != "y":
        print("Skipped")
        return sequence,None,None
    add_repeat(sequence, inner_name, inner_family, k=k)
    return sequence, inner_name, inner_family


def print_annotation(length, sequence, k, lyrebird, kmer2tf, annotation_cache):
    '''
    '''
    print("Annotation: (%s)" % (sequence))
    annotation = make_full_check(sequence, lyrebird, k, kmer2tf, annotation_cache)
    print("pos\tmonpos\tNuc\tKmer\tLB\tRB23\rRF\rRB13")
    pos = 0
    summary = defaultdict(int)

    min_tf, max_tf = 10000000000000, 0

    for i in range(len(sequence)-k+1):
        pos += 1
        if pos > length:
            pos = 1
        kmer = sequence[i:i+k]
        tf = kmer2tf[kmer]
        min_tf = min(min_tf, tf)
        max_tf = max(max_tf, tf)
        if annotation["lyrebird"][kmer]["inner_name"] == 'Unknown' == annotation["repbase"][kmer] == annotation["rfam"][kmer]:
            print(i+i, pos, kmer, tf)
            summary[("Unknown", "Unknown")] += 1
            continue
        print(i+1, pos, sequence[i], annotation["kmers"][kmer], annotation["lyrebird"][kmer]["inner_name"], end=" ")
        summary[(annotation["lyrebird"][kmer]["inner_name"], annotation["lyrebird"][kmer]["inner_family"])] += 1
        if k == 23:
            print(annotation["repbase"][kmer], annotation["rfam"][kmer], annotation["repbase_k13"][kmer])
        else:
            print("-","-","-")

    for x in summary:
        print(x, summary[x])

    print("Min tf: %s | max tf: %s" % (min_tf, max_tf))
    print(sequence)

    return annotation, summary


def submit_fragment_dialog(sequence, length, k, prefix, lyrebird):
    '''
    '''
    print("Found sequence:", sequence)
    inner_name, inner_family = ask_name()
    return add_repeat(sequence, inner_name, inner_family, k=k)
     

def manual_submission_trs(sequence, name, family, taxon, prefix):
    '''
    1) submit sequence to LyrebirdSeqs
    2) submit kmers
    3) enrich kmers over ed<3
    '''
    add_lyrebird_seq(lyrebird, sequence, taxon, name, family, priority=300, repbase_name=None, proofs=None)
    upload_sequences([sequence*2], family, name, k=23, kid=-1)


from trseeker.tools.sequence_tools import get_revcomp
from trseeker.seqio.tab_file import sc_iter_tab_file
from trseeker.models.trf_model import TRModel

def manual_submission_from_trf(trf_file, kmer, name, taxon, family, k=23):
    '''
    1) grep trf_file with kmer
    2) get arrays
    3) compile kmers, get 80% fraction
    4) submit sequences
    5) submit kmers
    '''
    rkmer = get_revcomp(kmer)
    sequences = []
    for i, trf_obj in enumerate(sc_iter_tab_file(trf_file, TRModel)):
        if kmer.lower() in trf_obj.trf_array or rkmer.lower() in trf_obj.trf_array:
            print("Found match:", trf_obj)
            sequence = get_revcomp(trf_obj.trf_array.upper())
            if rkmer in trf_obj.trf_array:
                sequences.append(sequence)
            else:
                sequences.append(sequence)
            proofs = {
                "trf": trf_obj.get_as_dict(),
            }
            add_lyrebird_seq(lyrebird, trf_obj.trf_consensus, taxon, name, family, priority=300, repbase_name=None, proofs=proofs)
    upload_sequences(sequences, family, name, k=23, kid=-1)


def manual_submission_from_te(sequence, name, taxon, family, k=23):
    '''
    '''
    
    proofs = {
        "manual": True,
    }
    add_lyrebird_seq(lyrebird, sequence, taxon, name, family, priority=300, repbase_name=None, proofs=proofs)
    upload_sequences([sequence], family, name, k=23, kid=-1)


def monomer_to_sequence(monomer, k):
    ''' Convert monomer to sequence formed by duplicated 
    monomer for tandem repeats with cases for short monomers.
    '''
    sequence = monomer
    length = len(monomer)
    if length < k:
        sequence = monomer * (k//length + 1)
    sequence *= 2
    return sequence


if __name__ == '__main__':
    # sequence = raw_input("Insert sequence: ").upper().strip()
    # family = raw_input("Family (e.g. SatDNA): ").strip()
    # name = raw_input("Name: ").strip()
    # taxon = raw_input("Taxon: ").strip()
    # prefix = raw_input("Prefix: ").strip()

    # manual_submission_trs(sequence, name, family, taxon, prefix)
    tdata = raw_input("TR|TE: ").strip() or "TR"
    if tdata == "TR":
        trf_file = raw_input("TRF file: ").strip() or "trf_all.trf"
        kmer = raw_input("Kmer: ").upper().strip()
        family = raw_input("Family (e.g. SatDNA): ").strip() or "SatDNA"
        name = raw_input("Name: ").strip()
        taxon = raw_input("Taxon: ").strip() or "Anopheles_stephensi"
        manual_submission_from_trf(trf_file, kmer, name, taxon, family, k=23)
    else:
        sequence = raw_input("Sequence: ").upper().strip()
        family = raw_input("Family (e.g. SINE): ").strip() or "SINE"
        name = raw_input("Name: ").strip()
        taxon = raw_input("Taxon: ").strip() or "Anopheles_stephensi"
        manual_submission_from_te(sequence, name, taxon, family, k=23)



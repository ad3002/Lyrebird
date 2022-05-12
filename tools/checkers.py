#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 21.05.2014
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com

from collections import defaultdict
# from trseeker.tools.jellyfish_tools import query_kmers
from trseeker.tools.sequence_tools import get_revcomp

def get_assembly_proofs(sequence, settings):
    raise NotImplemented


def get_short_se_proofs(sequence, settings):
    raise NotImplemented


def get_short_pe_proofs(sequence, settings):
    raise NotImplemented


def get_pacbio_proofs(sequence, settings):
    raise NotImplemented


def check_with_kmers(sequence, kmer2tf, k, verbose=False, annotation=None, jellyfish_dat=None):
    """ Add kmer frequencies annotation in give jellyfish database.
    :param jellyfish_db: jellyfish db file
    :param verbose: print annotation
    :param annotation: default dictionary with external annotation
    :param jellyfish_dat: dump of jellyfish database
    """
    kmers = set()
    result = defaultdict(int)
    for i in range(len(sequence)-k+1):
        kmer = sequence[i:i+k].upper()
        result[kmer] = kmer2tf[kmer]
        if verbose and annotation:
            for i in range(len(sequence)-k+1):
                kmer = sequence[i:i+k]
                print("\t\t%s\t%s\t%s" % (kmer, result[kmer], annotation[kmer]))
    return result


def check_with_rfam(sequence, lyrebird):
    """ Get information about kmer found in rfam database.
    """
    k = 23
    kmers = set()
    for i in range(len(sequence)-k+1):
        kmer = sequence[i:i+k]
        kmers.add(kmer)
    kmers = list(kmers)
    results = {}
    for kmer in kmers:
        result = lyrebird.query_rfam(kmer)
        if result:
            name = list(set([x.split(":")[0] for x in result["inner_name"]]))
            name.sort()
            name = ":".join(name)
            results[kmer] = name
        else:
            results[kmer] = "Unknown"
    return results


def check_with_repbase(sequence, lyrebird):
    """ Get information about kmer found in Repbase database with k=23.
    """
    k = 23
    kmers = set()
    for i in range(len(sequence)-k+1):
        kmer = sequence[i:i+k]
        kmers.add(kmer)
    kmers = list(kmers)
    results = {}
    for kmer in kmers:
        result = lyrebird.query_repbase(kmer)
        if result:
            name = list(set([x.split("&")[2] for x in result["inner_name"]]))
            name.sort()
            name = ":".join(name)
            results[kmer] = name
        else:
            results[kmer] = "Unknown"
    return results


def check_with_repbase_k13(sequence, lyrebird):
    """ Get information about kmer found in Repbase database with k=13.
    """
    k = 23
    kmers = set()
    for i in range(len(sequence)-k+1):
        kmer = sequence[i:i+k]
        kmers.add(kmer)
    kmers = list(kmers)
    results = {}
    for kmer in kmers:
        result = lyrebird.query_repbase_k13(kmer)
        if result:
            name = list(set([x.split("&")[2] for x in result["inner_name"]]))
            name.sort()
            name = ":".join(name)
            results[kmer] = name
        else:
            results[kmer] = "Unknown"
    return results


def check_with_lyrebird(sequence, lyrebird, k, annotation_cache):
    """ Get information about kmer found in Lyrebird database with k=13.
    """
    kmers = set()
    for i in range(len(sequence)-k+1):
        kmer = sequence[i:i+k]
        kmers.add(kmer)
    kmers = list(kmers)
    results = {}
    for kmer in kmers:
        if kmer in annotation_cache:
            result = annotation_cache[kmer]
        else:
            result = lyrebird.query_lyrebird(kmer)
        if result:
            results[kmer] = result
        else:
            results[kmer] = {"inner_name": "Unknown",
                             "inner_family": "Unknown",
            }
    return results


def make_full_check(sequence, lyrebird, k, kmer2tf, annotation_cache):
    """ Get all available annotations.
    """
    # Get unique kmers

    repbase = defaultdict(str)
    repbase_k13 = defaultdict(str)
    rfam = defaultdict(str)

    # repbase = check_with_repbase(sequence, lyrebird)
    # repbase_k13 = check_with_repbase_k13(sequence, lyrebird)
    # rfam = check_with_rfam(sequence, lyrebird)

    lb = check_with_lyrebird(sequence, lyrebird, k, annotation_cache)
    km = check_with_kmers(sequence, kmer2tf, k)

    annotation = {
        "repbase": repbase,
        "repbase_k13": repbase_k13,
        "rfam": rfam,
        "lyrebird": lb,
        "assembly": defaultdict(str),
        "short_se": defaultdict(str),
        "short_mp": defaultdict(str),
        "pacbio": defaultdict(str),
        "kmers": km,
    }
    return annotation

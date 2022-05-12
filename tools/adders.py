#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 21.05.2014
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com

from trseeker.tools.sequence_tools import get_revcomp, get_dust_score


def add_lyrebird_seq(lyrebird, sequence, taxon, inner_name, inner_family, priority=100, repbase_name=None, proofs=None):
    """ Simply add a sequence to lyrebird db.

    :param lyrebird: lyrebird object
    :param sequence: sequence
    :param taxon: taxon
    :param inner_name: inner_name
    :param inner_family: inner_family
    :param priority: priority, default 100
    :param repbase_name: repbase name, default None
    :param proofs: proofs dictionary
    """
    full_name = ":".join((sequence, inner_name, inner_family, taxon))

    if not proofs:
        proofs = {}

    data = {'sequence': sequence,
            'full_name': full_name,
            'priority': priority,
            'family_name': inner_family,
            'file_name': "lyrebird",
            'inner_name': inner_name,
            'repbase_name': repbase_name,
            'taxons': taxon,
            'name_taxon': "%s:%s" % (inner_name, taxon),
            'meta': proofs,
            }

    obj = lyrebird.add_lyrebird_seq(data)
    return obj


def add_kmers(lyrebird, sequence, inner_name, inner_family, priority, k=23, prev_known=None, rewrite_name=False, delta=0):
    """ Add kmers from sequence to lyrebird database.

    :param lyrebird: lyrebird object
    :param sequence: sequence
    :param inner_name: inner_name
    :param inner_family: inner_family
    :param k: k for kmer, default 23
    :param priority: priority
    :param rewrite_name: rewrite_name, default False
    :param prev_known: list of previously added kmers for speed up
    """
    kmers = set()
    sequence = sequence.upper()
    for i in range(len(sequence)-k+1):
        kmer = sequence[i:i+k]
        rev_kmer = get_revcomp(kmer)

        if prev_known and (kmer in prev_known or rev_kmer in prev_known):
            continue

        if kmer > rev_kmer:
            kmer = rev_kmer

        data = {'kid': -1,
                'pos': i,
                'kmer': kmer,
                'inner_name': inner_name,
                'inner_family': inner_family,
                'ref_kmer': kmer,
                'priority': priority,
                'ed': 0,
                'hd': 0,
                'dust': get_dust_score(kmer),
                'delta': delta,
                }
        lyrebird.add_update(data, rewrite_name=rewrite_name, skip_if_known=False, db=None)

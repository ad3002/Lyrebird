#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 21.05.2014
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com
"""
Databases:

data = {
        'kid': kid,
        'kmer': kmer,
        'inner_name': names,
        'inner_family': family_name,
        'ref_kmer': kmer,
        'ed': 0,
        'hd': 0,
        'dust': get_dust_score(kmer),
    }

Inner names:
    priority_family_name

    pr010_Microsatellite_(AC)n

Database with sequences:

    full_name
    sequence
    priority
    family_name
    inner_name
    repbase_name
    taxons

priority list

010 perfect Microsatellites
015 Microsatellites
020 SSR


"""

import os
import simplejson
from trseeker.tools.sequence_tools import get_shifts_variants 
from trseeker.tools.sequence_tools import get_revcomp
from trseeker.tools.sequence_tools import get_dust_score
from PyBioSnippets.Lyrebird.tools.connectors import *
from trseeker.seqio.fasta_file import sc_iter_fasta
from trseeker.tools.ngrams_tools import process_list_to_kmer_index
from PyExp import runner


settings = {
    "rDNA_fasta": "/home/akomissarov/Dropbox/PyBioSnippets/Lyrebird/fasta_libraries/rDNA.fa",
    "mtDNA_fasta": "/home/akomissarov/Dropbox/PyBioSnippets/Lyrebird/fasta_libraries/mtDNA.fa",
    "rfam_fasta": "/home/akomissarov/Dropbox/PyBioSnippets/Lyrebird/fasta_libraries/Rfam.fasta",
    "repbase_fasta": "/home/akomissarov/repbase/repbase.fasta",
    "repbase_index": "/home/akomissarov/repbase/repbase.index",
    "sinebase_fasta": "/home/akomissarov/Dropbox/PyBioSnippets/Lyrebird/fasta_libraries/SINEbase_21_04_14.fa",
    "tRNA_fasta": "/home/akomissarov/Dropbox/PyBioSnippets/Lyrebird/fasta_libraries/GtRNAdb-all-tRNAs.fa",
    "illumna_kmer_tab_file": "?",
}


def get_kmers_from_monomer(monomer, k):
    """ Get kmers based on monomer cyclic rotations.
    """
    result = set()
    n = len(monomer)
    for i in xrange(n):
        if i:
            monomer = monomer[1:] + monomer[0]
        if k % n == 0:
            kmer = monomer * (k/n)
        else:
            kmer = monomer * (k/n) + monomer[:k%n]
        if i == 0:
            ref_kmer = kmer
        rkmer = get_revcomp(kmer)
        result.add(kmer)
        result.add(rkmer)
    return ref_kmer, result


def upload_perfect_microsatellite(connector, k):
    """
    """
    print("Generate monomers...")
    monomers = set()
    for a in 'ACTG':
        monomers.add(a)
        for b in 'ACTG':
            monomers.add(a+b)
            for c in 'ACTG':
                monomers.add(a+b+c)
                for d in 'ACTG':
                    monomers.add(a+b+c+d)
    type_name = "Microsatellites"
    seen_monomers = set()
    kid = 0
    n = len(monomers)
    for i, monomer in enumerate(monomers):
        print("Process %s monomer from %s" % (i, n), end=" ")
        if monomer in seen_monomers:
            continue
        variants = set(get_shifts_variants(monomer) + get_shifts_variants(get_revcomp(monomer)))
        lex_consensus = min(variants)
        family_name = "(%s)n" % lex_consensus
        for variant in variants:
            ref_kmer, kmers = get_kmers_from_monomer(variant, k)
            for kmer in kmers:
                rkmer = get_revcomp(kmer)
                if kmer > rkmer:
                    kmer = rkmer
                data = {
                    'kid': kid,
                    'kmer': kmer,
                    'inner_name': family_name,
                    'inner_family': type_name,
                    'ref_kmer': kmer,
                    'ed': 0,
                    'hd': 0,
                    'dust': get_dust_score(kmer),
                }
                result = connector.add_update(data)
                if result == "added":
                    kid += 1
    print


def upload_ribosomal_dna(connector, k):
    """
    """ 
    fasta_file = settings["rDNA_fasta"]
    id2doc = {}
    data = []
    for i, seq_obj in enumerate(sc_iter_fasta(fasta_file)):
        data.append(seq_obj.sequence)
        id2doc[i] = seq_obj.header
        connector.add_seq(seq_obj)
    kid = 1000000
    kmer_specter = process_list_to_kmer_index(data, k, docids=True, cutoff=None, verbose=False)
    result = {}
    family_name = "rDNA"
    for index in kmer_specter:
        names = set([id2doc[x] for x in index[4]])
        kmer = index[0]
        rkmer = index[1]
        if kmer > rkmer:
            kmer = rkmer
        data = {
                'kid': kid,
                'kmer': kmer,
                'inner_name': names,
                'inner_family': family_name,
                'ref_kmer': kmer,
                'ed': 0,
                'hd': 0,
                'dust': get_dust_score(kmer),
            }
        result = connector.add_update(data)
        if result == "added":
            kid += 1


def upload_mitochondrial_dna(connector, k):
    """
    """ 
    fasta_file = settings["mtDNA_fasta"]
    id2doc = {}
    data = []
    print("Load sequences...")
    for i, seq_obj in enumerate(sc_iter_fasta(fasta_file)):
        data.append(seq_obj.sequence)
        id2doc[i] = seq_obj.header
    kid = 2000000
    kmer_specter = process_list_to_kmer_index(data, k, docids=True, cutoff=None, verbose=False)
    result = {}
    family_name = "mtDNA"
    n = len(kmer_specter)
    for i, index in enumerate(kmer_specter):
        print("Upload kmers %s from %s..." % (i, n), end=" ")
        names = set([id2doc[x] for x in index[4]])
        kmer = index[0]
        rkmer = index[1]
        if kmer > rkmer:
            kmer = rkmer
        data = {
                'kid': kid,
                'kmer': kmer,
                'inner_name': names,
                'inner_family': family_name,
                'ref_kmer': kmer,
                'ed': 0,
                'hd': 0,
                'dust': get_dust_score(kmer),
            }
        result = connector.add_update(data, rewrite_name=False)
        if result == "added":
            kid += 1
    print


def upload_rfam(connector, k):
    """
    """ 
    fasta_file = settings["rfam_fasta"]
    id2doc = {}
    seqs = []
    print("Load sequences...")
    for i, seq_obj in enumerate(sc_iter_fasta(fasta_file)):
        seqs.append(seq_obj.sequence.replace("u", "t"))
        
        print(seq_obj.header, end=" ")
        famid, inner_name, pos = seq_obj.header.split(";")
        pos, taxon = pos.split(":", 1)

        data = {
                'sequence': seq_obj.sequence.replace("u", "t"),
                'full_name': seq_obj.header,
                'priority': 500,
                'family_name': famid,
                'file_name': "rfam",
                'inner_name': inner_name,
                'repbase_name': None,
                'taxons': taxon,
                'name_taxon': "%s:%s" % (inner_name, taxon),
            }

        id2doc[i] = "%s:%s" % (inner_name, taxon)       
        connector.add_rfam_sequence(data)
    print

    kid = 20000000
    kmer_specter = process_list_to_kmer_index(seqs, k, docids=True, cutoff=None, verbose=False)
    result = {}
    family_name = "Rfam"
    n = len(kmer_specter)
    objs = []
    print("Compute JSON dump...")
    for i, index in enumerate(kmer_specter):
        
        names = list(set([id2doc[x] for x in index[4]]))
        kmer = index[0]

        print("Process kmer %s (%s from %s) ..." % (kmer, i, n), end=" ")
        
        rkmer = index[1]
        if kmer > rkmer:
            kmer = rkmer
        data = {
                'kid': kid,
                'kmer': kmer,
                'inner_name': names,
                'inner_family': family_name,
                'ref_kmer': kmer,
                'ed': 0,
                'hd': 0,
                'dust': get_dust_score(kmer),
            }
        kid += 1

        objs.append(simplejson.dumps(data))

    print("Join data")
    data = "\n".join(objs)

    print("Save data")
    dataset = {
        "db_name": "Lyrebird",
        "collection_name": "RfamIndex",
    }
    dataset["json_file_to_upload"] = "upload_rfam_database.temp"
    with open(dataset["json_file_to_upload"], "w") as fh:
        fh.write(data)
    upload_command = "mongoimport -d %(db_name)s -c %(collection_name)s %(json_file_to_upload)s" % dataset
    runner.run(upload_command)
    

def upload_repbase_database(k):
    """
    """
    fasta_file = settings["repbase_fasta"]
    index_file = settings["repbase_index"]
    id2doc = {}
    id2meta = {}
    seqs = []
    print("Load sequences...")
    for i, seq_obj in enumerate(sc_iter_fasta(fasta_file)):
        seqs.append(seq_obj.sequence.replace("u", "t"))
    print("Load index...")
    with open(index_file) as fh:
        for line in fh:
            rep_id, rep_file, rep_name = line.strip().split("\t")
            rep_file = rep_file.split(".")[0]
            repbase_name = rep_name[1:]
            rep_name = rep_name[1:]
            try:
                inner_name, family_name, taxon = rep_name.split("    ")
                rep_name = "%s&%s&%s&%s" % (rep_file, inner_name, family_name, taxon)
                id2meta[int(rep_id)] = (rep_file, inner_name, family_name, taxon)
            except:
                rep_name = "%s&%s&%s&%s" % (rep_file, rep_name, rep_file, "None")
                id2meta[int(rep_id)] = (rep_file, inner_name, rep_file, "None")
            print(rep_name)
            id2doc[int(rep_id)] = rep_name

            data = {
                'sequence': seqs[int(rep_id)],
                'full_name': rep_name,
                'priority': 900,
                'family_name': id2meta[int(rep_id)][2],
                'file_name': id2meta[int(rep_id)][0],
                'inner_name': id2meta[int(rep_id)][1],
                'repbase_name': repbase_name,
                'taxons': id2meta[int(rep_id)][3],
            }

            connector.add_repbase_sequence(data)


    print("Process kmers...")
    kmer_specter = process_list_to_kmer_index(seqs, k, docids=True, cutoff=None, verbose=False)
    kid = 80000000
    objs = []
    family_name = "Repbase"
    for index in kmer_specter:

        names = list(set([id2doc[x] for x in index[4]]))
        kmer = index[0]
        rkmer = index[1]
        if kmer > rkmer:
            kmer = rkmer
        data = {
                'kid': kid,
                'kmer': kmer,
                'inner_name': names,
                'inner_family': family_name,
                'ref_kmer': kmer,
                'ed': 0,
                'hd': 0,
                'dust': get_dust_score(kmer),
            }
        kid += 1
   
        objs.append(simplejson.dumps(data))

    print("Join data")
    data = "\n".join(objs)

    print("Save data")
    dataset = {
        "db_name": "Lyrebird",
        "collection_name": "RepbaseIndexV2",
    }
    dataset["json_file_to_upload"] = "upload_repbase_database.temp"
    with open(dataset["json_file_to_upload"], "w") as fh:
        fh.write(data)
    upload_command = "mongoimport -d %(db_name)s -c %(collection_name)s %(json_file_to_upload)s" % dataset
    runner.run(upload_command)


def upload_repbase_database_13(connector):
    """
    """
    fasta_file = settings["repbase_fasta"]
    index_file = settings["repbase_index"]
    id2doc = {}
    id2meta = {}
    seqs = []
    k = 13
    print("Load sequences...")
    for i, seq_obj in enumerate(sc_iter_fasta(fasta_file)):
        seqs.append(seq_obj.sequence.replace("u", "t"))
    print("Load index...")
    with open(index_file) as fh:
        for line in fh:
            rep_id, rep_file, rep_name = line.strip().split("\t")
            rep_file = rep_file.split(".")[0]
            repbase_name = rep_name[1:]
            rep_name = rep_name[1:]
            try:
                inner_name, family_name, taxon = rep_name.split("    ")
                rep_name = "%s&%s&%s&%s" % (rep_file, inner_name, family_name, taxon)
                id2meta[int(rep_id)] = (rep_file, inner_name, family_name, taxon)
            except:
                rep_name = "%s&%s&%s&%s" % (rep_file, rep_name, rep_file, "None")
                id2meta[int(rep_id)] = (rep_file, inner_name, rep_file, "None")
            print(rep_name)
            id2doc[int(rep_id)] = rep_name

            data = {
                'sequence': seqs[int(rep_id)],
                'full_name': rep_name,
                'priority': 900,
                'family_name': id2meta[int(rep_id)][2],
                'file_name': id2meta[int(rep_id)][0],
                'inner_name': id2meta[int(rep_id)][1],
                'repbase_name': repbase_name,
                'taxons': id2meta[int(rep_id)][3],
            }

            connector.add_repbase_sequence(data)


    print("Process kmers...")
    kmer_specter = process_list_to_kmer_index(seqs, k, docids=True, cutoff=None, verbose=False)
    kid = 80000000
    objs = []
    family_name = "Repbase"
    for index in kmer_specter:

        names = list(set([id2doc[x] for x in index[4]]))
        kmer = index[0]
        rkmer = index[1]
        if kmer > rkmer:
            kmer = rkmer
        data = {
                'kid': kid,
                'kmer': kmer,
                'inner_name': names,
                'inner_family': family_name,
                'ref_kmer': kmer,
                'ed': 0,
                'hd': 0,
                'dust': get_dust_score(kmer),
            }
        kid += 1
   
        objs.append(simplejson.dumps(data))

    print("Join data")
    data = "\n".join(objs)

    print("Save data")
    dataset = {
        "db_name": "Lyrebird",
        "collection_name": "RepbaseIndexV2_k13",
    }
    dataset["json_file_to_upload"] = "upload_repbase_database_k13.temp"
    with open(dataset["json_file_to_upload"], "w") as fh:
        fh.write(data)
    upload_command = "mongoimport -d %(db_name)s -c %(collection_name)s %(json_file_to_upload)s" % dataset
    runner.run(upload_command)


def upload_illumina(connector, k=23, name="illumina", family_name="adapters"):
    """
    """
    illumna_kmer_tab_file = settings["illumna_kmer_tab_file"]

    obj = []
    print("Compute JSON dump...")
    with open(illumna_kmer_tab_file) as fh:
        for i, line in enumerate(fh):

            kid = -1000000

            kmer, tf = line.strip().split("\t")
            
            print("Process kmer %s (%s from ?) ..." % (kmer, i), end=" ")
            
            rkmer = get_revcomp(kmer)
            if kmer > rkmer:
                kmer = rkmer
            data = {
                    'kid': kid,
                    'kmer': kmer,
                    'inner_name': name,
                    'inner_family': family_name,
                    'ref_kmer': kmer,
                    'ed': 0,
                    'hd': 0,
                    'dust': get_dust_score(kmer),
                }
            result = connector.add_update(data)
            if result == "added":
                kid += 1
    print


def upload_sinebase(k=23):
    """
    1. Join repbase
        names -> sequence
    2. Compite kmers
    3. Upload to DB

        kmer
        rkmer
        repbase_index
        repbase_pos
        dust_score
        family
        superfamily

    """
    fasta_file = settings["sinebase_fasta"]
    id2doc = {}
    data = []
    for i, seq_obj in enumerate(sc_iter_fasta(fasta_file)):
        data.append(seq_obj.sequence.upper())
        id2doc[i] = seq_obj.header
    kmer_specter = process_list_to_kmer_index(data, k, docids=True, cutoff=None, verbose=False)


    result = {}
    for index in kmer_specter:
        names = [id2doc[x] for x in index[4]]
        data = {
            "kmer": index[0],
            "rkmer": index[1],
            "repbase_index": index[4],
            "repbase_names": names,
            "repbase_pos": None,
            "dust_score": get_dust_score(index[0]),
            "family": None,
            "superfamily": "SINE",
        }
        result[index[0].upper()] = data
        result[index[1].upper()] = data
    return result


def upload_trna(k=23):
    """
    """
    fasta_file = settings["tRNA_fasta"]
    id2doc = {}
    data = []
    for i, seq_obj in enumerate(sc_iter_fasta(fasta_file)):
        data.append(seq_obj.sequence.upper())
        id2doc[i] = seq_obj.header
    kmer_specter = process_list_to_kmer_index(data, k, docids=True, cutoff=None, verbose=False)

    result = {}
    for index in kmer_specter:
        names = [id2doc[x] for x in index[4]]
        data = {
            "kmer": index[0],
            "rkmer": index[1],
            "repbase_index": index[4],
            "repbase_names": names,
            "repbase_pos": None,
            "dust_score": get_dust_score(index[0]),
            "family": None,
            "superfamily": "tRNA",
        }
        result[index[0].upper()] = data
        result[index[1].upper()] = data
    return result


def upload_sequences(sequences, family, name, k=23, kid=-1):
    """
    """
    if isinstance(sequences, str):
        sequences = [sequences]
    connector = LyrebirdConnector()
    kmer_specter = process_list_to_kmer_index(sequences, k, docids=True, cutoff=None, verbose=False)
    result = {}
    for index in kmer_specter:
        kmer = index[0]
        rkmer = index[1]
        if kmer > rkmer:
            kmer = rkmer
        data = {
                'kid': kid,
                'kmer': kmer,
                'inner_name': name,
                'inner_family': family,
                'ref_kmer': kmer,
                'ed': 0,
                'hd': 0,
                'dust': get_dust_score(kmer),
            }
        result = connector.add_update(data)


if __name__ == '__main__':
    
    connector = LyrebirdConnector(human=True)
    k = 23
    
    upload_perfect_microsatellite(connector, k)
    # upload_ribosomal_dna(connector, k)
    # upload_mitochondrial_dna(connector, k)
    # upload_rfam(connector, k)
    # upload_repbase_database_13(connector)
    # upload_repbase_database(k)
    # upload_illumina(connector, k=23, name="illumina", family_name="adapters")
    # upload_sinebase(k=k)
    # upload_trna(k=k)


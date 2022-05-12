#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 21.05.2014
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com

import aindex
import argparse
# import jellyfish
import sys, re

sys.path.append(".")
sys.path.append("./trseeker")
sys.path.append("./aindex")


# from trseeker.tools.jellyfish_tools import Kmer2tfAPI

from collections import defaultdict, Counter
# from Levenshtein import distance, hamming
# from trseeker.tools.sequence_tools import get_dust_score
from tools.sequence_tools import get_revcomp

# from trseeker.tools.ngrams_tools import get_most_freq_kmer
# from trseeker.tools.raw_reads_tools import raw_reads_continue_kmer_right

# from PyBioSnippets.Lyrebird.tools.lyrebird_io import load_kmers

from tools.connectors import LyrebirdConnector
from tools.extenders import extend_kmer_right, extend_kmer_left
from tools.submitters import learn_add_te, learn_add_trs
from tools.adders import add_lyrebird_seq, add_kmers
# from PyBioSnippets.Lyrebird.tools.repbase import *

# from PyBioSnippets.Lyrebird.denovo_repeats import *


def solution_known_and_tr(c):

    length = len(c["sequence"])
    copy_number, size_bp = get_copy_number(kmer2tf, length, coverage, kmer, kmer2tf[kmer])
    status =  "%s:%sbp" % (status_right, length)
    print(i, tf, status, kmer, status_right, sequence, copy_number, round(size_bp/1000000.,3), "Mb")
    inner_name, inner_family, specific_kmers = learn_add_trs(i, sequence, kmer2tf, kmer, jellyfish_db, all_kmers, k, taxon, prefix, lyrebird)
    skipped = enrich_tr(sequence, k, skipped, specific_kmers)
    report[inner_name] = (i, tf, kmer, inner_family, inner_name, copy_number, round(size_bp/1000000.,3))

def solution_known_and_te(context):
    pass

def enrich_tr(monomer, k, skipped):
    """
    """
    sequence = monomer_to_sequence(monomer, k)
    kmer_for_adding = set()
    extend_hd0 = 0
    extend_hd2 = 0
    for i in xrange(len(sequence)-k+1):
        kmer = sequence[i:i+k]
        kmer_for_adding.add(kmer)
        kmer_for_adding.add(get_revcomp(kmer))
    # extend_hd0 = len(kmer_for_adding)
    # for kmer in kmer_for_adding:
    #     kmer = str(kmer)
    #     for other_kmer in all_kmers:
    #         other_kmer = str(other_kmer)
    #         if hamming(kmer, other_kmer) < 3:
    #             skipped.add(other_kmer)
    #             skipped.add(get_revcomp(other_kmer))
    #             extend_hd2 += 2
    for kmer in kmer_for_adding:
        skipped.add(kmer)
    print("\t\t", "Extended hd0=%s and hd2=%s" % (extend_hd0, extend_hd2))
    return skipped

def enrich_kmer(kmer, skipped):
    """
    """
    kmer_for_adding = set()
    extend_hd0 = 2
    extend_hd2 = 0    
    kmer_for_adding.add(kmer)
    kmer_for_adding.add(get_revcomp(kmer))
    for kmer in kmer_for_adding:
        kmer = str(kmer)
        for other_kmer in all_kmers:
            other_kmer = str(other_kmer)
            if hamming(kmer, other_kmer) < 3:
                skipped.add(other_kmer)
                skipped.add(get_revcomp(other_kmer))
                extend_hd2 += 2
    for kmer in kmer_for_adding:
        skipped.add(kmer)
    print("\t\t", "Extended hd0=%s and hd2=%s" % (extend_hd0, extend_hd2))
    return skipped

def enrich_te(sequence, k, skipped):
    """
    """
    kmer_for_adding = set()
    extend_hd0 = 0
    extend_hd2 = 0
    for i in xrange(len(sequence)-k+1):
        kmer = sequence[i:i+k]
        kmer_for_adding.add(kmer)
        kmer_for_adding.add(get_revcomp(kmer))
    extend_hd0 = len(kmer_for_adding)
    for kmer in kmer_for_adding:
        kmer = str(kmer)
        for other_kmer in all_kmers:
            other_kmer = str(other_kmer)
            if hamming(kmer, other_kmer) < 3:
                skipped.add(other_kmer)
                skipped.add(get_revcomp(other_kmer))
                extend_hd2 += 2
    for kmer in kmer_for_adding:
        skipped.add(kmer)
    print("\t\t", "Extended hd0=%s and hd2=%s" % (extend_hd0, extend_hd2))
    return skipped

def get_copy_number(length, coverage, kmer, tf):
    c = (len(kmer) - length)
    if c < 0:
        c = 1
    cn = (float(tf)/coverage) * c
    return cn, cn*length 


def report_update_with_known(i, report, known, kmer2tf, kmer, tr_family_to_letter, coverage):
    tf = kmer2tf[kmer]
    #TODO: implement for all other families other than microsatellites
    if known["inner_family"] == "SatDNA":
        length, letter = re.findall("(\d+)(\w)", known["inner_name"])[0]
        if not int(length) in tr_family_to_letter:
            tr_family_to_letter[int(length)] = 0
        else:
            tr_family_to_letter[int(length)] = max([ord(letter)-65, tr_family_to_letter[int(length)]])

        # print tr_family_to_letter
        # raw_input("?")

    if known["inner_family"].lower() in ["microsatellite", "telomeric", "satdna"]:
        if "(" in known["inner_name"]:
            length = len(known["inner_name"]) - 3
        else:
            print(known)
            length = int(re.findall("\d+", known["inner_name"])[0])
        copy_number, size_bp = get_copy_number(length, coverage, kmer, tf)
        if known["inner_name"].startswith("("):
            # skipped = enrich_tr(known["inner_name"][1:-2], k, skipped)
            pass
        # report[known["inner_name"]] = (i, tf, kmer, known["inner_family"], known["inner_name"], copy_number, round(size_bp/1000000.,3))
    elif known["inner_family"].lower() in ["satdna", "sine", "line"]:
        # skipped = enrich_kmer(kmer, skipped)
        pass
    else:
        print("?-->Unknown family %s" % known["inner_family"])
        # skipped = enrich_kmer(kmer, skipped)
    # raw_input("Go to the next?")






def extract_satellite_dna_candidates(settings):
    """ The main function.
    """
    statistics = { 
        "families_tf": defaultdict(int),
    }
    report = {}
    tr_family_to_letter = {}
    annotation_cache = {}
    te_id = 0

    force_kmer = settings["force_kmer"]
    min_tf = settings["min_tf"]
    verbose = settings["verbose"]

    coverage = settings["coverage"]
    k = settings["k"]
    taxon = settings["taxon"]
    prefix = settings["prefix"]
    
    # 1. Load kmers with aindex

    pf_settings = {
      "index_prefix": settings["aindex"],
      "aindex_prefix": None,
      "reads_file": None,
    }
    kmer2tf = aindex.load_aindex(pf_settings, skip_reads=True, skip_aindex=True)

    # 2. Iterate over kmers
    skipped = set()
    lyrebird = LyrebirdConnector(local=settings["local_lb"])

    kmer_fh = open(settings["kmer_file"])

    inner_name = None
    with open(settings["output_file"], "w") as fh:
        pass
    for i, line in enumerate(kmer_fh):
        if inner_name:
            with open(settings["output_file"], "a") as fh:
                (i, tf, kmer, inner_family, inner_name, copy_number, amount, sequence) = report[inner_name]
                header = ":".join(map(str, ("gid%s" % gid, i, tf, kmer, inner_family, inner_name, copy_number, amount)))
                fh.write(">%s\n" % header)
                fh.write("%s\n" % sequence)
            inner_name = None
        if i < settings["start"]:
            continue
        if force_kmer == "Done":
            break
        if force_kmer:
            kmer = force_kmer
            tf = kmer2tf[kmer]
            force_kmer = "Done"
        else:
            if "\t" in line:
                gid, kmer, tf = line.strip().split("\t")
            else:
                gid, kmer, tf = line.strip().split()
        tf = int(tf)
        if tf < min_tf:
            print("Stopped by tf %s (%s)" % (tf, min_tf))
            break
        if verbose:
            print("Process %s kmer: %s with tf=%s (cov=%s) kmer=%s" % (i, kmer, tf, tf/coverage, kmer))
        
        # 2.1. check previously seen
        if kmer in skipped:
            if verbose > 2:
                print("\t", "previously seen, skipped")
            continue

        # 2.2. check previously known
        known = lyrebird.query(kmer)
        if known:
            annotation_cache[kmer] = known
            if verbose > 10:
                print("\t", "previously known: %s | %s" % (known["inner_family"], known["inner_name"]))
            report_update_with_known(i, report, known, kmer2tf, kmer, tr_family_to_letter, coverage)
            continue

        # 2.3. check previously unknown
        if verbose > 10:
            print("\tpreviously unknown, computing...")
            print("\t\tExtended right:", end=" ")
        status_right, sequence, length, path, known_meta_right = extend_kmer_right(kmer, kmer2tf, tf_cutoff=settings["tf_cutoff"], fraq_cutoff=settings["fraq_cutoff"], lyrebird=lyrebird)
        # print status_right, sequence, length
        
        if status_right == "TRs":

            length = len(sequence)
            copy_number, size_bp = get_copy_number(length, coverage, kmer, tf)
            status =  "%s:%sbp" % (status_right, length)
            
            if verbose > 10:
                print("\t\t", "Copy number:", copy_number, round(size_bp/1000000.,3), "Mb")

            if length < 100:
                inner_name, inner_family, specific_kmers = learn_add_trs(i, sequence, kmer2tf, kmer, k, taxon, prefix, lyrebird, tr_family_to_letter, annotation_cache, fast=settings["fast_trs_adding"])
            else:
                inner_name, inner_family, specific_kmers = learn_add_trs(i, sequence, kmer2tf, kmer, k, taxon, prefix, lyrebird, tr_family_to_letter, annotation_cache, fast=settings["fast_satdna_adding"])

            # skipped = enrich_tr(sequence, k, skipped)
            report[inner_name] = (i, tf, kmer, inner_family, inner_name, copy_number, round(size_bp/1000000.,3), sequence)
            continue

        if verbose > 10:
            print("\t\t", "Extended left: ", sequence)
        status_left, sequence, length, left_path, known_meta_left = extend_kmer_left(sequence, status_right, kmer2tf, k, tf_cutoff=settings["tf_cutoff"], fraq_cutoff=settings["fraq_cutoff"], lyrebird=lyrebird)
        length = len(sequence)
        if verbose > 10:
            print(sequence, length)
        print("STATUSES:", status_left, status_right)

        copy_number = tf / coverage
        size = len(sequence) * copy_number
        
        def add_default_dont_ask(lyrebird, sequence, taxon, inner_name, inner_family, priority, delta, status_left, status_right, k=23):
            
            if verbose > 10:
                print("New length: %s" % (len(sequence)-2*k))
                print("Sequence was added as: %s %s %s (%s %s) [add]" % (sequence, inner_name, inner_family, status_left, status_right))
            add_lyrebird_seq(lyrebird, sequence, taxon, inner_name, inner_family)
            length = len(sequence)
            add_kmers(lyrebird, sequence, inner_name, inner_family, priority, k=k, prev_known=None, rewrite_name=False, delta=delta)
            return True

        def add_default_or_ask(lyrebird, sequence, taxon, inner_name, inner_family, priority, delta, status_left, status_right, k=23):

            if verbose > 10:
                print("New length: %s" % (len(sequence)-2*k))
            if not settings["fast_te_adding"]:
                answer = raw_input("Sequence was added as: %s %s %s (%s %s) [add]" % (sequence, inner_name, inner_family, status_left, status_right)) or None
            else:
                answer = None

            if answer is None:
                add_lyrebird_seq(lyrebird, sequence, taxon, inner_name, inner_family)
                length = len(sequence)
                add_kmers(lyrebird, sequence, inner_name, inner_family, priority, k=k, prev_known=None, rewrite_name=False, delta=delta)
                return True
            return False
            
        size_bp = len(sequence)
        
        if len(sequence) > 500:
            inner_name, inner_family, specific_kmers = learn_add_te(i, sequence, kmer2tf, kmer, k, taxon, prefix, lyrebird, annotation_cache, default_name="TE%s" % te_id, default_family="TE", fast_te_adding=settings["fast_te_adding"])
            te_id += 1

            report[inner_name] = (i, tf, kmer, inner_family, inner_name, copy_number, round(size_bp/1000000.,3), sequence)
            continue

        if status_right == "Known" and status_left == "Known":

            if verbose < 100:
                add_default_dont_ask(lyrebird, sequence, taxon, known_meta_right["inner_name"], known_meta_right["inner_family"], known_meta_right["priority"]-10, length-2*k, status_left, status_right, k=23)
            else:
                add_default_or_ask(lyrebird, sequence, taxon, known_meta_right["inner_name"], known_meta_right["inner_family"], known_meta_right["priority"]-10, length-2*k, status_left, status_right, k=23)

            report[inner_name] = (i, tf, kmer, inner_family, inner_name, copy_number, round(size_bp/1000000.,3), sequence)
            continue
            
        
        if (status_right == "Known" and status_left == "Loop") or (status_right == "Loop" and status_left == "Known") or (status_right == "Known" and status_left == "Known"):

            if not known_meta_left:
                add_default_dont_ask(lyrebird, sequence, taxon, known_meta_right["inner_name"], known_meta_right["inner_family"], known_meta_right["priority"]-10, length-2*k, status_left, status_right, k=23)
                # if not add_default_or_ask(sequence, taxon, known_meta_right["inner_name"], known_meta_right["inner_family"], known_meta_right["priority"]-10, length-2*k, status_left, status_right, k=23):
                #     inner_name, inner_family, specific_kmers = learn_add_te(i, sequence, kmer2tf, kmer, k, taxon, prefix, lyrebird, annotation_cache)
                report[inner_name] = (i, tf, kmer, inner_family, inner_name, copy_number, round(size_bp/1000000.,3), sequence)
                continue
            elif not known_meta_right:
                add_default_dont_ask(lyrebird, sequence, taxon, known_meta_left["inner_name"], known_meta_left["inner_family"], known_meta_left["priority"]-10, length-2*k, status_left, status_right, k=23)
                # if not add_default_or_ask(sequence, taxon, known_meta_left["inner_name"], known_meta_left["inner_family"], known_meta_left["priority"]-10, length-2*k, status_left, status_right, k=23):
                #     inner_name, inner_family, specific_kmers = learn_add_te(i, sequence, kmer2tf, kmer, k, taxon, prefix, lyrebird, annotation_cache)
                report[inner_name] = (i, tf, kmer, inner_family, inner_name, copy_number, round(size_bp/1000000.,3), sequence)
                continue

            if known_meta_left["inner_name"] == known_meta_right["inner_name"]:
                
                add_default_dont_ask(lyrebird, sequence, taxon, known_meta_right["inner_name"], known_meta_right["inner_family"], known_meta_right["priority"]-10, length-2*k, status_left, status_right, k=23)
                
                # else:
                #     if not add_default_or_ask(sequence, taxon, known_meta_right["inner_name"], known_meta_right["inner_family"], known_meta_right["priority"]-10, length-2*k, status_left, status_right, k=23):
                #         inner_name, inner_family, specific_kmers = learn_add_te(i, sequence, kmer2tf, kmer, k, taxon, prefix, lyrebird, annotation_cache)
                report[inner_name] = (i, tf, kmer, inner_family, inner_name, copy_number, round(size_bp/1000000.,3), sequence)
                continue
            
            else:

                if "fusion" in known_meta_left["inner_family"]:
                    if known_meta_right["inner_name"] in known_meta_left["inner_name"]:
                        add_default_dont_ask(lyrebird, sequence, taxon, known_meta_left["inner_name"], known_meta_left["inner_family"], known_meta_left["priority"]-10, length-2*k, status_left, status_right, k=23)
                        # if not add_default_or_ask(sequence, taxon, known_meta_left["inner_name"], known_meta_left["inner_family"], known_meta_left["priority"]-10, length-2*k, status_left, status_right, k=23):
                        #     inner_name, inner_family, specific_kmers = learn_add_te(i, sequence, kmer2tf, kmer, k, taxon, prefix, lyrebird, annotation_cache)
                        report[inner_name] = (i, tf, kmer, inner_family, inner_name, copy_number, round(size_bp/1000000.,3), sequence)
                        continue
                    else:
                        print(status_right, status_left)
                        if known_meta_left:
                            print(known_meta_left["inner_name"], known_meta_left["inner_family"])
                        if known_meta_right:
                            print(known_meta_right["inner_name"], known_meta_right["inner_family"])
                        report[inner_name] = (i, tf, kmer, inner_family, inner_name, copy_number, round(size_bp/1000000.,3), sequence)    
                        continue
                        raw_input("what R in L fusion?")
                elif "fusion" in known_meta_right["inner_family"]:
                    if known_meta_left["inner_name"] in known_meta_right["inner_name"]:
                        add_default_dont_ask(lyrebird, sequence, taxon, known_meta_right["inner_name"], known_meta_right["inner_family"], known_meta_right["priority"]-10, length-2*k, status_left, status_right, k=23)
                        # if not add_default_or_ask(sequence, taxon, known_meta_right["inner_name"], known_meta_right["inner_family"], known_meta_right["priority"]-10, length-2*k, status_left, status_right, k=23):
                        #     inner_name, inner_family, specific_kmers = learn_add_te(i, sequence, kmer2tf, kmer, k, taxon, prefix, lyrebird, annotation_cache)
                        report[inner_name] = (i, tf, kmer, inner_family, inner_name, copy_number, round(size_bp/1000000.,3), sequence)
                        continue
                    else:
                        print(status_right, status_left)
                        if known_meta_left:
                            print(known_meta_left["inner_name"], known_meta_left["inner_family"])
                        if known_meta_right:
                            print(known_meta_right["inner_name"], known_meta_right["inner_family"])
                        report[inner_name] = (i, tf, kmer, inner_family, inner_name, copy_number, round(size_bp/1000000.,3), sequence)
                        continue
                        raw_input("what L in R fusion?")
                else:
                    inner_name = "%s:%s" % (known_meta_left["inner_name"], known_meta_right["inner_name"])
                    inner_family = "%s:%s:fusion" % (known_meta_left["inner_family"], known_meta_right["inner_family"])

                    if not add_default_or_ask(lyrebird, sequence, taxon, inner_name, inner_family, 50, 0, status_left, status_right, k=23):
                        inner_name, inner_family, specific_kmers = learn_add_te(i, sequence, kmer2tf, kmer, k, taxon, prefix, lyrebird, annotation_cache)
                    report[inner_name] = (i, tf, kmer, inner_family, inner_name, copy_number, round(size_bp/1000000.,3), sequence)
                    continue


        if (status_right == "Zero_cov" and status_left == "Zero_cov") or (status_right == "Zero_cov" and status_left == "Loop") or (status_right == "Loop" and status_left == "Zero_cov"):

            length = len(sequence)
            copy_number, size_bp = get_copy_number(length, coverage, kmer, tf)
            inner_name, inner_family, specific_kmers = learn_add_te(i, sequence, kmer2tf, kmer, k, taxon, prefix, lyrebird, annotation_cache, default_name="TE%s" % te_id, default_family="TE", fast_te_adding=settings["fast_te_adding"])
            te_id += 1
            print("\t\t", "Copy number:", copy_number, round(size_bp/1000000.,3), "Mb")
            report[inner_name] = (i, tf, kmer, inner_family, inner_name, copy_number, round(size_bp/1000000.,3), sequence)
            continue
            
        elif (status_left == "Zero_cov" and status_right == "Known"):

            inner_name = "left_end:%s" % (known_meta_right["inner_name"])
            inner_family = known_meta_right["inner_family"]

            if not add_default_or_ask(lyrebird, sequence, taxon, inner_name, inner_family, known_meta_right["priority"]-50, 0, status_left, status_right, k=23):
                inner_name, inner_family, specific_kmers = learn_add_te(i, sequence, kmer2tf, kmer, k, taxon, prefix, lyrebird, annotation_cache, default_name="TE%s" % te_id, default_family="TE", fast_te_adding=settings["fast_te_adding"])
                te_id += 1
            report[inner_name] = (i, tf, kmer, inner_family, inner_name, copy_number, round(size_bp/1000000.,3), sequence)
            continue

        elif (status_left == "Known" and status_right == "Zero_cov"):

            inner_name = "%s:right_end" % (known_meta_left["inner_name"])
            inner_family = known_meta_left["inner_family"]

            if not add_default_or_ask(lyrebird, sequence, taxon, inner_name, inner_family, known_meta_left["priority"]-50, 0, status_left, status_right, k=23):
                inner_name, inner_family, specific_kmers = learn_add_te(i, sequence, kmer2tf, kmer, k, taxon, prefix, lyrebird, annotation_cache, default_name="TE%s" % te_id, default_family="TE", fast_te_adding=settings["fast_te_adding"])
                te_id += 1
            report[inner_name] = (i, tf, kmer, inner_family, inner_name, copy_number, round(size_bp/1000000.,3), sequence)
            continue

        elif (status_left == "Loop" and status_right == "Loop"):

            inner_family = "Unknown"
            inner_name = "loop:loop:%s" % (sequence)
            size_bp = 0

            print(status_right, status_left)
            if known_meta_left:
                print(known_meta_left["inner_name"], known_meta_left["inner_family"])
            if known_meta_right:
                print(known_meta_right["inner_name"], known_meta_right["inner_family"])
            report[inner_name] = (i, tf, kmer, inner_family, inner_name, copy_number, round(size_bp/1000000.,3), sequence)

            if verbose > 1000:
                print(report[inner_name])

            continue

        else:
            print(status_right, status_left)
            if known_meta_left:
                print(known_meta_left["inner_name"], known_meta_left["inner_family"])
            if known_meta_right:
                print(known_meta_right["inner_name"], known_meta_right["inner_family"])
            raw_input("what?")

        while True:
            d = raw_input("Sequence found, check it with blast and giri (c - to continue)").strip()
            if d == "c":
                break
            continue

        if status_right == "Loop":
            loop_kmer = length[1]
            if loop_kmer in skipped:
                print("Loop with prev seen", kmer, tf, loop_kmer)
                continue
            else:
                loop_kmer = length[1]

                print(kmer, tf, status_right, sequence, length, path)

                print(left_status, left_sequence, left_length, left_path)
                for x in path:
                    if x[2][-1][-2] in skipped:
                        print("SKIP HIT:")
                    if x[2][-1][-2] != loop_kmer:
                        print(x[1], x[2][-1])
                    else:
                        print(x[1], x[2][-1], "<--")

                raise Exception("Unknown loop structure")
        else:
            print(kmer, tf, status_right, sequence, length, path)

            raise Exception("Other status: %s" % status_right)

            print(i, kmer, status, sequence, length)

            enrich_te(sequence, k, skipped)

            status =  "%s:%sbp" % (status, length)

        print(i, tf, kmer, known)
        raw_input("Next?")
        continue

    
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Extract tandem repeats candidates from raw reads.')
    parser.add_argument('-a', '--aindex', help='Aindex prefix', required=True)
    parser.add_argument('-i', '--input', help='Kmer sorted by freq from jellyfish dump -c -t command', required=True)
    parser.add_argument('-o', '--output', help='Output prefix', required=True)
    parser.add_argument('-p', '--prefix', help='Satellite naming suffix, e.g. Hsap or HS for hyman', required=True)
    parser.add_argument('-k', help='Kmer length, default 23', required=False, default=23)
    parser.add_argument('--kmer', help='Compute for exact kmer', required=False, default=None)
    parser.add_argument('-n', help='Number of families in output, default 31', required=False, default=31)
    parser.add_argument('-v', '--verbose', help='Verbose level', required=False, default=0)
    parser.add_argument('--coverage', help='Raw reads coverage for copy number estimation, default 1', required=False, default=1)
    args = vars(parser.parse_args())

    settigns = {
        "kmer_file": args["input"],
        "aindex": args["aindex"],
        "output_file": args["output"],
        "prefix": args["prefix"],
        "coverage": float(args["coverage"]),
        "n_families": int(args["n"]),
        "k": int(args["k"]),
        "v": int(args["verbose"]),
        "min_tf": float(args["coverage"]) * 100,
        "force_kmer": None,
        "start": 0,
        "tf_cutoff": float(args["coverage"]) * 100,
        "use_jf_db": False,
        "fraq_cutoff": 0.2,
        "fast_trs_adding": True,
        "fast_te_adding": True,
        "fast_satdna_adding": True,
        "local_lb": True,
        "dont_skip": True,
        "taxon": "NA",
        "prefix": "NA",
        "verbose": True,
    }

    extract_satellite_dna_candidates(settigns)

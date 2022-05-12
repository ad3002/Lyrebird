#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 07.09.2011
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com 

from PyExp import AbstractModel
from trseeker.seqio.tab_file import sc_iter_tab_file


class BlastResultModel(AbstractModel):
    """ Class for blast result data.

    Attributes:

    - "query_id" (int),
    - "query_gi" (int),
    - "query_ref",
    - "subject_id",
    - "subject_gi"(int),
    - "subject_ref",
    - "query_start" (int),
    - "query_end" (int),
    - "subject_start" (int),
    - "subject_end" (int),
    - "evalue"  (float),
    - "bit_score" (flaot),
    - "score" (int),
    - "alignment_length" (int),
    - "proc_identity" (float),
    - "identical" (int),
    - "mismatches" (int),
    - "positives" (int),
    - "gap_opens" (int),
    - "gaps" (int),
    - "proc_positives" (float),
    - "frames",
    - "query_frame" (int),
    - "subject_frame" (int),
    - "fraction_of_query" (float),    

    """

    dumpable_attributes = ["query_id",
                           "query_gi",
                           "query_ref",
                           "subject_id",
                           "subject_gi",
                           "subject_ref",
                           "query_start",
                           "query_end",
                           "subject_start",
                           "subject_end",
                           "evalue",
                           "bit_score",
                           "score",
                           "alignment_length",
                           "proc_identity",
                           "identical",
                           "mismatches",
                           "positives",
                           "gap_opens",
                           "gaps",
                           "proc_positives",
                           "frames",
                           "query_frame",
                           "subject_frame",
                           "fraction_of_query",
                           ]

    int_attributes = ["query_start",
                       "query_end",
                       "subject_start",
                       "subject_end",
                       "score",
                       "alignment_length",
                       "identical",
                       "mismatches",
                       "positives",
                       "gap_opens",
                       "gaps",
                       "query_frame",
                       "subject_frame",

                       ]

    float_attributes = [ "evalue",
                         "proc_identity",
                         "proc_positives",
                         "fraction_of_query",
                         "bit_score",
                       ]


def read_blast_file(blast_file, length):
    """ Read blast file. Return subject_ref -> list of matches (BlastResultModel models)."""
    #TODO: move it to readers
    
    gi_to_results = {}

    # remove # lines from blast file.
    data = []
    with open(blast_file) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            data.append(line)
    if not data:
        return gi_to_results
    data[-1] = data[-1].strip()
    with open(blast_file, "w") as fh:
        fh.writelines(data)

    # parser data
    for blast_obj in sc_iter_tab_file(blast_file, BlastResultModel):
        if blast_obj.query_end is None or blast_obj.query_end is None:
            print "Error parsing blast results:", blast_obj
            continue
        if length:
          blast_obj.fraction_of_query = abs(blast_obj.query_start - blast_obj.query_end) / float(length)
        subject_ref = blast_obj.subject_ref
        gi_to_results.setdefault(subject_ref, [])
        gi_to_results[subject_ref].append(blast_obj)
    return gi_to_results

#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 05.06.2011
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com
'''
Working with BLAST tab-delimited output files.
'''
import os
from trseeker.models.blast_model import read_blast_file
from PyExp import core_logger, AbstractModel
from collections import defaultdict
from trseeker.seqio.tab_file import sc_iter_tab_file


def _get_blast_result_intervals(blast_file, length, min_score):
    ''' Return gi->sorted list of blast objects. '''

    # read data
    gi_to_results = read_blast_file(blast_file, length)

    # create correct start/end positions
    for gi in gi_to_results:
        for i, blast_obj in enumerate(gi_to_results[gi]):
            # skip by score value
            if blast_obj.score < min_score:
                gi_to_results[gi][i] = None
                continue
            if blast_obj.query_start > blast_obj.query_end:
                temp = blast_obj.query_start
                gi_to_results[gi][i].query_start = blast_obj.query_end
                gi_to_results[gi][i].query_end = temp
        gi_to_results[gi] = [x for x in gi_to_results[gi] if x]

    # sort data
    for gi in gi_to_results:
        gi_to_results[gi].sort(key=lambda x: (x.query_start, x.query_end))

    return gi_to_results

def _remove_nested_join_overlap(blast_dataset):
    ''' Remove nested matches. 
        A: skip inners
            # -------
            #  ----
            # 
            # -------
            # ----
            # 
            # -------
            #    ---- 
            # 
            # -------
            # -------
        B: join overlapped
            # -------      -----
            #   ---------       -----
            
    '''
    for gi in blast_dataset:
        previous_item = None
        for i, item in enumerate(blast_dataset[gi]):
            if item is None:
                continue
            # check data
            assert item.query_start < item.query_end
            if not previous_item:
                previous_item = item
                pi = i
                continue
            # A: skip inners
            if item.query_start >= previous_item.query_start \
                    and item.query_end <= previous_item.query_end:
                blast_dataset[gi][i] = None
                continue

            # B: join overlapped
            if item.query_start <= previous_item.query_end \
                    and item.query_end > previous_item.query_end:

                blast_dataset[gi][pi].query_end = item.query_end
                previous_item.query_end = item.query_end

                blast_dataset[gi][pi].changed = True
                blast_dataset[gi][i] = None
                continue
            previous_item = item
            pi = i
        blast_dataset[gi] = [x for x in blast_dataset[gi] if x]
    return blast_dataset

def _join_gapped(blast_dataset, gap_size, min_align,):
    ''' Join gapped matches.
        C: remove short alignments 
        D: join gapped
        #  ------ 
        #            ------
    '''

    for gi in blast_dataset:
        previous_item = None
        for i, item in enumerate(blast_dataset[gi]):
            if item is None:
                continue

            # C: remove short alignments 
            length = item.query_end - item.query_start
            if length < min_align:
                blast_dataset[gi][i] = None
                continue

            if not previous_item:
                previous_item = item
                pi = i
                continue

            # D: join gapped
            if abs(previous_item.query_end - item.query_start) <= gap_size:
                blast_dataset[gi][pi].query_end = item.query_end
                previous_item.query_end = item.query_end
                blast_dataset[gi][pi].changed = True
                blast_dataset[gi][pi].gapped = True
                blast_dataset[gi][i] = None
                continue

            previous_item = item
            pi = i

        blast_dataset[gi] = [x for x in blast_dataset[gi] if x]
    return blast_dataset

def _filter_blast_dataset(blast_dataset, min_length):
    ''' Filter by minimal alignment length.

    - E: skip short aligns
    '''
    # E: skip short aligns 
    for gi in blast_dataset:
        for i, item in enumerate(blast_dataset[gi]):
            length = item.query_end - item.query_start
            if length < min_length:
                blast_dataset[gi][i] = None
                continue
        blast_dataset[gi] = [x for x in blast_dataset[gi] if x]
    blast_dataset = dict([(k, v) for k, v in blast_dataset.items() if v])
    return blast_dataset

def _format_output(blast_dataset, format_function):
    ''' Format blast_dataset output.'''
    assert hasattr(format_function, "__call__")
    return format_function(blast_dataset)

def format_function_repbase(blast_dataset):
    ''' Join with keys with ;.'''
    return ";".join(blast_dataset.keys())

def format_function_self(blast_dataset):
    ''' Join keys with ,.'''
    return ",".join(blast_dataset.keys())

def get_blast_result(blast_file, length, gap_size=1000, min_align=500, min_length=2400, format_function=None, min_score=90):
    """ Get blast result."""

    blast_dataset = _get_blast_result_intervals(blast_file, length, min_score)
    blast_dataset = _remove_nested_join_overlap(blast_dataset)
    blast_dataset = _join_gapped(blast_dataset, gap_size, min_align)
    blast_dataset = _filter_blast_dataset(blast_dataset, min_length)
    if not format_function:
        format_function = format_function_self
    result = _format_output(blast_dataset, format_function)
    return result

def update_with_repbase_blast_result(trs_dataset, annotation_self_folder, filters):
    ''' Add Repbase blast result to trs_dataset.'''
    for i, trf_obj in enumerate(trs_dataset):

        blast_output_file = os.path.join(annotation_self_folder, "%s.blast" % trf_obj.trf_id)
        result = get_blast_result(blast_output_file,
                                  trf_obj.trf_array_length,
                                  gap_size=filters["blast_gap_size"],
                                  min_align=filters["min_align"],
                                  min_length=filters["min_length"],
                                  format_function=format_function_repbase)
        trs_dataset[i].trf_repbase = result

    return trs_dataset

def update_with_self_blast_result(trs_dataset, annotation_self_folder, filters_obj):
    ''' Add vs_self blast result to trs_dataset.'''
    n = len(trs_dataset)
    for i, trf_obj in enumerate(trs_dataset):
        print i, n, "\r",
        blast_output_file = os.path.join(annotation_self_folder, "%s.blast" % trf_obj.trf_id)
        with open(blast_output_file) as fh:
            data = fh.read()
        if data.startswith("ALPHA"):
            data = data.split()
            ref = int(data[1])
            trs_dataset[i].trf_family_self = ('ALPHA', ref)
            continue
        filters = filters_obj.get_filters(trf_obj.trf_array_length)
        result = get_blast_result(blast_output_file,
                                  trf_obj.trf_array_length,
                                  gap_size=filters["blast_gap_size"],
                                  min_align=filters["min_align"],
                                  min_length=filters["min_length"],
                                  format_function=format_function_self,
                                  min_score=filters["min_blast_score"]
                                  )
        trs_dataset[i].trf_family_self = result
    print
    return trs_dataset

def update_with_ref_blast_result(trs_dataset, annotation_self_folder, filters):
    ''' Add vs_reference blast result to trs_dataset.'''
    n = len(trs_dataset)
    for i, trf_obj in enumerate(trs_dataset):
        print i, n, "\r",

        blast_output_file = os.path.join(annotation_self_folder, "%s.blast" % trf_obj.trf_id)
        result = get_blast_result(blast_output_file,
                                  trf_obj.trf_array_length,
                                  gap_size=filters["blast_gap_size"],
                                  min_align=filters["min_align"],
                                  min_length=filters["min_length"],
                                  format_function=format_function_self)
        trs_dataset[i].trf_family_ref = result
    print
    return trs_dataset




def parse_blast_tsv_file_wo_taxonomy(file_name, ref2length=None, add_func=None):
    ''' Parse blast tsv-file w/o or without taxonomy data.
    '''
    core_logger.info("Parsing: %s" % file_name)
    with open(file_name) as fh:
        for line in fh:
            if line.startswith("# Fields:"):
                fields = [x.strip().replace(" ","_").replace(".","").replace("%","p").replace("/","_") 
                            for x in line.split("# Fields:")[1].strip().split(",")]
                break
    
    class BlastData(AbstractModel):
        dumpable_attributes = fields

    def check_consistency(x):
        return len(x.split("\t")) == len(fields)

    hits = defaultdict(list)
    for i, tab_obj in enumerate(sc_iter_tab_file(file_name, BlastData, remove_starts_with="#", check_function=check_consistency)):
        print i, "\r",
        hit = {
            "hid": tab_obj.subject_id,
            "full_hid": tab_obj.subject_id,
            "S1": int(tab_obj.q_start),
            "E1": int(tab_obj.q_end),
            "S2": int(tab_obj.s_start),
            "E2": int(tab_obj.s_end),
            "evalue": float(tab_obj.evalue),
            "alignment_length": int(tab_obj.alignment_length),
            "p_identity": float(tab_obj.p_identity),
            "bit_score": int(round(float(tab_obj.bit_score),0)),
            "score": int(tab_obj.score),
            "identical": int(tab_obj.identical),
            "mismatches": int(tab_obj.mismatches),
            "positives": int(tab_obj.positives),
            "gap_opens": int(tab_obj.gap_opens),
            "gaps": int(tab_obj.gaps),
            "p_positives": float(tab_obj.p_positives),
            "query_frame": int(tab_obj.query_frame),
            "sbjct_frame": int(tab_obj.sbjct_frame),
        }
        if ref2length:
            hit["coverage"] = float(hit["alignment_length"])/ref2length[tab_obj.query_id]

        if add_func and hasattr(add_func, "__call__"):
            hit = add_func(hit, tab_obj)

        hits[tab_obj.query_id].append(hit)
    print
    return hits

def check_consistency(x):
    return len(x.split("\t")) == len(fields)

def iter_blast_tsv_file_wo_taxonomy(file_name, ref2length=None, add_func=None, check_consistency=None):
    ''' Parse blast tsv-file w/o or without taxonomy data.
    '''
    core_logger.info("Parsing: %s" % file_name)
    fields =None
    with open(file_name) as fh:
        for line in fh:
            if line.startswith("# Fields:"):
                fields = [x.strip().replace(" ","_").replace(".","").replace("%","p").replace("/","_") 
                            for x in line.split("# Fields:")[1].strip().split(",")]
                break
    if fields is None:
        raise Exception("Missed field in blast output.")
    
    class BlastData(AbstractModel):
        dumpable_attributes = fields

    hits = []
    last_obj = None
    for i, tab_obj in enumerate(sc_iter_tab_file(file_name, BlastData, skip_starts_with="#", check_function=check_consistency)):
        
        if last_obj and last_obj != tab_obj.query_id:
            yield last_obj, hits
            hits = []
            last_obj = tab_obj.query_id

        if last_obj is None:
            last_obj = tab_obj.query_id

        hit = {
            "query": tab_obj.query_id,
            "hid": tab_obj.subject_id,
            "full_hid": tab_obj.subject_id,
            "S1": int(tab_obj.q_start),
            "E1": int(tab_obj.q_end),
            "S2": int(tab_obj.s_start),
            "E2": int(tab_obj.s_end),
            "evalue": float(tab_obj.evalue),
            "alignment_length": int(tab_obj.alignment_length),
            "p_identity": float(tab_obj.p_identity),
            "bit_score": int(round(float(tab_obj.bit_score),0)),
            "score": int(tab_obj.score),
            "identical": int(tab_obj.identical),
            "mismatches": int(tab_obj.mismatches),
            "positives": int(tab_obj.positives),
            "gap_opens": int(tab_obj.gap_opens),
            "gaps": int(tab_obj.gaps),
            "p_positives": float(tab_obj.p_positives),
            "query_frame": int(tab_obj.query_frame),
            "sbjct_frame": int(tab_obj.sbjct_frame),
        }
        if ref2length:
            hit["coverage"] = float(hit["alignment_length"])/ref2length[tab_obj.query_id]

        if add_func and hasattr(add_func, "__call__"):
            hit = add_func(hit, tab_obj)

        hits.append(hit)

        
    if hits:
        yield last_obj, hits


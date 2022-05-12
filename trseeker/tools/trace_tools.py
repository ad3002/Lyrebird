#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 28.11.2011
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com
'''
Functions related to Trace data.
'''

import os
from trseeker.seqio.tab_file import sc_iter_simple_tab_file
from trseeker.seqio.fasta_file import sc_iter_fasta
import gzip

def get_clip_data(clip_file_name):
    ''' Read clip dictionary from clip_file_name gz archive.
    '''
    id2clip = {}
    fh = gzip.open(clip_file_name, 'rb')
    file_content = fh.readlines()
    for line in file_content:
        if not line:
            continue
        if line.startswith("TI"):
            continue
        sid, start, end = map(int, line.strip().split())
        id2clip[int(sid)] = (start, end)
    fh.close()
    return id2clip

def unclip_trace_file(fasta_file, clip_file, uncliped_file, trash_file):
    ''' Unclip. Remove vector flanks from Trace data.'''
    id2clip = get_clip_data(clip_file)
    result = []
    with open(trash_file, "a") as fh:
        for seq_obj in sc_iter_fasta(fasta_file):
            id = int(seq_obj.seq_gi)
            seq = seq_obj.sequence
            if id in id2clip:
                left_flank = seq_obj.sequence[id2clip[id][1]:]
                right_flank = seq_obj.sequence[:id2clip[id][0]]
                seq_obj.seq_sequence = seq_obj.sequence[id2clip[id][0]: id2clip[id][1]]
                if left_flank:
                    fh.write(">%s_LF\n%s\n" % (id, left_flank))
                if right_flank:
                    fh.write(">%s_RF\n%s\n" % (id, right_flank))
            result.append(seq_obj.fasta)

            

    with open(uncliped_file, "w") as fh:
        fh.writelines(result)

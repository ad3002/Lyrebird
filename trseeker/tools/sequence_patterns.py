#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (C) 2007-2009 Aleksey Komissarov ( ad3002@gmail.com )
# All rights reserved.
#
# This software is licensed as described in the file COPYING, which
# you should have received as part of this distribution.

"""
Function collection for work with sequence pattrens.

- __dna - dictionary for DNA15 <-> DNA5 alphabets translation
- PATTERNS - collection of common DNA patterns
- re_translate(sequence) -> DNA15
- re_get_plus_minus(sequence) -> revcom DNA15
- re_get_mutation(sequence) -> [sequence,...]
- get_double_pattern(pattern_static, pattern_dynamic) -> [patterns,...]
- *get_mutated_pattern(pattern, rate)* 
- get_mutated_pattern_twice(pattern, rate=2) -> [patterns, ...]
- get_mutated_pattern_trice(pattern, rate=3) -> [patterns, ...]
- remove_redundancy(list) -> [patterns,...]
- pattern_search(name, sequence, pattern_function) -> [(name, match start, match end),...]
"""
import re
from trseeker.tools.sequence_tools import get_revcomp
from trseeker.tools.other_tools import remove_redundancy

__dna = {'A':'[AN]',
         'C':'[CN]',
         'T':'[TN]',
         'G':'[GN]',
         'R':'[GAN]',
         'Y':'[TCN]',
         'M':'[ACN]',
         'K':'[GTN]',
         'S':'[GCN]',
         'W':'[ATN]',
         'H':'[ACTN]',
         'B':'[GCTN]',
         'V':'[GCAN]',
         'D':'[GATN]',
         'N':'[ACTGN]',
 }


DNA2IUPAC = {

     'A':'A',
     'C':'C',
     'T':'T',
     'G':'G',
     'GA':'R',
     'AG':'R',
     'TC':'Y',
     'CT':'Y',
     'AC':'M',
     'CA':'M',
     'GT':'K',
     'TG':'K',
     'GC':'S',
     'CG':'S',
     'AT':'W',
     'TA':'W',

     'ACT':'H',
     'CAT':'H',
     'CTA':'H',
     'ATC':'H',
     'TAC':'H',
     'TCA':'H',

     'GCT':'B',
     'GTC':'B',
     'TGC':'B',
     'TCG':'B',
     'CGT':'B',
     'CTG':'B',
     
     'GCA':'V',
     'GAC':'V',
     'AGC':'V',
     'ACG':'V',
     'CGA':'V',
     'CAG':'V',

     'GTA':'D',
     'GAT':'D',
     'AGT':'D',
     'ATG':'D',
     'TGA':'D',
     'TAG':'D',
     
     'ACGT': 'N',
     }    


PATTERNS = {'MARS1' : 'AATAAYAA',
            'MARS2' : 'AWWRTAANNWWGNNNC',
            'CENPB' : 'TTCGNNNNANNCGGG', # TTCG.{4}A.{2}CGGG   # CCCG.{2}T.{4}CGAA
            'PRDB9' : 'CCNCCNTNNCCNC',
            }

def re_translate(sequence):
    ''' Translate sequence in reg exp
    
    - sequence: sequence in DNA15
    '''
    return ''.join(__dna[nucleotide] for nucleotide in sequence)

def re_get_plus_minus(sequence):
    ''' Return +/- sequences in reg exp 
    
    - sequence: sequence in DNA15
    '''
    return [re_translate(sequence), get_revcomp(re_translate(sequence))]

def re_get_mutations(sequence):
    ''' Return list of sequences where one letter changed to N
    
    >>>get_mutation("ACNT")
    ['NCNT','ANNT','ACNN']
    
    '''
    # TODO: check out unique.
    out = set()
    seq = list(sequence)
    for i in range(0, len(seq)):
        temp = seq[:]
        if temp[i] == 'N':
            continue
        temp[i] = 'N'
        out.add("".join(temp))
    return list(out)

def get_double_pattern(pattern_static, pattern_dynamic):
    ''' Return list of patterns one not mutated and second mutated
    
    - pattern_static: pattern_static name
    - pattern_dynamic: pattern that mutated name
    '''
    patterns = set()
    patterns = patterns.union(re_get_plus_minus(pattern_static))
    patterns = patterns.union(re_get_plus_minus(pattern_dynamic))
    for variant in re_get_mutation(pattern_dynamic):
        patterns = patterns.union(re_get_plus_minus(variant))
    patterns.discard('')
    return list(patterns)

def get_mutated_pattern_twice(pattern):
    ''' Return list of twice mutated patterns
    
    - pattern: pattern name
    - rate: dummy parameter
    '''
    patterns = []
    patterns += re_get_plus_minus(pattern)
    for variant in re_get_mutation(pattern):
        for mutated_variant in re_get_mutation(variant):
            patterns += re_get_plus_minus(mutated_variant)
    return remove_redundancy(patterns)

def get_mutated_pattern_trice(pattern):
    ''' Return list of trice mutated patterns
    - pattern: pattern name
    - rate: dummy parameter
    '''
    patterns = []
    patterns += re_get_plus_minus(pattern)
    for variant in re_get_plus_minus(pattern):
        for mutated_variant in re_get_plus_minus(variant):
            patterns += re_get_plus_minus(mutated_variant)
            for more_mutated_variant in re_get_plus_minus(mutated_variant):
                patterns += re_get_plus_minus(more_mutated_variant)
    return remove_redundancy(patterns)

def get_motif_positions(sequence, motif):
    ''' Return list of motif positions in sequence.
    '''
    pos = 0
    result = []
    while True:
        r = sequence.find(motif, pos)
        if r == -1:
            return result
        result.append(r)
        pos = r+1


def pattern_search(name, sequence, pattern_function, pattern_function_params):
    ''' Return list triples (name, match start, match end) for given sequence
    that were found with given pattern function with given parameters
    
    - name: pattern name
    - sequence: sequence
    - pattern_function: pattern_function
    - pattern_function_params: parameters for pattern function
    '''
    result = set()
    m = False
    for reg_exp in pattern_function(*pattern_function_params):
        res = re.finditer(reg_exp, sequence, re.I)
        for match in res:
            result.add((name, match.start() + 1, match.end() + 1))
            m = True
        if m:
            break
    result = list(result)
    result.sort()
    return result
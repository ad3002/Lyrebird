#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (C) 2007-2009 Aleksey Komissarov ( ad3002@gmail.com )
# All rights reserved.
#
# This software is licensed as described in the file COPYING, which
# you should have received as part of this distribution.
"""
Function for simple sequence analysis (reduced part of Trseeker).

- get_revcomp(sequence)
- get_gc(sequence)
- check_gapped(sequence)
- clear_sequence(sequence)
- restriction(sequence, pattern, end="")
- restriction_fragments_n(sequence, pattern)
- get_subseq(seq, start, end)
- check_synonims(seq_a, seq_b)
- random_mutation(seq, n, alphabet="actgn ")
    
TODO: Not implemented:

- get_bin_sequence(sequence)
    
"""
import re
import random
from collections import defaultdict

def get_dust_score(sequence, k=4):
    ''' Compute dust score for givent score with given window size.
    '''
    d = defaultdict(int)
    for i in range(0, len(sequence)-k+1):
        d[sequence[i:i+k]] += 1
    score = 0.
    total = 0.
    for v in d.values():
        score += v*(v-1)/2.
        total += score
    return total/(len(sequence)-k+1)
    

def get_shifts_variants(sequence):
    '''
    '''
    shifts = set()
    for i in range(len(sequence)):
        shifts.add(sequence[i:]+sequence[:i])
    return list(shifts)
    

REVCOMP_DICTIONARY = dict(zip('ATCGNatcgn~[]', 'TAGCNtagcn~]['))

def get_revcomp(sequence):
    '''Return reverse complementary sequence.

    >>> complementary('AT CG')
    'CGAT'

    '''
    return ''.join(REVCOMP_DICTIONARY.get(nucleotide, '') for nucleotide in reversed(sequence))

def get_translation(sequence, from_string, to_string):
    '''Translate sequence from one alphabet to another.

    '''
    c = dict(zip(from_string, to_string))
    return ''.join(c.get(nucleotide, '') for nucleotide in sequence)

def get_comp(sequence):
    '''Return complementary sequence.

    >>> complementary('AT CG')
    'TAGC'

    '''
    c = dict(zip('ATCGNatcgn[]', 'TAGCNtagcn]['))
    return ''.join(c.get(nucleotide, '') for nucleotide in sequence)

def fix_strand(sequence):
    ''' Return normalized sequence with following rules:
    1) if T > A then return reverse complement
    2) if A == T and G > C  then return reverse complement
    3) else return sequence
    '''
    A = sequence.count('a')
    T = sequence.count('t')
    if T > A:
        return get_revcomp(sequence)
    if T == A:
        C = sequence.count('c')
        G = sequence.count('g')
        if G > C:
            return get_revcomp(sequence)
    return sequence

def get_gc(sequence):
    ''' Count GC content.
    '''
    length = len(sequence)
    if not length:
        return 0
    count_c = sequence.count('c')+sequence.count('C')
    count_g = sequence.count('g')+sequence.count('G')
    gc = float(count_c + count_g) / float(length)
    return float(gc)
    
def get_int_gc(sequence):
    ''' Get GC content from 0 to 100.
    '''
    gc = get_gc(sequence)
    return int(100*round(gc, 2))

def check_gapped(sequence):
    """ Check for N in sequence. Return n(with N) or w(whole).
    """
    w_regexp = re.compile('n|N')
    regexp_obj = w_regexp.search(sequence)
    if (regexp_obj):
        return True
    else:
        return False
    
def get_subseq(seq, start, end):
    ''' Return subsequence.
    '''
    return seq[start-1 : end]
    
def clear_sequence(sequence, lower=True):
    """ Clear sequence (full alphabet):
    
    - lower case
    - \s -> ""
    - [^actgn] -> ""
    """
    sequence = sequence.strip()
    if lower:
        sequence = sequence.lower().strip()
    sequence = re.sub("\s+", "", sequence)
    return re.sub("[^actgnuwsmkrybdhvACTGNUWSMKRYBDHV\-]", "", sequence)

def restriction(sequence, pattern, end=None):
    ''' Return list of sequences splited by pattern with added end.
    
    - patterm: reg_exp
    '''
    fragments = [ fragment for fragment in re.split(pattern, sequence) ]
    if len(fragments) == 1 and \
        fragments[0] == sequence:
            return fragments
    n = len(fragments)
    if not end:
        end = ""
    for i in range(0,n):
        if i == 0:
            fragments[i] += end
        elif i == n-1:
            fragments[i] = end + fragments[i]
        else:
            fragments[i] = end + fragments[i] + end
    return fragments
            
def restriction_fragments_n(sequence, pattern):
    '''Return number of fragments after restriction of sequence with given pattern.
    
    - patterm: reg_exp
    '''
    fragments = [ fragment for fragment in re.split(pattern, sequence) ]
    n = len(fragments)
    m = [fragment for fragment in fragments if fragment == ""]    
    return n - len(m)

def check_cyclic_repeats(seq_a, seq_b):
    ''' Check tandem repeat synonims  between two sequence with same length.

    >>> check_synonims("acc", "accacc")
    True 
    
    @param seq_a: sequence with length n
    '''
    if len(seq_a) != len(seq_b):
        print("Need sequences with same length")
        return False
    if seq_a in seq_b*2 or seq_b in seq_a*2:
        return True
    else:
        return False

def random_mutation(seq, n, alphabet="actgn +"):
    ''' Return sequence with n mutations.
    
    - alphabet: letters and ' ' - del, '+' - inserts
    '''
    N = len(seq)-1
    seq = [ x[1] for x in enumerate(seq) ]
    for mut in range(0, n):
        i = random.randint(0, N)
        j = random.randint(0, len(alphabet)-1)
        seq[i] = alphabet[j]
    mutated = "".join(seq)
    mutated.replace(" ","")
    mutated.replace("+",alphabet[random.randint(0, len(alphabet)-1)])
    return mutated

def get_consensus(strs):
    ''' Return consensus string for given list of strings.'''
    n = len(strs[0])
    for s in strs:
        assert len(s) == n
    result = {}
    result['A'] = [0]*n
    result['C'] = [0]*n
    result['T'] = [0]*n
    result['G'] = [0]*n
    for i in range(n):
        for s in strs:
            result[s[i]][i] += 1
    return result

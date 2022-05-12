#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 07.09.2010
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com
""" Edit distance functions.

- _delta(a,b) -> True/False
- get_pos(L,d) -> [positions]
- get_ed_similarity(s1,s2, monomer_mode=False, full_info=False, verbose=False) -> %sim
- get_edit_distance_info(s1,s2, verbose=False, monomer_mode=False) -> (r_str, pos, d, d*100/len1, n, len1, len2)
"""

import numpy
import re


def _delta(a, b):
    ''' Compare two elements, if one element equal "-", than return False.'''
    if a == "-" or b == "-":
        return False
    elif a == b:
        return True
    else:
        return False

def get_pos(L, d):
    ''' Return list of element (*d*) positions in given list.'''
    l = L[:]
    pos = []
    i = 0
    while d in l:
        pos.append(l.index(d) + i)
        l.remove(d)
        i += 1
    return pos

def get_ed_similarity(s1, s2, monomer_mode=False, full_info=False, verbose=False):
    ''' Return edit distance between two given strings.

    Edit distance dynamic programming implementation
    Length of second sequence does not change similarity value.
    
    - monomer mode: double short monomer sequence if len2/len1 < 2
    - full_info: boolean, return full information
    
    Return: percent of ED similarity, or if full_info (distance, number of d, positions of d, S matrix)
    '''
    len1 = len(s1)
    len2 = len(s2)
    if len1 == 0 or len2 == 0:
        return 0.0
    if len1 > len2:
        (s1, s2) = (s2, s1)
        (len1, len2) = (len2, len1)
    if monomer_mode:
        r = len2 / float(len1)
        if r < 2:
            s2 = s2 * 2
            len2 = len2 * 2
    S = numpy.zeros(shape=(len1 + 1, len2 + 1))
    for j in range(1, len2 + 1):
        S[0, j] = S[0, j - 1] + _delta("-", s2[j - 1])
    for i in range(1, len1 + 1):
        if verbose:
            print("%.2f" % (float(i) / len1))
        S[i, 0] = S[i - 1, 0] + _delta(s1[i - 1], "-")
        for j in range(1, len2 + 1):
            S[i, j] = max(S[i - 1, j - 1] + _delta(s1[i - 1], s2[j - 1]),
                          S[i - 1, j] + _delta(s1[i - 1], "-"),
                          S[i - 1, j - 1] + _delta("-", s2[j - 1])
                          )
    result = [ S[len1, j] for j in range(1, len2 + 1) ]
    d = max(result)
    if full_info:
        n = result.count(d)
        pos = get_pos(result, d)
        for i in range(0, len(pos)):
            pos[i] = pos[i] + 1 - len1
        return (d, n, pos, S)
    return float(d) * 100 / len1

def get_edit_distance_info(s1, s2, verbose=False, monomer_mode=False):
    ''' Return two sequences edit distance full information.
    
    - s1: sequence of tandem repeat monomer
    - s2: sequence of tandem repeat monomer
    - verbose: boolean
    - monomer mode: double short monomer sequence if len2/len1 < 2
    
    Return: ("max_sim Nmax_sim %sim pos_list", pos_list, all_result, length seq1, length seq2)
    '''
    len1 = len(s1)
    len2 = len(s2)
    if len1 == 0 or len2 == 0:
        return ("", [], 0, 0.0, 0, 0, 0)
    if len1 > len2:
        (s1, s2) = (s2, s1)
        (len1, len2) = (len2, len1)
    (d, n, pos, S) = get_ed_similarity(s1, s2, full_info=True, verbose=verbose)
    # Create output string
    # (max_sim, N max_sim, %sim, pos_list)
    r_str = "%s\t%s\t%s\t%s\t" % (d, n, d * 100 / len1, str(pos))
    # Return ("max_sim\tNmax_sim\t%sim\tpos_list", pos_list, sim, %sim, n, length seq1, length seq2)
    return (r_str, pos, d, d * 100 / len1, n, len1, len2)

def get_edit_distance(s1, s2):
    ''' Return ED valie for two sequences.''' 
    return int(get_ed_similarity(s1, s2, monomer_mode=False, full_info=True)[0])

def get_edit_distance_row(s1, s2):
    ''' Get last row of ED matrix between two strings.'''
    return get_ed_similarity(s1, s2, monomer_mode=False, full_info=True)[-1][-1]

def hamming_distance(s1, s2):
    """ Get Hamming distance: the number of corresponding symbols that differs in given strings.
    """
    return sum(i != j for i, j in zip(s1, s2))

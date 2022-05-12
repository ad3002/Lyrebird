#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 07.09.2011
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com 
'''
Functions related to TRs group analysis.
'''
from trseeker.tools.statistics import get_mean

def get_index(i):
    ''' Get next family index.'''
    s = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    if i<len(s):
        return s[i]
    else:
        return "X%s" % i

def get_popular(s):
    ''' Get most frequent element of list.'''
    r = {}
    for x in s:
        r.setdefault(x, 0)
        r[x] += 1
    c = [(v, k) for k, v in r.items()]
    c.sort(reverse=True)
    return c[0][1]

def get_family_name(trf_objs, seen_units):
    ''' Get unit and letter for family.'''
    
    # read units to pmatch
    units = {}
    for trf_obj in trf_objs:
        units.setdefault(trf_obj.trf_period, [])
        units[trf_obj.trf_period].append(trf_obj.trf_pmatch)

    # find optimal unit
    max_pmatch = 0
    for u in units:
        m = get_mean(units[u]) * len(units[u])
        if m > max_pmatch:
            max_pmatch = m
            right_u = u
        if m == max_pmatch:
            right_u = min(u, right_u)

    # change unit
    if not right_u in seen_units:
        seen_units[right_u] = 0
        letter = get_index(seen_units[right_u])
    else:
        seen_units[right_u] += 1
        letter = get_index(seen_units[right_u])
    return right_u, letter, seen_units

def join_families_with_common(families):
    ''' Join families with common members.'''
    n = len(families)
    for i in xrange(0, n):
        if not families[i]:
            continue
        for j in xrange(i + 1, n):
            if not families[j]:
                continue
            if families[i].intersection(families[j]):
                families[i] = families[i].union(families[j])
                families[j] = None
    families = [x for x in families if x]
    for i in xrange(len(families)):
        families[i] = list(families[i])
    families.sort(key=lambda x: len(x), reverse=True)
    return families
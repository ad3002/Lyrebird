#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 16.02.2011
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com 
'''
Other usefull functions.
'''
import time
import math

def sort_dictionary_by_key(d, reverse=False):
    ''' Sort dictionary by key. Retrun list of (k, v) pairs.'''
    result = [(k, v) for k, v in d.items()]
    result.sort(reverse=reverse)
    return result

def as_strings(alist):
    ''' Cast list elements to string.'''
    return [str(x) for x in alist]

def dict2list(adict):
    ''' Return list from dictionary, return list of (key, value) pairs.'''
    return [(k, v) for k, v in adict.items()]

def sort_dictionary_by_value(d, reverse=False):
    ''' Sort dictionary by value. Retrun list of (v, k) pairs.'''
    result = [(v, k) for k, v in d.items()]
    result.sort(reverse=reverse, key=lambda x: x[0])
    return result

def remove_duplicates_and_sort(data):
    ''' Remove duplicates from list and sort it.'''
    data = list(set(data))
    data.sort()
    return data

def remove_duplicates(data):
    ''' Remove duplicates from list.'''
    return list(set(data))

def remove_redundancy(alist):
    ''' Remove redundancy or empty elements in given list.
    '''
    s = set()
    for item in alist:
        s.add(item)
    s.discard('')
    return list(s)

def clear_fragments_redundancy(data, extend=False, same_case_func=None):
    ''' Remove nested fragments, ata format [start, end, ...].
    '''
    data.sort()

    last = None
    last_start = None
    last_end = None
    extended = False
    for i in xrange(0, len(data)):

        if not last:
            last = i
            last_start = data[i][0]
            last_end = data[i][1]
            continue

        start = data[i][0]
        end = data[i][1]

        # -----
        # ----- x

        if start == last_start and end == last_end:
            if same_case_func and hasattr(same_case_func, '__call__'):
                if same_case_func(data[last], data[i]):
                    data[i] = None
                else:
                    data[last] = None
            else:
                data[i] = None
            continue

        # -----
        #  ---  x 
        # ---   x
        #   --- x
        if last_start <= start and last_end > end:
            data[i] = None
            continue
        if last_start < start and last_end >= end:
            data[i] = None
            continue

        # -----   
        # -------
        if  end > last_end and start == last_start:
            last_end = end
            last_start = start
            data[last] = None
            last = i
            continue

        # -----
        #    -----
        if last_start < start < last_end and last_end < end:
            if extend:
                extended = True
                last_end = end
                data[last] = None
                last = i
                continue
            else:
                yield data[last], False
                last_end = end
                last_start = start
                last = i
                continue

        yield data[last], extended

        last_end = end
        last_start = start
        last = i

        extended = False
    yield data[last], extended
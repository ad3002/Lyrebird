#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (C) 2007-2009 Aleksey Komissarov ( ad3002@gmail.com )
# All rights reserved.
#
# This software is licensed as described in the file COPYING, which
# you should have received as part of this distribution.
'''
Varios usefull functions. Again.

TODO: join with other_tools
'''
def save_list(file, l):
    ''' A function saves python list *l* in *file*.
    '''
    with open(file, "ab") as fw:
        [ fw.write("%s\n" % x) for x in l ]

def save_dict(file, d):
    ''' A function saves dictionatu *d* in *file*.
    '''
    with open(file, "ab") as fw:
        [ fw.write("%s\t%s\n" % (x, d[x])) for x in d.keys() ]

def save_sorted_dict(d, file, by_value=True, reverse=True, min_value=None, key_type=None):
    ''' A function saves dictionatu *d* in *file*.
    '''
    if min_value is None:
        d = [(k, v) for k, v in d.items()]
    else:
        d = [(k, v) for k, v in d.items() if v>min_value]
    if key_type:
        d = [(key_type(k), v) for (k, v) in d]
    if by_value:
        d.sort(key=lambda x: x[1], reverse=reverse)
    else:
        d.sort(reverse=reverse)

    with open(file, "w") as fw:
        [fw.write("%s\t%s\n" % (x, y)) for x, y in d]

def count_lines(file):
    ''' Return number of line in file
    '''
    n = 0
    with open(file, "rb") as fh:
        for line in fh:
            if line:
                n += 1
    return n

def sort_file_by_int_field(file_name, field):
    ''' Sort tab-delimited file data by give int field.'''
    data = []
    with open(file_name, "r") as fh:
        data = fh.readlines()
        data = [x.split("\t") for x in data]

    data.sort(reverse=True, key=lambda x: int(x[field]))

    with open(file_name, "w") as fh:
        data = ["\t".join(x) for x in data]
        fh.writelines(data)
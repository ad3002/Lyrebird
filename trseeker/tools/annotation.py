#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 07.09.2010
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com
"""
Functions related to annotation data.
"""

def compute_intersection_intervals(data_a, data_b):
    """ Compute intersection between two datasets.
    Intersection types: match, nested, overlap, and adjacent
    @param data_a: list of tuples (chrName, startPos, endPos, feature_id)
    @param data_b: list of tuples (chrName, startPos, endPos, feature_id)
    @return: (a_chr, a_start, a_end, a_feature, b_chr, b_start, b_end, b_feature, intersection_type) tuple
    """
    data_a.sort()
    data_b.sort()
    if not data_a or not data_b:
        return []
    if data_a[0][0] > data_b[0][0]:
        data_a, data_b = data_b, data_a
    (a_chr, a_start, a_end, a_feature) = data_a.pop(0)
    (b_chr, b_start, b_end, b_feature) = data_b.pop(0)
    hits = []
    while True:
        if b_chr != a_chr:
            if not data_b:
                break
            (b_chr, b_start, b_end, b_feature) = data_b.pop(0)
            continue
        # -----
        # -----
        if a_start == b_start and a_end == b_end:
            hits.append((a_chr, a_start, a_end, a_feature, b_chr, b_start, b_end, b_feature, "match"))
            if not data_a or not data_b:
                break
            (b_chr, b_start, b_end, b_feature) = data_b.pop(0)
            (a_chr, a_start, a_end, a_feature) = data_a.pop(0)
            continue
        # -----
        #  ---
        if a_start <= b_start and a_end >= b_end:
            hits.append((a_chr, a_start, a_end, a_feature, b_chr, b_start, b_end, b_feature, "nested"))
            if not data_b:
                break
            (b_chr, b_start, b_end, b_feature) = data_b.pop(0)
            continue
        #  ---
        # -----
        if a_start >= b_start and a_end <= b_end:
            hits.append((a_chr, a_start, a_end, a_feature, b_chr, b_start, b_end, b_feature, "nested"))
            if not data_a:
                break
            (a_chr, a_start, a_end, a_feature) = data_a.pop(0)
            continue
        #     -----
        # -----
        if a_start > b_start and a_end > b_end and b_end > a_start:
            hits.append((a_chr, a_start, a_end, a_feature, b_chr, b_start, b_end, b_feature, "overlap"))
            if not data_b:
                break
            (b_chr, b_start, b_end, b_feature) = data_b.pop(0)
            continue
        # -----
        #     -----
        if a_start < b_start and a_end < b_end and a_end > b_start:
            hits.append((a_chr, a_start, a_end, a_feature, b_chr, b_start, b_end, b_feature, "overlap"))
            if not data_a:
                break
            (a_chr, a_start, a_end, a_feature) = data_a.pop(0)
            continue
        # -----
        #         -----
        if a_end == b_start:
            hits.append((a_chr, a_start, a_end, a_feature, b_chr, b_start, b_end, b_feature, "adjacent"))
            if not data_a:
                break
            (a_chr, a_start, a_end, a_feature) = data_a.pop(0)
            continue
        #         -----
        # -----
        if b_end == a_start:
            hits.append((a_chr, a_start, a_end, a_feature, b_chr, b_start, b_end, b_feature, "adjacent"))
            if not data_b:
                break
            (b_chr, b_start, b_end, b_feature) = data_b.pop(0)
            continue
        #           -----
        # -----
        if b_end < a_start:
            if not data_b:
                break
            (b_chr, b_start, b_end, b_feature) = data_b.pop(0)
            continue
        # -----
        #            -----
        if a_end < b_start:
            if not data_a:
                break
            (a_chr, a_start, a_end, a_feature) = data_a.pop(0)
            continue
    return hits


def test_compute_intersection_intervals():
    """ Test for compute_intersection_intervals function.
    @todo: move it to testing framework
    @return: None
    """
    data_a = [
        ("chrA", 0, 10, 0),
        ("chrA", 25, 40, 1),
        ("chrA", 45, 50, 2),
        ("chrA", 60, 70, 3),
        ("chrA", 80, 90, 4),
        ("chrA", 100, 150, 5),
        ("chrA", 250, 270, 6),
        ("chrA", 290, 340, 7),
    ]

    data_b = [
        ("chrA", 25, 40, 10),
        ("chrA", 47, 49, 20),
        ("chrA", 59, 90, 30),
        ("chrA", 80, 90, 40),
        ("chrA", 125, 175, 50),
        ("chrA", 275, 280, 60),
    ]
    data = compute_intersection_intervals(data_a, data_b)
    assert data[0] == ('chrA', 25, 40, 1, 'chrA', 25, 40, 10, 'match')
    assert data[1] == ('chrA', 45, 50, 2, 'chrA', 47, 49, 20, 'nested')
    assert data[2] == ('chrA', 60, 70, 3, 'chrA', 59, 90, 30, 'nested')
    assert data[3] == ('chrA', 80, 90, 4, 'chrA', 59, 90, 30, 'nested')
    assert data[4] == ('chrA', 100, 150, 5, 'chrA', 125, 175, 50, 'overlap')
    assert [] == compute_intersection_intervals([], [])
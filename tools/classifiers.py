#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 21.05.2014
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com

from trseeker.tools.sequence_tools import get_revcomp
from trseeker.tools.sequence_tools import get_shifts_variants


def ask_name(summary=None, default_name=None, default_family=None):
    """ Ask and return family name and inner name.
    """

    if summary:
        for x in summary:
            if x[0] == "Unknown":
                continue
            default_name = x[0]
            default_family = x[1]
            break

    inner_name = raw_input("Please provide inner name (%s):" % default_name).strip() or default_name
    if not inner_name:
        print("Skipped")
        return None
    inner_family = raw_input("Please provide inner family (%s):" % default_family).strip() or default_family
    if not inner_family:
        print("Skipped")
        return None
    return inner_name, inner_family


def ask_trs_name(sequence, length, prefix, tr_family_to_letter, fast=False):
    """ Ask and return family name and inner name for tandem repeat.
    """

    monomer = sequence[:length]
    variants = set(get_shifts_variants(monomer) + get_shifts_variants(get_revcomp(monomer)))
    monomer = min(variants)


    if length in tr_family_to_letter:
        if tr_family_to_letter[length]+65 < 91:
            letter = chr(tr_family_to_letter[length]+65)
        else:
            letter = ":%s" % tr_family_to_letter[length]
        tr_family_to_letter[length] += 1
    else:
        tr_family_to_letter[length] = 1
        letter = 'A'


    if length < 10:
        default_name = "(%s)n" % monomer
        default_family = "Microsatellite"
    else:
        default_name = "%s%s%s" % (prefix, length, letter)
        default_family = "SatDNA"

    if fast:
        return monomer, default_name, default_family

    inner_name = raw_input("Please check name (%s):" % default_name).strip() or default_name
    if not inner_name:
        print("Skipped")
        return monomer, None, None
    inner_family = raw_input("Please check family (%s):" % default_family).strip() or default_family
    if not inner_family:
        print("Skipped")
        return monomer, None, None
    return monomer, inner_name, inner_family

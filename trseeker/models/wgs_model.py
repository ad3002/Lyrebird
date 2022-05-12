#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 25.12.2010
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com 

from PyExp import AbstractModel

class WGSModel(AbstractModel):
    '''Model for reading WGS meta data.

    Dumpable attributes:

    - "wgs_prefix",
    - "wgs_taxon",
    - "wgs_gpid",
    - "wgs_acronym",
    - "wgs_contigs" (int),
    - "wgs_length" (int),
    - "wgs_mean_gc" (float),
    - "wgs_trs_all (int)",
    - "wgs_trs_3000" (int),
    - "wgs_trs_1000" (int),
    - "wgs_trs_500" (int),
    - "wgs_trs_all_proc" (float),
    - "wgs_trs_3000_proc" (float),
    - "wgs_trs_1000_proc" (float),
    - "wgs_trs_500_proc" (float),
    - "wgs_trs_all_length" (int),
    - "wgs_trs_3000_length" (int),
    - "wgs_trs_1000_length" (int),
    - "wgs_trs_500_length" (int),
    - "wgs_sum_gc" (float)

    Usage:

    >>> wgs_model = WGSModel(prefix, taxon, gpid, acronym=None)

    '''

    dumpable_attributes = ["wgs_prefix",
                           "wgs_taxon",
                           "wgs_gpid",
                           "wgs_acronym",
                           "wgs_contigs",
                           "wgs_length",
                           "wgs_mean_gc",
                           "wgs_trs_all",
                           "wgs_trs_3000",
                           "wgs_trs_1000",
                           "wgs_trs_500",
                           "wgs_trs_all_proc",
                           "wgs_trs_3000_proc",
                           "wgs_trs_1000_proc",
                           "wgs_trs_500_proc",
                           "wgs_trs_all_length",
                           "wgs_trs_3000_length",
                           "wgs_trs_1000_length",
                           "wgs_trs_500_length",
                           "wgs_sum_gc",
                           ]
    int_attributes = ["wgs_contigs",
                      "wgs_length",
                      "wgs_trs_all",
                      "wgs_trs_3000",
                      "wgs_trs_1000",
                      "wgs_trs_500",
                      "wgs_trs_all_length",
                      "wgs_trs_3000_length",
                      "wgs_trs_1000_length",
                      "wgs_trs_500_length",
                      ]
    float_attributes = ["wgs_mean_gc",
                        "wgs_trs_all_proc",
                        "wgs_trs_3000_proc",
                        "wgs_trs_1000_proc",
                        "wgs_trs_500_proc",
                        "wgs_sum_gc",

                        ]

    def __init__(self, prefix, taxon, gpid, acronym=None):
        
        super(WGSModel, self).__init__()
        self.wgs_prefix = prefix
        self.wgs_taxon = taxon
        self.wgs_gpid = gpid
        self.wgs_acronym = acronym
        
        
    def preprocess_data(self):

        if self.wgs_trs_all_length:
            self.wgs_trs_all_proc = self.wgs_trs_all_length / float(self.wgs_length)
        if self.wgs_trs_3000_length:
            self.wgs_trs_3000_proc = self.wgs_trs_3000_length / float(self.wgs_length)
        if self.wgs_trs_1000_length:
            self.wgs_trs_1000_proc = self.wgs_trs_1000_length / float(self.wgs_length)
        if self.wgs_trs_500_length:
            self.wgs_trs_500_proc = self.wgs_trs_500_length / float(self.wgs_length)
        if not self.wgs_mean_gc:
            self.wgs_mean_gc = self.wgs_sum_gc / self.wgs_contigs

    def clear_trf(self):
        ''' Clear trf information.'''
        self.wgs_trs_all = 0
        self.wgs_trs_all_length = 0
        self.wgs_trs_3000 = 0
        self.wgs_trs_3000_length = 0
        self.wgs_trs_1000 = 0
        self.wgs_trs_1000_length = 0
        self.wgs_trs_500 = 0
        self.wgs_trs_500_length = 0



#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 28.11.2010
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com 

from PyExp import AbstractModel

class DatasetModel(AbstractModel):
    """ Class for dataset wrapping.

    Dumpable attributes:

    - dataset_taxon
    - dataset_gi
    - dataset_sources
    - dataset_description
    - dataset_gc (float)
    - dataset_length (int)
    - dataset_trs_n (int)
    - dataset_length (int)
    - dataset_trs_mean_gc (float)
    - dataset_trs_fraq (float)

    """

    dumpable_attributes = ["dataset_taxon",
                           "dataset_gi",
                           "dataset_sources",
                           "dataset_description",
                           "dataset_gc",
                           "dataset_length",
                           "dataset_trs_n",
                           "dataset_trs_length",
                           "dataset_trs_mean_gc",
                           "dataset_trs_fraq",
                           ]
                           
    int_attributes = ["dataset_length",
                       "dataset_trs_n",
                       "dataset_trs_length"]

    float_attributes = ["dataset_gc",
                        "dataset_trs_mean_gc",
                       "dataset_trs_fraq", 
                       ]
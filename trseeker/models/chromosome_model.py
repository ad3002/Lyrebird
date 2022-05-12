#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 08.09.2011
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com 
from PyExp import AbstractModel

class ChomosomeModel(AbstractModel):
    ''' Chromosome model.

    Dumpable attributes:

    - "chr_genome",
    - "chr_number",
    - "chr_taxon",
    - "chr_prefix",
    - "chr_gpid",
    - "chr_acronym",
    - "chr_contigs",
    - "chr_length",
    - "chr_mean_gc",
    - "chr_trs_all",
    - "chr_trs_3000",
    - "chr_trs_all_proc",
    - "chr_trs_3000_proc",
    - "chr_trs_all_length",
    - "chr_trs_3000_length",
    - "genome_gaps",
    - "chr_sum_gc",
    
    '''


    dumpable_attributes = [
                    "chr_genome",
                    "chr_number",
                    "chr_taxon",
                    "chr_prefix",
                    "chr_gpid",
                    "chr_acronym",
                    "chr_contigs",
                    "chr_length",
                    "chr_mean_gc",
                    "chr_trs_all",
                    "chr_trs_3000",
                    "chr_trs_all_proc",
                    "chr_trs_3000_proc",
                    "chr_trs_all_length",
                    "chr_trs_3000_length",
                    "genome_gaps",
                    "chr_sum_gc",
                    ]

    def preprocess_data(self):
        if self.chr_trs_all_length:
            self.chr_trs_all_proc = self.chr_trs_all_length / float(self.chr_length)
        if self.chr_trs_3000_length:
            self.chr_trs_3000_proc = self.chr_trs_3000_length / float(self.chr_length)
        if not self.chr_mean_gc:
            self.chr_mean_gc = self.chr_sum_gc / self.chr_contigs

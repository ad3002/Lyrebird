#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 
#@author:
#@contact:

from PyExp import AbstractModel


class GenbankData(AbstractModel):
    
    dumpable_attributes = [
        "gb_locus",
        "gb_moltype",
        "gb_divison",
        "gb_moltype",
        "gb_moddata",
        "gb_length",
        "gb_definition",
        "gb_accessions",
        "gb_versions",
        "gb_keywords",
        "gb_source",
        "gb_organism",
        "gb_dblink",
        "gb_references",
        "gb_comments",
        "gb_features",
        "gb_contig",
        
        
    ]
    
    list_attributes = [
        "gb_accessions",
        "gb_versions",
        "gb_keywords",
        "gb_dblink",
        "gb_references",
        "gb_comments",
        "gb_features",
    ]
    
    int_attributes = [
        "gb_length",
    ]

    
    def __init__(self):
        super(GenbankData, self).__init__()
        
        self.gb_dblink = []
        self.gb_references = []
        self.gb_comments = []
        self.gb_features = []
        
    def parse_locus(self, line):
        d = [x for x in line.split() if x]
         # ['LOCUS', 'NZ_JSZC01000016', '295444', 'bp', 'DNA', 'linear', 'CON', '23-DEC-2019']
        if len(d) == 7:
            _, self.gb_locus, self.gb_length, _, self.gb_moltype, self.gb_divison, self.gb_moddata = d
        else:
            _, self.gb_locus, self.gb_length, _, self.gb_moltype, gb_moltype, self.gb_divison, self.gb_moddata = d
        
    
#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 20.06.2011
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com 

from PyExp import AbstractModel

class NgramModel(AbstractModel):
    ''' Ngram model.

    TODO: check it.

    Dumpable attributes:

    - seq_r
    - seq_f
    - tf (int)
    - df (int)
    - taxons
    - trs
    - families

    Usage:

    >>> ng = NgramModel(seq_f, seq_r)
    >>> for trf_obj in tr_lisT:
    ...     ng.add(trf_obj, tf)

    
    String representation: **fseq rseq tf df proj fams**
            
    - forward seq
    - reverse seq
    - global tf
    - global df
    - projects
    - families

    '''

    dumpable_attributes = [
                           "seq_f", 
                           "seq_r", 
                           "tf", 
                           "df", 
                           "taxons", 
                           "trs", 
                           "families"
                           ]

    int_attributes = ["df"]
    float_attributes = ["tf"]

    def __init__(self, seq_f, seq_r):

        super(NgramModel, self).__init__()

        self.seq_f = seq_f
        self.seq_r = seq_r
        self.taxons = set()
        self.trs = set()
        self.families = {}

    def add_tr(self, trf_obj, tf):
        ''' Add tandem repeat to model.'''

        self.tf += tf
        self.df += 1

        id = "%s_%s" % (trf_obj.project, trf_obj.trf_id)

        self.taxons.add(trf_obj.project)
        self.trs.add(id)

        family = "%s_%s" % (trf_obj.project, trf_obj.trf_family)

        self.families.setdefault(family, {"tf":0, "df":0, "ids":set()})
        self.families[family]["tf"] += tf
        self.families[family]["df"] += 1
        self.families[family]["ids"].add(id)

    def __str__(self):
        """ String representation.
            fseq\trseq\tf\tdf\tproj\tfams\n
            
            forward seq
            reverse seq
            global tf
            global df
            projects
            families
        """

        data = [self.seq_f,
                self.seq_r,
                self.tf,
                self.df,
                len(self.taxons),
                len(self.families),
                ]

        return "%s\n" % "\t".join([str(x) for x in data])

    def get_families(self):
        ''' Get family data.

        Return:

        ???
        '''
        result = []
        for family in self.families:
            if self.families[family]["tf"] < 10 or len(self.taxons) < 2:
                continue
            data = "%s,%s,%s" % (family, self.families[family]["tf"], self.families[family]["df"])
            result.append(data)
        fresult = ";   ".join(result)
        return len(result), "%s\t%s\n" % (self.seq_f, fresult)
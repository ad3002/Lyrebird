#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 07.09.2011
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com 

from PyExp import AbstractModel
from trseeker.seqio.tab_file import sc_iter_tab_file

class GFF3Model(AbstractModel):
    ''' Class for gff3 data.
    '''

    dumpable_attributes = ["seqid",
                           "source",
                           "ftype",
                           "start",
                           "end",
                           "score",
                           "strand",
                           "phase",
                           "attributes",
                           ]

    int_attributes = [ "start",
                       "end",
                       ]

    float_attributes = [ "score",
                       ]

    def set_with_bed_obj(self, bed_obj, **kwargs):
      ''' Set object with BED object.
      '''
      self.seqid = bed_obj.seqid
      self.start = bed_obj.start
      self.end = bed_obj.end
      if not "source" in kwargs:
        kwargs["source"] = "."
      self.source = kwargs["source"]
      if not "ftype" in kwargs:
        kwargs["ftype"] = "."
      self.ftype = kwargs["ftype"]
      if not "score" in kwargs:
        kwargs["score"] = "."
      self.score = kwargs["score"]
      if not "strand" in kwargs:
        kwargs["strand"] = "."
      self.strand = kwargs["strand"]
      if not "phase" in kwargs:
        kwargs["phase"] = "."
      self.phase = kwargs["phase"]
      if not "attributes" in kwargs:
        kwargs["attributes"] = "."
      self.attributes = kwargs["attributes"]


class BEDModel(AbstractModel):
    ''' Class for BED data.
    '''

    dumpable_attributes = ["seqid",
                           "start",
                           "end",
                           ]

    int_attributes = [ "start",
                       "end",
                       ]

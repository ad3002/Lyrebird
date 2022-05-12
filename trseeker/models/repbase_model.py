#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 20.08.2014
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com 

from PyExp import AbstractModel

class RepbaseModel(AbstractModel):
    """ Class for repbase data.
    """

    dumpable_attributes = [
      "score",
      "pdivergence",
      "pdeletions",
      "pinsertions",
      "query",
      "qstart",
      "qend",
      "qrem",
      "strand",
      "name",
      "family",
      "rrem",
      "rstart",
      "rend",
    ]
                           
    int_attributes = [
      "score",
      "qstart",
      "qend",
      "qrem",
      "rrem",
      "rstart",
      "rend",
    ]

    float_attributes = [
      "pdivergence",
      "pdeletions",
      "pinsertions",
    ]

    def set_element_coverage(self):
      if self.strand == "C":
        self.rlength = self.rstart + self.rrem
      else:
        self.rlength = self.rstart + self.rend
      self.pcov = round(float(abs(self.qstart-self.qend))/self.rlength, 2)

    def get_sequence(self, genome):
      '''
      '''
      if not self.query in genome:
        return None
      sequence = genome[self.query][self.qstart:self.qend]
      return sequence


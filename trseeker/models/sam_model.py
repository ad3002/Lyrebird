#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 28.11.2010
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com 

from PyExp import AbstractModel

class SAMModel(AbstractModel):
    """ Class for SAM data wrapping.

    Dumpable attributes:

    - QNAME
    - FLAG
    - RNAME
    - POS
    - MAPQ
    - CIGAR
    - RNEXT
    - PNEXT
    - TLEN
    - SEQ
    - QUAL
    - features
    
    Example:
    FCC2H05ACXX:5:2308:18685:69859#/2       4       *       0       0       *       *       0       0       CCCTGCTTATGCTCTGTCTCTCTCTGTCTCAAAAATAAATAAAAACATTAAAGAAATTAAAAAAAAAAAAAGAAAAA  CCCFFFFFHGHHHIIJIJJJJIJIJJGGGIJIJJJJJJJJJIJJHGHGJIHGHHIIIJJJJJIJIHFDD==&8ACCD   YT:Z:UU


    """

    dumpable_attributes = [
        "QNAME",
        "FLAG",
        "RNAME",
        "POS",
        "MAPQ",
        "CIGAR",
        "RNEXT",
        "PNEXT",
        "TLEN",
        "SEQ",
        "QUAL",
        "features",

    ]
    
    int_attributes = [
        "FLAG",
        "POS",
        "MAPQ",
        "PNEXT",
        "TLEN",
    ]

    @property
    def fragment_length(self):
        return self.POS - self.PNEXT

    def save_original(self, line):
        self.original = line

    def as_sam(self):
        '''
        '''
        s = []
        for attr in self.dumpable_attributes:
            if not attr == "features":
                s.append(str(getattr(self, attr)))
        s.append(self.raw_features)
        s = "\t".join(s)
        return "%s\n" % s

    def as_fastq(self):
        return "@%s\n%s\n+\n%s\n" % (self.QNAME, self.SEQ, self.QUAL)


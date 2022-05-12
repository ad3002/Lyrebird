#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 20.11.2014
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com 

from PyExp import AbstractModel

class RepeatModel(AbstractModel):
    """ Class for repeat data.
    """

    dumpable_attributes = [
      "rid",
      "rsid",
      "consensus_sequence",
      "consensus_length",
    ]
                           
    int_attributes = [
      "consensus_length",
    ]

    float_attributes = [
    ]

    

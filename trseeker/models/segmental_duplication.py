#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 13.03.2014
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com

import os, sys


from PyExp import AbstractModel

class LastzModel(AbstractModel):
    ''' Container for information about alignments.

    score,length1,length2,identity,name1,strand1,zstart1,end1,name2,strand2,zstart2,end2
    5501	65	65	60/65	92.3%	2L	+	142938	143003	2L	+	142940	143005
    '''

    dumpable_attributes = [
            "score",
            "length1",
            "length2",
            "identity",
            "identityp",
            "name1",
            "strand1",
            "zstart1",
            "end1",
            "name2",
            "strand2",
            "zstart2",
            "end2",
    ]

    int_attributes = [
        "score",
        "length1",
        "length2",
        "zstart1",
        "end1",
        "zstart2",
        "end2",
    ]

    float_attributes = [
    	"identityp",
    ]

    def preprocess_pair(self, key, value):
    	if key == "identityp":
    		value = value[:-1]
    	return key, value

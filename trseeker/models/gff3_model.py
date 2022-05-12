#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 28.10.2014
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com 
"""
"""


from PyExp import AbstractModel
from trseeker.seqio.tab_file import TabDelimitedFileIO
import csv
import collections
from PyExp import trseeker_logger


class Gff3Model(AbstractModel):
    """ Class for gff3 data wrapping.
    """

    dumpable_attributes = [
        "seqid",
        "source",
        "type",
        "start",
        "end",
        "score",
        "strand",
        "phase",
        "attributes",
    ]
    
    int_attributes = [
        "start",
        "end",
    ]

    @property
    def length(self):
        return abs(self.end - self.start)

    def get_coordinates(self):
        if self.end < self.start:
            return (self.end, self.start)
        return (self.start, self.end)

    def save_original(self, line):
        self.original = line

    def as_gff3(self):
        '''
        '''
        s = []
        for attr in self.dumpable_attributes:
            if not attr == "attributes":
                s.append(str(getattr(self, attr)))
        if not hasattr(self, "attributes"):
            self.attributes = {}     
        if self.attributes:
            attr_keys = self.attributes.keys()
            attr_keys.sort()
            attr = []
            for k in attr_keys:
              attr.append("%s=%s" % (k, self.attributes[k]))  
            self.raw_features = ";".join(attr)
        s.append(self.raw_features)
        s = "\t".join(s)
        return "%s\n" % s


class Gff3FeatureDict(collections.MutableMapping):
    """A dictionary for gff3 features
    """

    def __init__(self, *args, **kwargs):
        self.store = dict()
        self.update(dict(*args, **kwargs))  # use the free update to set keys

    def __getitem__(self, key):
        return self.store[self.__keytransform__(key)]

    def __setitem__(self, key, value):
        self.store[self.__keytransform__(key)] = value

    def __delitem__(self, key):
        del self.store[self.__keytransform__(key)]

    def __iter__(self):
        return iter(self.store)

    def __len__(self):
        return len(self.store)

    def __keytransform__(self, key):
        return key

class Gff3FileIO(TabDelimitedFileIO):
    """ 
    """

    def __init__(self, *args, **kwargs):
        """ 
        """
        super(TabDelimitedFileIO, self).__init__(*args, **kwargs)
        self.headers = []

    def read_online(self, file_name):
        """ Overrided. Yield items online from data from input_file."""

        def skip_comments(iterable):
            for line in iterable:
                if not line.startswith('#'):
                    yield line

        with open(file_name) as fh:
            for i, line in enumerate(fh):
                if line.startswith("#"):
                    self.headers.append(line)
                else:
                    break

        fields = Gff3Model().dumpable_attributes

        with open(file_name) as fh:
            for data in csv.DictReader(skip_comments(fh), fieldnames=fields, delimiter='\t', quoting=csv.QUOTE_NONE):
                _features = {}
                if data["attributes"]:
                    data["raw_features"] = data["attributes"]
                    for item in data["attributes"].split(";"):
                        if not item.strip():
                            continue
                        k, v = item.strip().split("=")
                        _features[k] = v
                data["attributes"] = _features
                obj = Gff3Model()
                try:
                    obj.set_with_dict(data)
                except:
                    trseeker_logger.error("Can't parse features for %s" % str(data))
                    continue
                yield obj


def sc_gff3_reader(gff3_file):
    """
    """
    reader = Gff3FileIO()
    for gff3_obj in reader.read_online(gff3_file):
        yield gff3_obj


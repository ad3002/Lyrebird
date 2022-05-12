#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 05.06.2011
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com 

"""
Classes:

- TRFFileIO(AbstractBlockFileIO)

"""
import re, os
from trseeker.models.sam_model import SAMModel
from trseeker.seqio.tab_file import TabDelimitedFileIO
import csv
import collections
from PyExp import runner

class SAMFeatureDict(collections.MutableMapping):
    """A dictionary for SAM features
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

class SAMFileIO(TabDelimitedFileIO):
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
                if line.startswith("    ") or line.startswith(" ") or line.startswith("perl"):
                    continue
                if not line.startswith('@'):
                    yield line

        with open(file_name) as fh:
            for i, line in enumerate(fh):
                if line.startswith("@"):
                    self.headers.append(line)
                else:
                    break

        fields = SAMModel().dumpable_attributes

        with open(file_name) as fh:
            for data in csv.DictReader(skip_comments(fh), fieldnames=fields, delimiter='\t', quoting=csv.QUOTE_NONE):
                _features = {}
                if data["features"]:
                    data["raw_features"] = " ".join(data["features"])
                    k, t, v = data["features"].split(":")
                    if t == "i":
                        t = int
                    else:
                        t = str
                    _features[k] = t(v)
                if data.has_key(None):
                    data["raw_features"] = " ".join(data[None])
                    for line in data[None]:
                        k, t, v = line.split(":")
                        if t == "i":
                            t = int
                        else:
                            t = str
                        _features[k] = t(v)
                    del data[None]
                data["features"] = _features
                obj = SAMModel()
                try:
                    obj.set_with_dict(data)
                except:
                    continue
                yield obj

def sc_sam_reader(sam_file):
    """
    """
    reader = SAMFileIO()
    for sam_obj in reader.read_online(sam_file):
        yield sam_obj


def sc_sam_to_fastqs(sam_files):
    """
    """
    commands = []
    for sam_file in sam_files:
        assert sam_file.endswith("sam")
        data = {
            "sam_file": sam_file,
            "fastq_file_1": sam_file.replace(".sam", "_1.fastq"),
            "fastq_file_2": sam_file.replace(".sam", "_2.fastq"),
        }
        commands.append("""cat %(sam_file)s | grep -v ^@ | awk 'NR%%2==1 {print "@"$1"\\n"$10"\\n+\\n"$11}' > %(fastq_file_1)s""" % data)
        commands.append("""cat %(sam_file)s | grep -v ^@ | awk 'NR%%2==0 {print "@"$1"\\n"$10"\\n+\\n"$11}' > %(fastq_file_2)s""" % data)
    runner.run_parallel_no_output(commands)


def sc_sam_get_headers(sam_file):
    """
    """
    reader = SAMFileIO()
    for sam_obj in reader.read_online(sam_file):
        return reader.headers


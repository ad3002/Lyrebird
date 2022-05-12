#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 21.05.2014
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com

from pymongo import MongoClient
from trseeker.tools.sequence_tools import get_revcomp
from collections import defaultdict


class AbstractKmerDBConnector(object):
    """ Abstract connector to datastore.
    """

    def __init__(self, address="localhost", port=27017):
        """ Initiation:
        - it adds mongo client
        :param address: database address
        :param port: database port
        :return: None
        """
        self.client = MongoClient("mongodb://%s:%s/" % (address, port))


class LyrebirdConnector(AbstractKmerDBConnector):
    """ Connector to LyrebirdDB.

    Fields:

        - kmer (key, in uppercase)
        - inner_name
        - inner_family
        - names
        - families
        - ref_kmer
        - hamming_dist (default: 0)
        - edit_distance (default: 0)
        - rev_comp
        - dust (default: 0)

    Collections:

        self.db - main database
        self.index - index database
        self.seqs - seqs database
        self.rfam
        self.rfam_seq
        self.repbase
        self.repbase13
        self.repbase_seq

    For human genome init with human=True.
    """

    def __init__(self, human=False, local=False):
        """ Init databases.
        :param human: force usage of human-specific databases
        :return: None
        """
        super(LyrebirdConnector, self).__init__()
        self.db = self.client.Lyrebird_HS
        if human:
            self.index = self.db.HumanIndex
            self.seqs = self.db.HumanSequence
        else:
            self.index = self.db.MainIndex
            self.seqs = self.db.LyrebirdSeqs
        self.rfam = self.db.RfamIndex
        self.rfam_seq = self.db.RfamSeqs

        self.repbase = self.db.RepbaseIndexV2
        self.repbase13 = self.db.RepbaseIndexV2_k13
        self.repbase_seq = self.db.RepbaseSeqs

        self.local = False
        if local:
            self.local = True
            self.local_index = defaultdict(dict)
            self.local_seqs = defaultdict(dict)

    # Query section

    def query(self, kmer):
        return self.query_lyrebird(kmer)

    def query_repbase(self, kmer):
        """ Query repbase.
        :param kmer: kmer
        :return: repbase item.
        """
        kmer = kmer.lower()
        data = self.repbase.find_one({'kmer': kmer})
        if not data:
            rev_kmer = get_revcomp(kmer)
            data = self.repbase.find_one({'kmer': rev_kmer})
        return data

    def query_repbase_k13(self, kmer):
        """ Query repbase for 13-mers.
        :param kmer: 17-mer
        :return: repbase item.
        """
        kmer = kmer.lower()
        data = self.repbase13.find_one({'kmer': kmer})
        if not data:
            rev_kmer = get_revcomp(kmer)
            data = self.repbase13.find_one({'kmer': rev_kmer})
        return data

    def query_lyrebird(self, kmer):
        """ Query Lyrebird.
        :param kmer: kmer
        :return: lyrebird item.
        """
        if self.local:
            kmer = kmer.upper()
            data = self.local_index[kmer]
            if not data:
                rev_kmer = get_revcomp(kmer)
                data = self.local_index[rev_kmer]
            return data

        kmer = kmer.upper()
        data = self.index.find_one({'kmer': kmer})
        if not data:
            rev_kmer = get_revcomp(kmer)
            data = self.index.find_one({'kmer': rev_kmer})
        return data

    def query_seq(self, inner_name):
        """ Query sequences.
        :param inner_name: inner name of sequence, e.g. header
        :return: sequence item.
        """
        if self.local:
            return self.local_seqs[inner_name]
            
        data = self.seqs.find_one({'inner_name': inner_name})
        return data

    def query_rfam(self, kmer):
        """ Query Rfam.
        :param kmer: kmer
        :return: Rfam item.
        """
        kmer = kmer.lower()
        data = self.rfam.find_one({'kmer': kmer})
        if not data:
            rev_kmer = get_revcomp(kmer)
            data = self.rfam.find_one({'kmer': rev_kmer})
        return data

    # Adders section

    def add(self, kmer, ref_kmer, hd, ed, name, family, rev_comp, replace=False):
        """
        Add kmer to lyrebird database
        :param kmer: kmer
        :param ref_kmer: rev_kmer
        :param hd: hamming distance from existent kmer
        :param ed: edit distance from existent kmer
        :param name: repeat name
        :param family: repeat family
        :param rev_comp: is it rev_kmer?
        :param replace: replace flag
        :return: status and obj tuple
        """
        kmer = kmer.upper()
        ref_kmer = ref_kmer.upper()

        main = {
            'kid': -1,
            'kmer': kmer,
            'inner_name': name,
            'inner_family': family,
            'names': [name],
            'families': [family],
            'ref_kmer': ref_kmer,
            'hd': hd,
            'ed': ed,
            'rev_comp': rev_comp,
            'dust': 0,
            'priority': 100,
        }

        if self.local:
            self.local_index[kmer] = main
            return "local", main

        obj = self.index.find_one({'kmer': kmer})
        if not obj:
            self.index.save(main)
            return "added", main
        else:
            if replace:
                return "skipped", main
            self.index.remove(obj["_id"])
            self.index.save(main)
            return "updated", main

    def add_update(self, data, rewrite_name=False, skip_if_known=False, db=None):
        """ Updating existing kmer with more complex logic.

        data = {
                'kid': kid,
                'kmer': kmer,
                'inner_name': inner_name,
                'inner_family': inner_family,
                'ref_kmer': kmer,
                'ed': 0,
                'hd': 0,
                'dust': get_dust_score(kmer),
                'priority': 300,
            }

        "families" : [
                "Microsatellite"
        ],
        "inner_family" : "Microsatellite",
        "ed" : 1,
        "ref_kmer" : "AAAAAAAAAAAAAAAAAAAAAAA",
        "inner_name" : "(A)n",
        "names" : [
                "(A)n"
        ],
        "dust" : 17.1,
        "hd" : 1,
        "kid" : -1,
        "kmer" : "CAAAAAAAAAAAAAAAAAAAAAA"
        "priority": 300,

        :param data: data dictionary
        :param rewrite_name: rewrite flag
        :param skip_if_known: skip if known flag
        :param db: database name
        :return: status, obj tuple
        """
        data["kmer"] = data["kmer"].upper()
        data["ref_kmer"] = data["ref_kmer"].upper()
        if db == "rfam":
            db = self.rfam
        else:
            db = self.index
        if "names" not in data:
            data["names"] = []
        if "families" not in data:
            data["families"] = [] 
        if self.local:       
            obj = self.local_index[data["kmer"]]
        else:
            obj = db.find_one({'kmer': data["kmer"]})
        if obj:
            if skip_if_known:
                return "skipped", data
            for name in obj["names"]:
                if name not in data["names"]:
                    data["names"].append(name)
            for family in obj["families"]:
                if family not in data["families"]:
                    data["families"].append(family)
            if "priority" not in obj:
                obj["priority"] = 100
            if not rewrite_name or data["priority"] < obj["priority"]:
                data["inner_name"] = obj["inner_name"]
                data["inner_family"] = obj["inner_family"]
            if not self.local:
                db.remove(obj["_id"])
        if data["inner_name"] not in data["names"]:
            data["names"].append(data["inner_name"])
        if data["inner_family"] not in data["families"]:
            data["families"].append(data["inner_family"])
        if self.local:
            self.local_index[data["kmer"]] = data
        else:
            db.save(data)
        return "added", data

    def add_repbase_sequence(self, data):
        """ Add repbase sequence

        :param data: dictionary
        :return: object
        """
        db = self.repbase_seq
        obj = db.find_one({'full_name': data["full_name"]})
        if obj:
            db.remove(obj["_id"])
        obj = db.save(data)
        return obj

    def add_rfam_sequence(self, data):
        """ Add rfam sequence

        :param data: dictionary
        :return: object
        """
        db = self.rfam_seq
        obj = db.find_one({'full_name': data["full_name"]})
        if obj:
            db.remove(obj["_id"])
        obj = db.save(data)
        return obj

    def add_lyrebird_seq(self, data):
        """ Add lyrebird sequence

        :param data: dictionary
        :return: object
        """

        if self.local:
            self.local_seqs[data["full_name"]] = data
            return data
        db = self.seqs
        obj = db.find_one({'full_name': data["full_name"]})
        if obj:
            db.remove(obj["_id"])
        obj = db.save(data)
        return obj

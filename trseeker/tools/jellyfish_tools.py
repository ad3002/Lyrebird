#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 07.09.2010
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com
"""
Jellyfish python wrapper

Used settings:
settings["blast_settings"]["jellyfish_location"]
settings["jellyfish_settings"]["hash_size"]
settings["jellyfish_settings"]["threads"]
settings["jellyfish_settings"]["both_strands"]
"""
import sys
import os
try:
    import psutil
except:
    print("Not psutil installed")

from trseeker.settings import load_settings
from PyExp import sc_iter_filepath_folder
import subprocess
from trseeker.tools.seqfile import sort_file_by_int_field
from trseeker.tools.ngrams_tools import process_list_to_kmer_index
from trseeker.tools.sequence_tools import get_revcomp
from collections import defaultdict

jellyfish_available = True
try:
    from jellyfish import jellyfish
except:
    print("Failed: from jellyfish import jellyfish")
    try:
        import jellyfish
    except:
        print("Failed: import jellyfish")
        jellyfish_available = False



class Kmer2tfAPI2(object):

    def __init__(self, jf_api):
        self.jf_api = jf_api
        self.change = defaultdict(int)

    def __getitem__(self, kmer):
        rev_kmer = get_revcomp(kmer)
        if kmer < rev_kmer:
            return self.jf_api[jellyfish.MerDNA(kmer)]
        return self.jf_api[jellyfish.MerDNA(rev_kmer)]


class Kmer2tfAPI(object):

    def __init__(self, jf_api):
        self.jf_api = jf_api
        self.change = defaultdict(int)

    def get_item_ch(self, kmer):
        rev_kmer = get_revcomp(kmer)
        if kmer < rev_kmer:
            return self.jf_api[jellyfish.MerDNA(kmer)] + self.change[kmer] + self.change[rev_kmer]
        return self.jf_api[jellyfish.MerDNA(rev_kmer)] + self.change[kmer] + self.change[rev_kmer]

    def __getitem__(self, kmer):
        rev_kmer = get_revcomp(kmer)
        if kmer < rev_kmer:
            return self.jf_api[jellyfish.MerDNA(kmer)]
        return self.jf_api[jellyfish.MerDNA(rev_kmer)]

    def increment(self, seq, k):

        for i in  xrange(len(seq)-k+1):
            self.change[seq[i:i+k]] += 1

    def decrement(self, seq, k):
        for i in  xrange(len(seq)-k+1):
            self.change[seq[i:i+k]] -= 1

    def save_changes(self, file_name):
        with open(file_name, "w") as fh:
            for k,v in self.change.items():
                fh.write("%s\t%s\n" % (k, v))

    def iter_kmers(self):
        for kmer, tf in self.jf_api:
            yield kmer, tf



class Kmer2tfAPI_cache(Kmer2tfAPI):

    def __init__(self, jf_api):
        self.jf_api = jf_api
        self.change = defaultdict(int)
        self.cache = {}

    def __getitem__(self, kmer):
        if kmer in self.cache:
            return self.cache[kmer]
        rev_kmer = get_revcomp(kmer)
        if kmer < rev_kmer:
            if not kmer in self.cache:
                self.cache[kmer] = self.jf_api[jellyfish.MerDNA(kmer)]
            return self.cache[kmer]
        if not kmer in self.cache:
            self.cache[kmer] = self.jf_api[jellyfish.MerDNA(kmer)]
        return self.cache[kmer]


        

settings = load_settings()
location = settings["blast_settings"]["jellyfish_location"]
location_new = settings["blast_settings"]["jellyfish_location_new"]


def count_kmers(input_file, output_prefix, k, mintf=None, strand=False):
    """
    Count kmers with Jellyfish.
    @param input_file: input fasta or fastq file name
    @param output_prefix: output path with file prefix
    @param k: k-mer length
    @param mintf: count only kmers with frequency greater than mintf
    @return: None
    """
    params = {
        "location": location,
        "input_fasta": input_file,
        "k": k,
        "hash_size": settings["jellyfish_settings"]["hash_size"],
        "hash_bits": settings["jellyfish_settings"]["hash_bits"],
        "threads": settings["jellyfish_settings"]["threads"],
        "both_strands": settings["jellyfish_settings"]["both_strands"],
        "output_prefix": output_prefix,
        "mintf": "",
    }
    if strand:
        params["both_strands"] = ""

    # check memory in GB
    memory = psutil.phymem_usage().total/1000000000
    if memory < 100:
        params["hash_size"] = "1G"
    if mintf:
        params["mintf"] = "--lower-count=%s" % mintf
    command = "%(location)s count %(mintf)s -m %(k)s -o %(output_prefix)s -c %(hash_bits)s -s %(hash_size)s %(both_strands)s -t %(threads)s %(input_fasta)s" % params
    print("Execute:", command)
    os.system(command)


def count_kmers_new(input_file, output_prefix, k, mintf=None, strand=False):
    """
    Count kmers with Jellyfish 2.
    @param input_file: input fasta or fastq file name
    @param output_prefix: output path with file prefix
    @param k: k-mer length
    @param mintf: count only kmers with frequency greater than mintf
    @return: None
    """
    params = {
        "location": location_new,
        "input_fasta": input_file,
        "k": k,
        "hash_size": settings["jellyfish_settings"]["hash_size"],
        "hash_bits": settings["jellyfish_settings"]["hash_bits"],
        "threads": settings["jellyfish_settings"]["threads"],
        "both_strands": settings["jellyfish_settings"]["both_strands_new"],
        "output_prefix": output_prefix,
        "mintf": "",
    }
    if strand:
        params["both_strands"] = ""
    # check memory in GB
    memory = psutil.phymem_usage().total/1000000000
    if memory < 100:
        params["hash_size"] = "1G"
    if mintf:
        params["mintf"] = "--lower-count=%s" % mintf
    command = "%(location)s count %(mintf)s -m %(k)s -o %(output_prefix)s -c %(hash_bits)s -s %(hash_size)s %(both_strands)s -t %(threads)s %(input_fasta)s" % params
    print("Execute:", command)
    os.system(command)

def merge_kmers(folder, output_prefix, output_file):
    """
    Merge Jellyfish count output to one file.
    @param folder: folder with input files
    @param output_prefix: output prefix
    @param output_file: output file
    @return: None
    """
    if not output_file.endswith(".jf"):
        output_file += ".jf"
    params = {
        "location": location,
        "output_file": output_file,
        "output_prefix": output_prefix,
    }
    file_count = 0
    for file_name in sc_iter_filepath_folder(folder, mask="."):
        # print "Check prefix", file_name, output_prefix
        if output_prefix in file_name:
            file_count += 1
    assert file_count > 0
    if file_count == 1:
        command = "cp %(output_prefix)s_0 %(output_file)s" % params
    else:
        command = "%(location)s merge -o %(output_file)s %(output_prefix)s\_*" % params
    print("Execute:", command)
    os.system(command)
    command = "rm %(output_prefix)s_*" % params
    print(command)
    os.system(command)


def stats_kmers(db_file, stats_file, new=True):
    """
    Compute statistics for kmers.
    @param db_file: jf db file
    @param stats_file: output stats file
    @return: None
    """
    params = {
        "location": location,
        "db_file": db_file,
        "stats_file": stats_file,
    }
    if new:
        params["location"] = location_new
    command = "%(location)s stats --verbose -o %(stats_file)s %(db_file)s" % params
    print("Execute:", command)
    os.system(command)


def histo_kmers(db_file, histo_file, new=True):
    """
    Compute frequencies histogram.
    @param db_file: jf db file
    @param histo_file: histogram output file
    @return: None
    """
    params = {
        "location": location,
        "db_file": db_file,
        "histo_file": histo_file,
    }
    if new:
        params["location"] = location_new
    command = "%(location)s histo --verbose -o %(histo_file)s -h 1000000000 %(db_file)s" % params
    print("Execute:", command)
    os.system(command)


def dump_kmers(db_file, kmers_file, dumpmintf, new=False):
    """
    Dump and sort k-mers database to tab-delimited file.
    @param db_file: jf db file
    @param kmers_file: output kmer file
    @param dumpmintf: minimum tf to dump
    @return: None
    """
    params = {
        "location": location,
        "db_file": db_file,
        "kmers_file": kmers_file,
        "dumpmintf": dumpmintf,
    }
    if new:
        params["location"] = location_new
    command = "%(location)s dump --column -L %(dumpmintf)s --tab -o %(kmers_file)s %(db_file)s" % params
    print("Execute:", command)
    os.system(command)
    if os.path.getsize(kmers_file)/1000000000 <10:
        sort_file_by_int_field(kmers_file, 1)

def query_kmers(db_file, query_hashes, both_strands=True, verbose=False, new=False, batch_size=1000):
    """
    Query jellyfish database.
    @param db_file: jf db file
    @param query_hashes: kmers to query
    @param both_strands: use both strands
    @param verbose: verbose
    @return: dictionary hash to tf
    """
    if len(query_hashes[0]) > 23 or new:
        return query_kmers_new(db_file, query_hashes, both_strands=True, verbose=verbose)

    params = {
        "location": location,
        "db_file": db_file,
        "query_hashes": query_hashes,
        "both_strands": "",
    }
    if both_strands:
        params["both_strands"] = "-C"
    command = "%(location)s query %(both_strands)s %(db_file)s" % params
    if verbose > 1:
        print(command)
    final_result = defaultdict(int)
    n = len(query_hashes)
    step = batch_size
    for k in xrange(0,n,step):
        if verbose > 1:
            print(k, n)
        pp = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE, universal_newlines=True)
        # query = " ".join(query_hashes[k:k+100]).upper()
        for i, query in enumerate(query_hashes[k:k+step]):
            pp.stdin.write(query+" ")
        data = pp.communicate(input=query)
        # sys.stdout.flush()
        error = data[1]
        if "Can't open file" in error:
            return None
        if verbose > 1:
            print("Data size", len(data[0]))
        for item in data[0].strip().split("\n"):
            if item:
                if verbose > 2:
                    print(item)
                key, value = item.upper().strip().split()
                final_result[key] = int(value)
            else:
                if verbose > 2:
                    print(data)
                final_result[-1] = data[1]
        
    return(final_result)

def query_kmers_new(db_file, query_hashes, both_strands=True, verbose=False):
    """
    Query jellyfish database.
    @param db_file: jf db file
    @param query_hashes: kmers to query
    @param both_strands: use both strands
    @param verbose: verbose
    @return: dictionary hash to tf
    """
    params = {
        "location": location_new,
        "db_file": db_file,
        "query_hashes": query_hashes,
        "both_strands": "",
    }
    # if both_strands:
    #     params["both_strands"] = "-C"
    command = "%(location)s query -i %(both_strands)s %(db_file)s" % params
    if verbose:
        print(command)
    final_result = defaultdict(int)
    n = len(query_hashes)

    step = 1000
    for k in xrange(0,n,step):
        if verbose > 1:
            print(k, n)
        pp = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE, universal_newlines=True)
        # query = " ".join(query_hashes[k:k+100]).upper()
        for i, query in enumerate(query_hashes[k:k+step]):
            pp.stdin.write(query+" ")
        data = pp.communicate(input=query)
        error = data[1]
        if "Can't open file" in error:
            return None
        for item in data[0].strip().split("\n"):
            if item.strip():
                if verbose:
                    print(k, n, item)
                final_result[query] = int(item.strip().split()[0])
            else:
                if verbose:
                    print(k, n, data)
                final_result[-1] = data[1]
    return final_result


def query_and_write_coverage_histogram(db_file, query_sequence, output_file, k=23):
    """
    Save coverage histogram into output_file for given query_sequence.
    @param db_file:
    @param query_sequence:
    @param output_file:
    @param k:
    @return:
    """
    index = process_list_to_kmer_index([query_sequence], k, docids=False, cutoff=None, verbose=False)
    query_hashes = [x[0] for x in index]
    data =  query_kmers(db_file, query_hashes, both_strands=True)
    keys = data.keys()
    for kmer in keys:
        rkmer = get_revcomp(kmer)
        data[kmer] = int(data[kmer])
        if not rkmer in data:
            data[rkmer] = int(data[kmer])
    with open(output_file, "w") as fh:
        for i in xrange(0, len(query_sequence)-k + 1):
            kmer = query_sequence[i:i+k]
            p0 = 0
            p1 = 0
            p2 = 0
            rkmer = get_revcomp(kmer)
            if kmer in data:
                p1 = int(data[kmer])
            if rkmer in data:
                p2 = int(data[rkmer])
            p = max(p1,p2,p0)
            fh.write("%s\t%s\t%s\n" % (kmer, i, p))
    return data

def get_sequence_coverage(db_file, query_sequence, k=23):
    """
    Get coverage histogram for given query_sequence.
    @param db_file:
    @param query_sequence:
    @param k:
    @return:
    """
    query_sequence = query_sequence.upper()
    index = process_list_to_kmer_index([query_sequence], k, docids=False, cutoff=None, verbose=False)
    query_hashes = [x[0].upper() for x in index]
    data =  query_kmers(db_file, query_hashes, both_strands=True)
    keys = data.keys()
    for kmer in keys:
        rkmer = get_revcomp(kmer)
        data[kmer] = int(data[kmer])
        if not rkmer in data:
            data[rkmer] = int(data[kmer])
    
    coverage = []
    for i in xrange(0, len(query_sequence)-k + 1):
        kmer = query_sequence[i:i+k]
        p0 = 0
        p1 = 0
        p2 = 0
        rkmer = get_revcomp(kmer)
        if kmer in data:
            p1 = int(data[kmer])
        if rkmer in data:
            p2 = int(data[rkmer])
        p = max(p1,p2,p0)
        coverage.append(p)
    return coverage


def sc_count_and_dump_kmers_for_folder(folder, output_prefix, kmers_file, k=23, mintf=None):
    """
    Count kmers, merge them, and save to tab-delimited kmers_file.
    @param folder: folder with input files
    @param output_prefix: prefix for input files
    @param kmers_file: output dump file
    @param k: kmer length
    @param mintf: minimal tf for count
    @return: None
    """
    count_kmers(input_file, output_prefix, k, mintf=mintf)
    merge_kmers(folder, input_file, input_file)
    db_file = "%s.jf" % input_file
    dump_kmers(db_file, kmers_file)

def sc_count_and_dump_kmers_for_file(fasta_file, jellyfish_data_folder, jf_db, jf_dat, k, mintf, dumpmintf, strand=False, use_new=True):
    """
    Count kmers, merge them, and save to tab-delimited kmers_file.
    @param fasta_file: fasta file
    @param jellyfish_data_folder: output jellyfish folder
    @param jf_db: jf database output file
    @param jf_dat: jf dump output file
    @param k: kmer length
    @param mintf: minimal tf for count
    @param dumpmintf: minimal tf for dump
    @return: None
    """
    ouput_prefix = os.path.join(
            jellyfish_data_folder,
            "%s___" % jf_db,
        )

    jf_db = os.path.join(jellyfish_data_folder, jf_db)
    jf_dat = os.path.join(jellyfish_data_folder, jf_dat)
    
    count_kmers_new(fasta_file, jf_db, k, mintf=mintf, strand=strand)
    if dumpmintf >= 0:
        dump_kmers(jf_db, jf_dat, dumpmintf=dumpmintf, new=True)

    if strand or dumpmintf < 0:
        return

    print("Sort data...")
    temp_file = jf_dat+".temp"
    data = {
        "in": jf_dat,
        "out": temp_file,
    }
    command = "sort -k2nr %(in)s > %(out)s" % data
    print(command)
    os.system(command)
    command = "mv %(out)s %(in)s" % data
    print(command)
    os.system(command)
    for file_path in sc_iter_filepath_folder(jellyfish_data_folder):
        if "%s___" % jf_db in file_path:
            print("Removing", file_path)
            os.unlink(file_path)

def load_kmer2tf(jf_data_file):
    """ Load kmer2tf file from jellyfish dat file
    """
    print("Read file %s..." % jf_data_file)
    kmer2tf = {}
    with open(jf_data_file) as fh:
        print("Read data...")
        kmer2tf = fh.readlines()
        print("Format data...")
        kmer2tf = [x.split("\t") for x in kmer2tf]
        print("Convert to dict..")
        kmer2tf = dict(kmer2tf)
    print("Done.")
    return kmer2tf


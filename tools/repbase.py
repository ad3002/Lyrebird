#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 21.05.2014
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com

import sys
sys.path.append("/Users/akomissarov/Dropbox/workspace")

from trseeker.tools.ngrams_tools import process_list_to_kmer_index
from trseeker.seqio.fasta_file import sc_iter_fasta
from trseeker.tools.sequence_tools import get_dust_score
from PyBioSnippets.Lyrebird.tools.uploaders import *


def get_known_in_repbase_kmers(kmers, repbase_connector):
	''' Retrun a dictionary kmer to list of repbase hits.
	'''
	repbase_kmers = {}
	found = 0
	for i, kmer in enumerate(kmers):
		res = repbase_connector.query(kmer)
		if res:
			repbase_kmers[kmer] = res
			found += 1
			print found, i, kmer, res[0]
		else:
			print found, i, kmer, "-"
	return repbase_kmers


def get_repbase_coverage(repbase_name, jellyfish_db):
	''' Check coverage of given repeat in repbase.
	'''
	pass


def get_known_in_sine_kmers(kmers, k):
	''' Retrun a dictionary kmer to list of repbase hits.
	'''
	dataset = upload_sinebase(k=k)
	repbase_kmers = {}
	found = 0
	for i, kmer in enumerate(kmers):
		res = None
		if kmer in dataset:
			res = dataset[kmer]
		if res:
			repbase_kmers[kmer] = res
			found += 1
			print found, i, kmer, res["repbase_names"][0]
		# else: 
		# 	print found, i, kmer, "-"
	return repbase_kmers



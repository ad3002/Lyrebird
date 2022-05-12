#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 07.09.2013
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com
'''
'''
import os

def build_bowtie2_db(input_file, output_file):
	'''
	'''
	data = {
		"input_file": input_file,
		"output_file": output_file,
	}
	command = "bowtie2-build %(input_file)s %(output_file)s" % data
	print command
	os.system(command)

def map_with_bowtie2(fastq1, fastq2, unpaired, db, sam_output, local=None):
	'''
	'''

	if isinstance(fastq1, list):
		fastq1 = ",".join(fastq1)
	if isinstance(fastq2, list):
		fastq2 = ",".join(fastq2)
	if isinstance(unpaired, list):
		unpaired = ",".join(unpaired)

	data = {
		"fastq1": fastq1,
		"fastq2": fastq2,
		"db": db,
		"sam": sam_output,
		"unpaired": unpaired,
	}
	if not local:
		data["mode"] = "--very-sensitive"
	else:
		data["mode"] = "--very-sensitive-local"
	if os.path.isfile(sam_output):
		answer = raw_input("File exists (%s), overwrite? (skip)" % sam_output) or None
		if answer is None:
			return
	if sam_output.endswith("bam"):
		if unpaired:
			command = "bowtie2 -x %(db)s -p 50 -1 %(fastq1)s -2 %(fastq2)s -U %(unpaired)s %(mode)s  | samtools view -bS - > %(sam)s" % data
		else:
			command = "bowtie2 -x %(db)s -p 50 -1 %(fastq1)s -2 %(fastq2)s %(mode)s | samtools view -bS - > %(sam)s" % data
		print command
	if sam_output.endswith("sam"):
		if unpaired:
			command = "bowtie2 -x %(db)s -p 50 -1 %(fastq1)s -2 %(fastq2)s -U %(unpaired)s -S %(sam)s %(mode)s" % data
		else:
			command = "bowtie2 -x %(db)s -p 50 -1 %(fastq1)s -2 %(fastq2)s -S %(sam)s %(mode)s" % data
		print command
	os.system(command)
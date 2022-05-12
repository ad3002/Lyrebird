#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 26.09.2010
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com 
'''
Handle genome project data from NCBI
'''
import urllib, csv
from collections import defaultdict
# from PyExp.models.abstract_model import AbstractModel

# class NCBIGenomeProjectBacteria(AbstractModel):
    
#     dumpable_attributes = [
#         "Organism/Name",
#         "Kingdom",
#         "Group",
#         "SubGroup"
#         "Size (Mb)",
#         "Chrs",
#         "Organelles",
#         "Plasmids",
#         "BioProjects",
#     ]

# class NCBIGenomeProjectEukaryota(AbstractModel):
    
#     dumpable_attributes = [
#         "Organism/Name",
#         "BioProject Accession",
#         "BioProject ID",
#         "Kingdom",
#         "Group",
#         "SubGroup"
#         "Size (Mb)",
#         "GC%",
#         "Assembly",
#         "Chromosomes",
#         "Organelles",
#         "Plasmids",
#         "WGS",
#         "Scaffolds",
#         "Genes",
#         "Proteins",
#         "Release Date",
#         "Modify Date",
#         "Status",
#         "Center",
#     ]                       


def load_genome_projects(subgroup):
    ''' Scrap finished WGS genome projects wtih given SubGroup name.
    '''
    url = "ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/eukaryotes.txt"
    data = urllib.urlopen(url).readlines()
    data[0] = data[0][1:]
    genome_projects = []
    for line in csv.DictReader(data, delimiter='\t'):
        if line["SubGroup"]in subgroup and line["WGS"] != "-":
            name = line["Organism/Name"].replace(".","").replace("(","").replace(")","").replace("-","_")
            wgs = line["WGS"].replace("01","").replace("02","").replace("03","").replace("04","")
            genome_projects.append((name, line["Status"], line["Chromosomes"], wgs))
    return genome_projects

def load_genome_projects_exclude(notsubgroups):
    ''' Scrap finished WGS genome projects wtih given SubGroup name.
    '''
    url = "ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/eukaryotes.txt"
    data = urllib.urlopen(url).readlines()
    data[0] = data[0][1:]
    genome_projects = []
    r = defaultdict(int)
    for line in csv.DictReader(data, delimiter='\t'):
        if line["SubGroup"] not in notsubgroups and line["WGS"] != "-":
            #print "%s\t%s\t%s\t%s" % (line["Organism/Name"], line["Status"], line["Chromosomes"], line["WGS"])
            name = line["Organism/Name"].replace(".","").replace("(","").replace(")","").replace("-","_")
            wgs = line["WGS"].replace("01","").replace("02","").replace("03","").replace("04","")
            genome_projects.append((name, line["Status"], line["Chromosomes"], wgs))
            r[line["SubGroup"]] += 1
    return genome_projects

def print_add_project_data(genome_projects, pid_type, pr_type):
    ''' Print data for initiating PySatDNA projects.
    TODO: fixit
    '''
    name2wgs = {}
    for item in genome_projects:
        name = "_".join(item[0].lower().split()) + "_wgs"
        name2wgs.setdefault(name, [])
        name2wgs[name].append(item[-1])

    seen = set()
    for item in genome_projects:
        name = "_".join(item[0].lower().split()) + "_wgs"
        if name in seen:
            continue
        seen.add(name)
        wgs = item[-1]
        pid = "%s_%s" % (pid_type, name)
        words = item[0].split()
        name_short = words[0][0] + ". "+ " ".join(words[1:])
        suffix = "".join([x[0].upper() for x in item[0].split()])
        project = {
            'path_to': "/home/%s/%s" % (pr_type, name),
            'title': name,
            'taxon': item[0],
            'taxon_short': name_short,
            'fasta_suffix': "fsa_nt",
            'wgs_projects': list(set(name2wgs[name])),
            'large_tr_cutoff': 3000,
            'suffix': suffix,
            'type': 'wgs',
            'pid': pid,
        }

        print "pid_%s = '%s'" % (pid, pid)
        print "project_data_%s = " % pid, project

    print "all_%s = [" % pid_type
    seen = set()
    for item in genome_projects:
        name = "_".join(item[0].lower().split()) + "_wgs"
        if name in seen:
            continue
        seen.add(name)
        pid = "%s_%s" % (pid_type, name)
        print  "(pid_%s, project_data_%s),\n" % (pid, pid)
    print "]"
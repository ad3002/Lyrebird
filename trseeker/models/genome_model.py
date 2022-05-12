#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 25.12.2010
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com 

from PyExp import AbstractModel


class GenomeModel(AbstractModel):
    ''' Container for genome information according 
    to NCBI fields from ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/eukaryotes.txt

    '''
    dumpable_attributes = [
        "genome_taxon",
        "genome_taxid",
        "genome_bioproject",
        "genome_bpid",
        "genome_group",
        "genome_subgroup",
        "genome_size_mb",
        "genome_gc",
        "genome_assembly_ref",
        "genome_chromosomes",
        "genome_organelles",
        "genome_plasmids",
        "genome_wgs",
        "genome_scaffolds",
        "genome_genes",
        "genome_proteins",
        "genome_release_date",
        "genome_modify_date",
        "genome_status",
        "genome_center",
        "genome_biosample_ref",
        ]


class RepeatMaskerAnnotation(AbstractModel):
    ''' Container for RepeatMasker track.
    For example from http://hgdownload.cse.ucsc.edu/goldenPath/felCat5/database/rmsk.txt.gz

    Table fields:
      `bin` smallint(5) unsigned NOT NULL,
      `swScore` int(10) unsigned NOT NULL,
      `milliDiv` int(10) unsigned NOT NULL,
      `milliDel` int(10) unsigned NOT NULL,
      `milliIns` int(10) unsigned NOT NULL,
      `genoName` varchar(255) NOT NULL,
      `genoStart` int(10) unsigned NOT NULL,
      `genoEnd` int(10) unsigned NOT NULL,
      `genoLeft` int(11) NOT NULL,
      `strand` char(1) NOT NULL,
      `repName` varchar(255) NOT NULL,
      `repClass` varchar(255) NOT NULL,
      `repFamily` varchar(255) NOT NULL,
      `repStart` int(11) NOT NULL,
      `repEnd` int(11) NOT NULL,
      `repLeft` int(11) NOT NULL,
      `id` char(1) NOT NULL,
    '''
    
    dumpable_attributes = [
        'bin',
        'swScore',
        'milliDiv',
        'milliDel',
        'milliIns',
        'genoName',
        'genoStart',
        'genoEnd',
        'genoLeft',
        'strand',
        'repName',
        'repClass',
        'repFamily',
        'repStart',
        'repEnd',
        'repLeft',
        'id',
    ]
    int_attributes = [
        'bin',
        'swScore',
        'milliDiv',
        'milliDel',
        'milliIns',
        'genoStart',
        'genoEnd',
        'genoLeft',
        'repStart',
        'repEnd',
        'repLeft',
    ]
#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 14.06.2011
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com 

from PyExp import AbstractModel
from trseeker.tools.parsers import trf_parse_param, trf_parse_head, \
    parse_fasta_head, parse_chromosome_name, trf_parse_line
from trseeker.tools.sequence_tools import clear_sequence, get_gc


class TRModel(AbstractModel):
    """ Class for tandem repeat wrapping
    
    Public properties:

    Indexes:
    
    - id, 
    - giid, 
    - trf_id, 
    - project, 
    - trf_gi
    
    Coordinates:  
        
    - trf_l_ind, 
    - trf_r_indel
    
    TR:

    - trf_period, 
    - trf_n_copy, 
    - trf_pmatch,
    - trf_pvar,
    - trf_entropy,
    - trf_array_length, 
    - trf_joined, 
    - trf_chr

    Sequence:

    - trf_consensus, 
    - trf_array

    GC%:
    
    - trf_array_gc,
    - trf_consensus_gc
    
    Other:
    
    - trf_head, 
    - trf_params
    
    Annotation:
    
    - trf_family, 
    - trf_subfamily, 
    - trf_subsubfamily,
    - trf_hor, 
    - trf_n_chrun, 
    - trf_n_refgenome,
    - trf_chr_refgenome, 
    - trf_bands_refgenome,
    - trf_repbase, 
    - trf_strand
    
    Dumpable atributes:

    - "project",
    - "id" (int),
    - "trf_id" (int),
    - "trf_type",
    - "trf_family",
    - "trf_family_prob",
    - "trf_l_ind" (int),
    - "trf_r_ind" (int),
    - "trf_period" (int),
    - "trf_n_copy" (float),
    - "trf_pmatch" (float),
    - "trf_pvar" (float),
    - "trf_consensus",
    - "trf_array",
    - "trf_array_gc" (float),
    - "trf_consensus_gc" (float),
    - "trf_gi",
    - "trf_head",
    - "trf_param",
    - "trf_array_length" (int),
    - "trf_chr",
    - "trf_joined" (int),
    - "trf_superfamily",
    - "trf_superfamily_ref",
    - "trf_superfamily_self",
    - "trf_subfamily",
    - "trf_subsubfamily",
    - "trf_family_network",
    - "trf_family_self",
    - "trf_family_ref",
    - "trf_hor" (int),
    - "trf_n_chrun" (int),
    - "trf_chr_refgenome",
    - "trf_bands_refgenome",
    - "trf_repbase",
    - "trf_strand",

    """

    dumpable_attributes = ["project",
                           "id",
                           "trf_id",
                           "trf_type",
                           "trf_family", # set for alpha in cf_premask_alpha
                           "trf_family_prob",
                           "trf_l_ind",
                           "trf_r_ind",
                           "trf_period",
                           "trf_n_copy",
                           "trf_pmatch",
                           "trf_pvar",
                           "trf_entropy",
                           "trf_consensus",
                           "trf_array",
                           "trf_array_gc",
                           "trf_consensus_gc",
                           "trf_gi",
                           "trf_head",
                           "trf_param",
                           "trf_array_length",
                           "trf_chr",
                           "trf_joined",
                           "trf_superfamily",
                           "trf_superfamily_ref",
                           "trf_superfamily_self",
                           "trf_subfamily",
                           "trf_subsubfamily",
                           "trf_family_network",
                           "trf_family_self",
                           "trf_family_ref",
                           "trf_hor",
                           "trf_n_chrun",
                           "trf_chr_refgenome",
                           "trf_bands_refgenome",
                           "trf_repbase",
                           "trf_strand",
                           ]

    int_attributes = [
                           "trf_id",
                           "id",
                           "trf_l_ind",
                           "trf_r_ind",
                           "trf_period",
                           "trf_array_length",
                           "trf_joined",
                           "trf_hor",
                           "trf_n_chrun",
                       ]

    float_attributes = [
                           "trf_n_copy",
                           "trf_family_prob",
                           "trf_entropy",
                           "trf_pmatch",
                           "trf_pvar",
                           "trf_array_gc",
                           "trf_consensus_gc",
                        ]

    def set_project_data(self, project):
        ''' Add project data to self.project.'''
        self.project = project

    def set_raw_trf(self, head, body, line):
        """ Init object with data from parsed trf ouput.
        
        Parameters:
        
        - trf_obj: TRFObj instance
        - head: parsed trf head
        - body: parsed trf body
        - line: parsed trf line
        """
        self.trf_param = 0
        self.trf_head = trf_parse_head(head).strip()
        self.trf_gi = parse_fasta_head(self.trf_head)[0]
        self.trf_chr = parse_chromosome_name(self.trf_head)

        (self.trf_l_ind,
         self.trf_r_ind,
         self.trf_period,
         self.trf_n_copy,
         self.trf_l_cons,
         self.trf_pmatch,
         self.trf_indels,
         self.trf_score,
         self.trf_n_a,
         self.trf_n_c,
         self.trf_n_g,
         self.trf_n_t,
         self.trf_entropy,
         self.trf_consensus,
         self.trf_array) = trf_parse_line(line)

        self.trf_pvar = int(100 - float(self.trf_pmatch))

        try:
            self.trf_l_ind = int(self.trf_l_ind)
        except:
            print(self)

        self.trf_r_ind = int(self.trf_r_ind)
        self.trf_period = int(self.trf_period)
        self.trf_n_copy = float(self.trf_n_copy)

        self.trf_consensus = clear_sequence(self.trf_consensus)
        self.trf_array = clear_sequence(self.trf_array)

        self.trf_array_gc = get_gc(self.trf_array)
        self.trf_consensus_gc = get_gc(self.trf_consensus)
        self.trf_chr = parse_chromosome_name(self.trf_head)
        self.trf_array_length = len(self.trf_array)


    def get_index_repr(self):
        """ Get string for index file."""
        return "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (self.trf_id,
                                                         self.trf_period,
                                                         self.trf_array_length,
                                                         self.trf_array_gc,
                                                         self.trf_pvar,
                                                         self.trf_gi,
                                                         self.trf_l_ind,
                                                         self.trf_r_ind,
                                                         self.trf_chr)

    def get_numerical_repr(self):
        """ Get str for Mathematica."""
        return "%s\t%s\t%.2f\n" % (self.trf_period,
                                   self.trf_array_length,
                                   self.trf_array_gc)

    def get_fasta_repr(self, add_project=False):
        ''' Get array fasta representation, head - trf_id.'''
        if add_project:
          return ">%s_%s\n%s\n" % (self.trf_id, self.project, self.trf_array)
        else:
          return ">%s\n%s\n" % (self.trf_id, self.trf_array)

    def get_monomer_fasta_repr(self):
        ''' Get monomer fasta representation, head - trf_id.'''
        return ">%s\n%s\n" % (self.trf_id, self.trf_consensus)

    def get_family_repr(self):
        """ Get str for family index."""
        return "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (self.trf_id,
                                                     self.trf_period,
                                                     self.trf_array_length,
                                                     self.trf_array_gc,
                                                     self.trf_pvar,
                                                     self.trf_gi,
                                                     self.trf_l_ind,
                                                     self.trf_r_ind,
                                                     self.trf_chr,
                                                     self.trf_repbase,
                                                     self.trf_superfamily,
                                                     self.trf_family,
                                                     self.trf_subfamily)

    @property
    def fasta(self):
      return self.get_fasta_repr()

    def get_gff3_string(self, 
                        chromosome=True, 
                        trs_type="complex_tandem_repeat", 
                        probability=1000,
                        tool="PySatDNA",
                        prefix=None,
                        properties=None,
                        force_header=False):
      '''Return TR in gff format.
      '''
      if chromosome and self.trf_chr and self.trf_chr != "?":
          seqid = self.trf_chr
      elif self.trf_gi and self.trf_gi != "Unknown":
          seqid = self.trf_gi
      else:
          seqid = self.trf_head
      features = []
      if not properties:
          properties = {}
      for name, attr in properties.items():
          features.append("%s=%s" % (
                  name,
                  getattr(self, attr)
                )
          )
      features = ";".join(features)
      if prefix:
        seqid = prefix + seqid
      if self.trf_l_ind < self.trf_r_ind:
        strand = "+"
      else:
        strand = "-"
        self.trf_l_ind, self.trf_r_ind = self.trf_r_ind, self.trf_l_ind

      if force_header:
        seqid = self.trf_head
      d = (seqid, 
           tool, 
           trs_type, 
           self.trf_l_ind, 
           self.trf_r_ind, 
           probability, 
           strand, 
           ".", 
           features,
          )
      return "%s\n" % "\t".join(map(str, d))

    def get_bed_string(self):
      '''Return TR in bed format.
      '''
      if self.trf_l_ind < self.trf_r_ind:
        strand = "+"
      else:
        strand = "-"
        self.trf_l_ind, self.trf_r_ind = self.trf_r_ind, self.trf_l_ind

      seqid = self.trf_head
      d = (seqid, 
           self.trf_l_ind, 
           self.trf_r_ind, 
          )
      return "%s\n" % "\t".join(map(str, d))

class NetworkSliceModel(TRModel):
    """ Class for network slice data.
    """
    def __init__(self):
      self.dumpable_attributes = ["gid"] + self.dumpable_attributes
      self.int_attributes = ["gid"] + self.int_attributes
      super(TRModel, self).__init__()

class TRsClassificationModel(AbstractModel):
  ''' Model for keeping classification data.
  '''
  dumpable_attributes = ["project",
                           "id",
                           "trf_id",
                           "trf_period",
                           "trf_array_length",
                           "trf_array_gc",

                           "trf_type",
                           "trf_family",
                           "trf_subfamily",
                           "trf_family_prob",
                           "trf_family_kmer",
                           "trf_subfamily_kmer",
                           "trf_family_self",

                           "class_ssr", # from cf_find_ssr
                           "class_tssr", # cf_find_ssr
                           "class_sl", # cf_find_ssr, a number of loci
                           "class_good", # cf_trs_classififcaiton_by_type
                           "class_micro", # cf_trs_classififcaiton_by_type
                           "class_100bp", # cf_trs_classififcaiton_by_type
                           "class_perfect", # cf_trs_classififcaiton_by_type
                           "class_x4", # cf_trs_classififcaiton_by_type
                           "class_entropy", # cf_trs_classififcaiton_by_type
                           "class_gc", # cf_trs_classififcaiton_by_type
                           "trf_consensus", 
                ]

  alt_dumpable_attributes = ["project",
                           "id",
                           "trf_id",
                           "trf_period",
                           "trf_array_length",
                           "trf_array_gc",
                           "trf_type",
                           "trf_family",
                           "trf_subfamily",
                           "trf_family_prob",
                           "class_ssr",
                           "class_sl",
                           "class_good",
                           "class_micro",
                           "class_100bp",
                           "class_perfect",
                           "class_x4",
                           "class_entropy",
                           "class_gc",
                           "trf_consensus",
                ]

  int_attributes = [
                         "trf_id",
                         "id",
                         "trf_period",
                         "trf_array_length",
                     ]

  float_attributes = [
                         "trf_family_prob",
                      ]

  def set_with_trs(self, trf_obj):
    '''Set data with trf_obj.
    '''
    for attr in self.dumpable_attributes:
      if hasattr(trf_obj, attr):
        setattr(self, attr, getattr(trf_obj, attr))

  @property
  def network_head(self):
    '''Return TRs representation for network slices.
    '''
    dumpable_attributes = [
                           "trf_id",
                           "trf_period",
                           "trf_array_length",
                           "trf_array_gc",
                           "trf_type",
                           "trf_family",
                           "trf_subfamily",
                           "trf_family_prob",
                           "trf_family_kmer",
                           "trf_subfamily_kmer",
                           "trf_family_self",
                           "class_ssr",
                           "class_tssr",
                           "class_sl",
                           "class_good",
                           "class_micro",
                           "class_100bp",
                           "class_perfect",
                           "class_x4",
                           "class_entropy",
                           "class_gc",
                           "trf_consensus",
                ]
    return self.get_as_string(dumpable_attributes)






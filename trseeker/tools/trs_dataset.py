#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 04.06.2010
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com 
'''
Functions related to TRs annotation.
'''

import os
from trseeker.models.trf_model import TRModel
from trseeker.seqio.tab_file import sc_iter_tab_file

BUBBLE_CHART_TEMPLATE = "BubbleChart3D[{%s}]"
LEGENDED_CHART_TEMPLATE = 'Legended[Style[Import["%(file)s"], RGBColor[%(color_r)s, %(color_g)s, %(color_b)s]], "%(title)s"]'

COLORS_50 = {
  "denim": ("#1560bd", (21,96,189)),
}

class Legend(object):
    ''' Class for description color legend.'''
    file = ""
    color_r = 0.3
    color_g = 0.3
    color_b = 0.3
    title = "None"

    def as_dict(self):
        return {"file":self.file,
                "color_r":self.color_r,
                "color_g":self.color_g,
                "color_b":self.color_b,
                "title":self.title,
                }

COLORS_POS = [
              "ALPHA",
              "CEN",
              "TR1",
              "TR2",
              "TR3",
              "TR4",
              "TR5",
              "TR6",
              "TR7",
              "TR8",
              "TR9",
              "TR10",
              "TR11",
]

COLORS = {
              "ALPHA":(0.33, 0.36, 0.75),
              "CEN":(208 / 255.0, 32 / 255.0, 144 / 255.0),
              "TR1": (255 / 255.0, 127 / 255.0, 36 / 255.0),
              "TR2":(0 / 255.0, 0 / 255.0, 128 / 255.0),
              "TR3":(127 / 255.0, 255 / 255.0, 212 / 255.0),
              "TR4":(255 / 255.0, 193 / 255.0, 193 / 255.0),
              "TR5":(255 / 255.0, 64 / 255.0, 64 / 255.0),
              "TR6":(33 / 255.0, 144 / 255.0, 255 / 255.0),
              "TR7":(105 / 255.0, 139 / 255.0, 35 / 255.0),
              "TR8":(1, 215 / 255.0, 0),
              "TR9":(202 / 255.0, 255 / 255.0, 112 / 255.0),
              "TR10":(204 / 255.0, 102 / 255.0, 153 / 255.0),
              "TR11":(102 / 255.0, 0 / 255.0, 0 / 255.0),
              "TE":(50 / 255.0, 205 / 255.0, 50 / 255.0),
              "OTHER":(0.7, 0.7, 0.7),
}

MAP_DATASET_TO_GF = {
  "MaSat": "ALPHA",
  "Alpha": "ALPHA",
  "TR-92A-RN": "ALPHA",
  "alpha735": "ALPHA", # Bsp

  "MiSat": "CEN",
  "TR-172A-HS": "CEN",
  "TR-196A-RN": "CEN",
  "74A": "CEN", # SAT1_CF
  'TR-37B-CLF': "CEN",

  "TR-21A-MM": "TR1",
  "TR-21B-MM": "TR1",
  "HS2": "TR1",
  "TR-5A-HS": "TR1",
  "TR-84A-RN": "TR1",
  'TR-70A-CLF': "TR1",

  "TR-84A-MM": "TR2",
  "TR-84B-HS": "TR2",
  "TR-84B-RN": "TR2",
  "TR-84A-CLF": "TR2",

  "TR-31A-MM": "TR3",
  "TR-76A-HS": "TR3",
  'TR-41A-RN': "TR3",
  'TR-37A-CLF': "TR3",
  
  "TR-31B-MM": "TR4",
  "TR-68A-HS": "TR4",
  "TR-30A-RN": "TR4",
  "TR-10A-CLF": "TR4",
  
  "TR-22A-MM": "TR5",
  "TR-84A-HS": "TR5",
  "TR-100A-RN": "TR5",
  "TR-8A-CLF": "TR5",

  "TR-38B-MM": "TR6",
  "TR-217A-HS": "TR6",
  "TR-188A-HS": "TR6",
  "TR-211A-HS": "TR6",
  "TR-20A-RN": "TR6",
  "TR-24A-CLF": "TR6",

  "TR-54A-MM": "TR7",
  "TR-48A-HS": "TR7",
  "TR-10A-RN": "TR7",
  "TR-42A-CLF": "TR7",


  "TR-27B-MM": "TR8",
  "TR-40B-HS": "TR8",
  'TR-65A-RN': "TR8",
  'TR-6B-CLF': "TR8",

  "TR-24C-MM": "TR9",
  'TR-47A-HS': "TR9",
  'TR-11A-RN': "TR9",
  "TR-248A-CLF": "TR9", #CarSat1

  "TR-208A-MM": "TE",
  "TR-212A-MM": "TE",
  "TR-218A-MM": "TE",
  
  'TR-125A-HS': "TE",
  "TR-218A-HS": "TE",
  "TR-210A-RN": "TE",
  "TR-205A-RN": "TE",
}              



TRANSLATE_FAMS = {
  "TR-120A-MM": "MiSat",
  "TR-21B-MM": "TR-21A-MM",
  
  "alpha": "Alpha",
  "hsatII": "HS2",
}

def get_colors(family):
    ''' Return color according to family.
    TODO: fix it.
    '''
    if family in COLORS:
      return COLORS[family]
    return COLORS["OTHER"]

def get_colors_rgb(family):
    ''' Get RGB color according to family.'''
    colors = {"OTHER":"#ba4017",
              "ALPHA":"#535bc7",
              "CEN":"#e1004c",
              "BASE":"#ff8f00",
              "ML": "#291763",
              "SL": "#ffff40",
              "TE":"#2f9414",
              }
    return colors[family]


def join_legends(legends):
    ''' Join legend data for visualization.'''
    data = [LEGENDED_CHART_TEMPLATE % x.as_dict() for x in legends]
    return ",".join(data)

def get_family_name(trf_obj):
    family_network = trf_obj.trf_family_self
    if family_network is None:
        family_network = trf_obj.trf_family
    if family_network in TRANSLATE_FAMS:
        family_network = TRANSLATE_FAMS[family_network]
    if family_network in MAP_DATASET_TO_GF:
        family_network = MAP_DATASET_TO_GF[family_network]
    return family_network

def create_mathematice_dataset_by_family(trs_dataset, path_to_mathematica_folder, min_borders, max_borders, families):
    ''' Create mathematica dataset for visualization.
    TODO: check it.
    Create files with data for 3D plot groups with list_of_sat_dna. 
    And print matematica code for them.
    '''
    # clean folder
    default_file = os.path.join(path_to_mathematica_folder, "unseen.csv")
    default_file = default_file.replace("\\", "\\\\")
    path2fam = {}
    if os.path.isfile(default_file):
        os.remove(default_file)
    for trf_obj in trs_dataset:
        # family_network = get_family_name(trf_obj)
        # if trf_obj.trf_family_self in TRANSLATE_FAMS:
        #     trf_obj.trf_family_self = TRANSLATE_FAMS[trf_obj.trf_family_self]
        if trf_obj.trf_family_self in families:
          idx = families.index(trf_obj.trf_family_self)
          if idx < len(COLORS_POS):
            family_network = COLORS_POS[idx]
          else:
            family_network = None  
        else:
          family_network = None
        if family_network in COLORS:
          path = os.path.join(path_to_mathematica_folder, "%s.csv" % family_network)
          path = path.replace("\\", "\\\\")
          path2fam[family_network] = path
          if os.path.isfile(path):
              os.remove(path)
          with open(path, "w") as fh:
              # add mnimal and maximal values
              fh.write("%s\n%s\n" % (min_borders, max_borders))

    families_list = set()
    
    for trf_obj in trs_dataset:
        monomer = trf_obj.trf_period
        pmatch = trf_obj.trf_pmatch
        array_gc = trf_obj.trf_array_gc
        # family_network = get_family_name(trf_obj)
        # if trf_obj.trf_family_self in TRANSLATE_FAMS:
        #     trf_obj.trf_family_self = TRANSLATE_FAMS[trf_obj.trf_family_self]
        if trf_obj.trf_family_self in families:
          idx = families.index(trf_obj.trf_family_self)
          if idx < len(COLORS_POS):
            family_network = COLORS_POS[idx]
          else:
            family_network = None  
        else:
          family_network = None
        l = 3
        string = "%s,%s,%s,%s\n" % (monomer, pmatch, array_gc, l)
        families_list.add(family_network)
        
        with open(default_file, "a") as fd:
            if family_network in COLORS:
                with open(path2fam[family_network], "a") as fw:
                    fw.write(string)
            else:
                    fd.write(string)

    c = get_colors("OTHER")
    leg = Legend()
    leg.file = default_file
    leg.title = "Other"
    leg.color_r = c[0]
    leg.color_g = c[1]
    leg.color_b = c[2]
    legends = [leg]
    
    families_list = list(families_list)
    families_list.sort()

    for family in families_list:
        if family in path2fam:
          path = path2fam[family]
          c = get_colors(family)
          leg = Legend()
          leg.file = path
          idx = COLORS_POS.index(family)
          leg.title = families[idx]
          leg.color_r = c[0]
          leg.color_g = c[1]
          leg.color_b = c[2]
          legends.append(leg)
        else:
          path = default_file
          family = "OTHER"

    math_string = BUBBLE_CHART_TEMPLATE % join_legends(legends)
    print math_string
    return math_string
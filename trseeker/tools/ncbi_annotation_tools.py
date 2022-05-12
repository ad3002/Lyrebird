#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 26.09.2010
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com 
'''
Readers for various NCBI data.
'''

def get_ideogram_dict(idiogram_file, mode="NCBI"):
    ''' Read ideogram file and return return dict chr -> list: 

    TODO: add model.

    - 0 #chr
    - 1 arm
    - 2 band
    - 3 iscn_top
    - 4 iscn_bottom
    - 5 bases_top
    - 6 bases_bot
    - 7 stain
    - 8 density
    - 9 length, "length" = bases_top - bases_bot
    - 10 band letter: A
    - 11 band index: 1.1
    - 12 band short index: 1
    '''
    chrs = {}
    with open(idiogram_file, "r") as fh:
        for line in fh:
            if mode == "NCBI":
                if line.startswith("#chr"):
                    continue
                data = line.strip().split("\t")
                if int(data[5]) < int(data[6]):
                    data[5] = int(data[5])
                    data[6] = int(data[6])
                else:
                    data[5] = int(data[6])
                    data[6] = int(data[5])
                data[0] = data[0].upper()

                if len(data) == 8:
                    data.extend([0, "", ""])
                data.extend(["", "", "", ""])

                data[10] = data[2][:1]
                if len(data[2]) > 1:
                    data[11] = data[2][1:]
                    data[12] = data[11].split(".")[0]

                data.append(int(data[6]) - int(data[5]))
                if data[0] in chrs:
                    chrs[data[0]].append(data)
                else:
                    chrs[data[0]] = [data]
            else:
                '''
                keys: 
                    0 #chr
                    1 arm
                    2 band
                    3 iscn_top
                    4 iscn_bottom
                    5 bases_top
                    6 bases_bot
                    7 stain
                    8 density
                    9 length, "length" = bases_top - bases_bot
                    10 band letter: A
                    11 band index: 1.1
                    12 band short index: 1
                '''
                if line.startswith("#chr"):
                    continue
                data = line.strip().split("\t")

                if len(data) == 5:
                    data.extend([0, 0, 0, 0, "", "", "", ""])

                if int(data[1]) < int(data[2]):
                    data[5] = int(data[1])
                    data[6] = int(data[2])
                else:
                    data[5] = int(data[2])
                    data[6] = int(data[1])
                data[0] = data[0].replace("chr", "").upper()

                data[10] = data[3]
                data[11] = data[3]
                data[12] = data[3]

                data[12] = int(data[2]) - int(data[1])

                if data[4] == "gneg":
                    data[8] = 0
                elif data[4] == "acen":
                    data[8] = 100
                elif data[4] == "gpos100":
                    data[8] = 100
                elif data[4] == "gpos50":
                    data[8] = 50
                elif data[4] == "gpos75":
                    data[8] = 75
                elif data[4] == "gpos25":
                    data[8] = 25

                if data[0] in chrs:
                    chrs[data[0]].append(data)
                else:
                    chrs[data[0]] = [data]

    return chrs

def get_gene_list(file_gene_list, color='#000000', left_padding=30, gene_group_label='C57BL/6J'):
    ''' Read annotation file and returns dictionary:
        
    TODO: implement this.

    **chr->[ [chr_start, chr_end, gene_name, color, left_padding], ]**

    '''
    raise NotImplemented


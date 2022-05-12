#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 14.02.2014
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com 
"""
Functions related to repeats annotation in genomic data.
"""

import os
import re
from trseeker.seqio.tab_file import sc_iter_tab_file
from trseeker.seqio.tab_file import sc_write_model_to_tab_file
from trseeker.models.annotation_models import BEDModel, GFF3Model
from PyExp import trseeker_logger, runner
from trseeker.settings import load_settings
from trseeker.tools.repbase_tools import compute_repbase_stats

settings = load_settings()


def run_windowmasker(input_fasta, output_dir):
    ''' Wrapper for running WindowMasker.

    windowmasker [-h] [-help] [-xmlhelp] [-ustat unit_counts]
    [-in input_file_name] [-out output_file_name] [-checkdup check_duplicates]
    [-fa_list input_is_a_list] [-mem available_memory] [-unit unit_length]
    [-genome_size genome_size] [-window window_size] [-t_extend T_extend]
    [-t_thres T_threshold] [-set_t_high score_value] [-set_t_low score_value]
    [-parse_seqids] [-outfmt output_format] [-t_high T_high] [-t_low T_low]
    [-infmt input_format] [-exclude_ids exclude_id_list] [-ids id_list]
    [-text_match text_match_ids] [-use_ba use_bit_array_optimization]
    [-sformat unit_counts_format] [-smem available_memory] [-dust use_dust]
    [-dust_level dust_level] [-mk_counts] [-convert] [-version-full]
    '''
    data = {
        "input_fasta": input_fasta,
        "output_mk_counts": os.path.join(output_dir, "windowmasker.mk_counts"),
        "output_wm": os.path.join(output_dir, "windowmasker.txt"),
        "output_bed": os.path.join(output_dir, "windowmasker.bed"),
        "output_gff3": os.path.join(output_dir, "windowmasker.gff3"),
    }

    if not os.path.isfile(data["output_mk_counts"]):
        command = "windowmasker -mk_counts -in %(input_fasta)s -out %(output_mk_counts)s" % data
        runner.run(command)

    if not os.path.isfile(data["output_wm"]):
        command = "windowmasker -in %(input_fasta)s -ustat %(output_mk_counts)s -outfmt interval -out %(output_wm)s" % data
        runner.run(command)
    if not os.path.isfile(data["output_bed"]):
        command = """awk 'BEGIN { OFS="\t"; scaffold="" } { if (match($0, ">.*")) { gsub(">", "", $0); gsub(" .*", "", $0); scaffold=$0 } else { gsub(" - ", "\t", $0); print scaffold "\t" $1 "\t" $2+1 } }' %(output_wm)s> %(output_bed)s""" % data
        runner.run(command)

    # bedtools maskfasta -fi seabass.fasta -bed seabass.bed -fo seabass.masked.fa

    gff_objs = []
    kwargs = {
        "source": "WindowMasker",
        "ftype": "Repeat",
        "score": 0,
        "strand": "+",
        "phase": ".",
    }

    trseeker_logger.info("Convert to GFF3...")
    masked_bp = 0
    for bed_obj in sc_iter_tab_file(data["output_bed"], BEDModel):
        gff_obj = GFF3Model()
        gff_obj.set_with_bed_obj(bed_obj, **kwargs)
        gff_objs.append(gff_obj)
        masked_bp += abs(bed_obj.start - bed_obj.end)
    sc_write_model_to_tab_file(data["output_gff3"], gff_objs)
    return data["output_bed"], data["output_gff3"], masked_bp


def run_dustmasker(input_file, output_dir):
    ''' Wrapper for running dust.
    '''
    data = {
        'input_file': input_file,
        'output_file': os.path.join(output_dir, "dustmasker.txt"),
        'output_bed': os.path.join(output_dir, "dustmasker.bed"),
        'output_gff3': os.path.join(output_dir, "dustmasker.gff3"),
    }
    command = "dustmasker -in %(input_file)s -out %(output_file)s" % data
    runner.run(command)
    command = """awk 'BEGIN { OFS="\t"; scaffold="" } { if (match($0, ">.*")) { gsub(">", "", $0); gsub(" .*", "", $0); scaffold=$0 } else { gsub(" - ", "\t", $0); print scaffold "\t" $1 "\t" $2+1 } }' %(output_file)s > %(output_bed)s""" % data
    runner.run(command)

    gff_objs = []
    kwargs = {
        "source": "DUST",
        "ftype": "SimpleSequence",
        "score": 0,
        "strand": "+",
        "phase": "."
    }
    trseeker_logger.info("Convert to GFF3...")
    masked_bp = 0
    for bed_obj in sc_iter_tab_file(data["output_bed"], BEDModel):
        gff_obj = GFF3Model()
        gff_obj.set_with_bed_obj(bed_obj, **kwargs)
        gff_objs.append(gff_obj)
        masked_bp += abs(bed_obj.start - bed_obj.end)
    sc_write_model_to_tab_file(data["output_gff3"], gff_objs)
    return data["output_bed"], data["output_gff3"], masked_bp


def run_repeatmasker_repbase(input_fastas, output_dir, species, threads=20):
    ''' Wrapper for running RepeatMasker with Repbase database.
    TODO: remove hardcoded RM path.
    '''
    trseeker_logger.info("Running repeatmasker for %s library" % species)
    data = {
        "threads": threads,
        "species": species,
        "output_dir": output_dir,
        "input_fastas": input_fastas,
        "path": settings["repeatmasker"]["path"],
    }
    if os.listdir(output_dir): 
        d = raw_input("Folder %s in not empty (skip)" % output_dir) or None
        if not d:
            print "Skipped"
            return "Skipped"
    command = "%(path)s -s -pa %(threads)s -a -inv -xsmall -gff  -species %(species)s -dir %(output_dir)s  -nolow -source -html -u %(input_fastas)s" % data
    runner.run(command)
    return "ok"


def run_repeatmasker_local(input_fastas, output_dir, library_fasta):
    ''' Wrapper for running RepeatMasker with given database.
    TODO: remove hardcoded RM path.
    '''
    data = {
        "threads": 20,
        "library_fasta": library_fasta,
        "output_dir": output_dir,
        "input_fastas": input_fastas,
        "path": settings["repeatmasker"]["path"],
    }
    command = "%(path)s -s -pa %(threads)s -a -inv -xsmall -gff  -lib %(library_fasta)s -dir %(output_dir)s  -nolow -source -html -u %(input_fastas)s" % data
    runner.run(command)


def run_repeatscout(input_fasta, output_dir):
    '''
    '''
    file_name = os.path.split(input_fasta)[-1]
    result = {}
    data = {
        "input_fasta": input_fasta,
        "output_freq": os.path.join(output_dir, "scout.freq"),
        "library_file": os.path.join(output_dir, "scout.lib.fa"),
        "filtered_library": os.path.join(output_dir, "scout.lib.f1.fa"),
        "repeat_masker_output_dir": os.path.join(output_dir, "repeat_masker"),
        "repeat_masker_not_trf_output_dir": os.path.join(output_dir, "repeat_masker_no_trf"),
        "repeatmasker_tbl_file": os.path.join(output_dir, "repeat_masker", "%s.tbl" % file_name),
        "repeatmasker_ori_out_file": os.path.join(output_dir, "repeat_masker", "%s.ori.out" % file_name),
        "repeatmasker_gff_file": os.path.join(output_dir, "repeat_masker", "%s.out.gff" % file_name),
        "repeatmasker_masked_file": os.path.join(output_dir, "repeat_masker", "%s.masked" % file_name),
    }

    trseeker_logger.info("RepeatScout step 1: build kmer table")    
    if not os.path.isfile(data["output_freq"]):
        command = "/home/akomissarov/libs/RepeatScout-1/build_lmer_table -sequence %(input_fasta)s -freq %(output_freq)s -v" % data
        runner.run(command)
        
    trseeker_logger.info("RepeatScout step 2: build repeats families")
    if not os.path.isfile(data["library_file"]):
        command = "/home/akomissarov/libs/RepeatScout-1/RepeatScout -sequence %(input_fasta)s -output %(library_file)s -freq %(output_freq)s -vvvv" % data
        runner.run(command)

    trseeker_logger.info("RepeatScout step 3.1: run RepeatMasker")
    if not os.path.isdir(data["repeat_masker_output_dir"]):
        os.makedirs(data["repeat_masker_output_dir"])
        run_repeatmasker_local(data["input_fasta"], data["repeat_masker_output_dir"], data["library_file"])

    trseeker_logger.info("RepeatScout step 3.2: parse RepeatMasker results")
    if os.path.isfile(data["repeatmasker_tbl_file"]):
        with open(data["repeatmasker_tbl_file"]) as fh:
            content = fh.read()
            masked_percent = re.findall("bases masked:\s*(\d+)\s*bp\s*\(\s*(\d+.\d+)\s*\%\)", content, re.S)
            data["masked_bp"] = int(masked_percent[0][0])
            data["pmasked"] = float(masked_percent[0][1])
    else:
        data["masked_bp"] = 0
        data["pmasked"] = 0.0

    return data








        # if not os.path.isfile(data["filtered_library"]):
    #     command = "cat %(library_file)s |  /home/akomissarov/libs/RepeatScout-1/filter-stage-1.prl > %(filtered_library)s" % data
    #     print command
    #     os.system(command)



    # command = "cat %(filtered_library)s |  /home/akomissarov/libs/RepeatScout-1/filter-stage-2.prl --cat %(repeatmasker_out_file)s --thresh=%(rm_thresh)s" % data
    # print command
    # os.system(command)  

    # compile families
    

    # compare with known tandem repeats by blast

    # compare with known repeats by blast
    # for family_fasta_file in iter_folder()

    







def run_ltr_harvest(input_fasta, output_dir):
    ''' Run LTR_Harvest and mkvtree clusterization.
    '''
    file_name = os.path.split(input_fasta)[-1]
    data = {
        "input_fasta": input_fasta,
        "index_name": os.path.join(output_dir, "indexname.fa"),
        "pred_out_file": os.path.join(output_dir, "ltrharvest.out.fa"),
        "gff3_out_file": os.path.join(output_dir, "ltrharvest.out.gff3"),
        "out_out_file": os.path.join(output_dir, "ltrharvest.out.out"),
        "pred_out_inner_file": os.path.join(output_dir, "ltrharvest.out.inner.fa"),
        "lic_file": os.path.join(output_dir, "vmatch.lic"),
        "repeat_masker_output_dir": os.path.join(output_dir, "repeatmasker"),
        "repeatmasker_out_file": os.path.join(output_dir, "repeatmasker", "%s.tbl" % file_name),
    }
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    if not os.path.isfile(data["pred_out_file"]):
        command = "/home/akomissarov/libs/genometools-unstable/bin/gt suffixerator -db %(input_fasta)s -indexname %(index_name)s -tis -suf -lcp -des -ssp -sds -dna" % data
        print command
        os.system(command)
        # command = "/home/akomissarov/libs/genometools-unstable/bin/gt ltrharvest -index %(index_name)s -v -out %(pred_out_file)s -outinner %(pred_out_inner_file)s -seed 30 -xdrop 5 -mat 2 -mis -2 -ins -3 -del -3 -minlenltr 100 -maxlenltr 1000 -mindistltr 1000 -maxdistltr 15000 -similar 90.0 -overlaps all -mintsd 5 -maxtsd 20 -motif tgca -motifmis 0 -vic 60 -longoutput" % data
        command = "/home/akomissarov/libs/genometools-unstable/bin/gt ltrharvest -index %(index_name)s -v -out %(pred_out_file)s -outinner %(pred_out_inner_file)s -gff3 %(gff3_out_file)s -minlenltr 100 -maxlenltr 6000 -mindistltr 1500 -maxdistltr 25000 -mintsd 5 -maxtsd 5 -motif tgca -similar 99 -vic 10  > %(out_out_file)s" % data
        print command
        os.system(command)

    # if not os.path.isfile(data["lic_file"]):
    #     command = "cp /home/akomissarov/libs/vmatch.lic %s" % output_dir
    #     print command
    #     os.system(command)

    #     command = "/home/akomissarov/libs/vmatch-2.2.3-Linux_x86_64-64bit/mkvtree -db %(pred_out_file)s -dna -pl -allout -v" % data
    #     print command
    #     os.system(command)

    #     command = "/home/akomissarov/libs/vmatch-2.2.3-Linux_x86_64-64bit/vmatch -dbcluster 95 7 Cluster-pred-chrAll -p -d -seedlength 50 -l 1101 -exdrop 9 %(pred_out_file)s" % data
    #     print command
    #     os.system(command)

    if not os.path.isdir(data["repeat_masker_output_dir"]):
        os.makedirs(data["repeat_masker_output_dir"])
        run_repeatmasker_local(data["input_fasta"], data["repeat_masker_output_dir"], data["pred_out_file"])

    print data["repeatmasker_out_file"]
    if os.path.isfile(data["repeatmasker_out_file"]):
        with open(data["repeatmasker_out_file"]) as fh:
            content = fh.read()
            masked_percent = re.findall("bases masked:\s*(\d+)\s*bp\s*\(\s*(\d+.\d+)\s*\%\)", content, re.S)
            data["masked"] = int(masked_percent[0][0])
            data["pmasked"] = float(masked_percent[0][1])
    else:
        data["masked"] = 0
        data["pmasked"] = 0.0

    return data










def run_process_repeatscout_results(input_fasta, output_dir):
    '''
    '''
    only_file_name = os.path.split(input_fasta)[-1]
    data = {
        "repeatmasker_ori_out_file": os.path.join(output_dir, "repeat_masker", "%s.ori.out" % only_file_name),
        "output_families_folder": os.path.join(output_dir, "families"),
        "output_tab_file": os.path.join(output_dir, "scout.families.tsv"),
    }
    if not os.path.isdir(data["output_families_folder"]):
        os.makedirs(data["output_families_folder"])

    
    compute_repbase_stats(data["repeatmasker_ori_out_file"], data["output_tab_file"], input_fasta, data["output_families_folder"])
    

def run_repeat_modeler(db_name, input_file, output_dir):
    '''
    '''
    data = {
        'db_name': db_name,
        'input_file': input_file,
    }
    command = "/home/akomissarov/libs/RepeatModeler/BuildDatabase -name %(db_name)s -engine ncbi %(input_file)s" % data
    print command
    os.system(command)
    command = "/home/akomissarov/libs/RepeatModeler/RepeatModeler -database %(db_name)s -engine ncbi -pa 30" % data
    print command
    print os.system(command)







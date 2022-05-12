#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 07.09.2010
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com
"""
BLAST wrapper.
Used settings

>>> location = settings["blast_settings"]["blast_location"]
>>> b = settings["blast_settings"]["blast_b"]
>>> e = settings["blast_settings"]["blast_e"]
>>> v = settings["blast_settings"]["blast_v"]
>>> OS = settings["trseeker"]["os"]

"""
import sys
from trseeker.seqio.fasta_file import sc_iter_fasta

sys.path.append("/Users/akomissarov/Dropbox/workspace")

from collections import defaultdict

from trseeker.seqio.tab_file import sc_iter_tab_file, sc_iter_simple_tab_file
from trseeker.models.blast_model import BlastResultModel
from trseeker.models.trf_model import TRModel
import os

from trseeker.settings import load_settings
from trseeker.seqio.tr_file import get_all_trf_objs

settings = load_settings()
location = settings["blast_settings"]["blast_location"]
if not location:
    location = ""
b = settings["blast_settings"]["blast_b"]
e = settings["blast_settings"]["blast_e"]
v = settings["blast_settings"]["blast_v"]
OS = settings["trseeker"]["os"]


def blastn(database, query, output, e_value=None, dust=True, nhits=None, threads=30):
    """ Blast fasta file versus blast database.
    Available output format parameters:
        7 qseqid qgi qacc sseqid means sgi sacc qstart qend sstart send evalue bitscore score length pident nident mismatch positive gapopen gaps ppos frames qframe sframe
    Command examples:
        blastn -query HTT_gene -task megablast -db hs_chr -db_soft_mask 30 -outfmt 7 -out HTT_megablast_mask.out -num_threads 4
        $ echo 1786181 | ./blastn -db ecoli -outfmt "7 qacc sacc evalue qstart qend sstart send"
    @todo: do something with parameters
    @param database: database created with makeblastdb
    @param query: fasta file
    @param output: output file
    @param e_value: e value, e.g. 1e-20
    @return: None
    """
    if OS == "win":
        program_name = "blastn.exe"
    else:
        program_name = "blastn"
    format_string = '   '

    if not e_value:
        e_value = e

    data = {
        'location': location,
        'program_name': program_name,
        'query': query,
        'database': database,
        'output': output,
        'e_value': e_value,
        'format_string': format_string,
        'b': b,
        'v': v,
        'threads': threads,
    }

    if nhits:
        data['b'] = nhits

    if dust:
        data["dust"] = "-dust yes"
    else:
        data["dust"] = "-dust no"

    string = '%(location)s%(program_name)s -query %(query)s -task blastn -db %(database)s -out %(output)s -evalue %(e_value)s -word_size 10 -outfmt %(format_string)s %(dust)s -soft_masking "false" -max_target_seqs %(b)s -num_threads %(threads)s' % data

    print string
    os.system(string)


def blastp(database, query, output, e_value=None, threads=30, nhits=50):
    """ Blastp
    """
    if OS == "win":
        program_name = "blastp.exe"
    else:
        program_name = "blastp"
    format_string = '"7 qseqid qgi qacc sseqid means sgi sacc qstart qend sstart send evalue bitscore score length pident nident mismatch positive gapopen gaps ppos frames qframe sframe"'

    if not e_value:
        e_value = e

    data = {
        'location': location,
        'program_name': program_name,
        'query': query,
        'database': database,
        'output': output,
        'e_value': e_value,
        'format_string': format_string,
        'b': nhits,
        'v': v,
        'threads': threads,
    }

    string = "%(location)s%(program_name)s -query %(query)s -db %(database)s -out %(output)s -evalue %(e_value)s -outfmt %(format_string)s -max_target_seqs %(b)s -num_threads %(threads)s" % data

    print string
    os.system(string)


def create_db(fasta_file, output, verbose=False, title=None, prot=False):
    """  Create BLAST database.
    Command example:
        "%s%s -in %s -hash_index -dbtype nucl -parse_seqids -out %s" % (location, progr_name, fasta_file, output)
        formatdb.exe -i E:\home\ad3002\work\mouse_wgs\fa\name.fa -p F -o T -V T -n mouse_wgs -t mouse_wgs
        makeblastdb -in hs_chr.fa    hs_chr -title "Human chromosomes, Ref B37.1"
    @todo: do something with parameters
    @param fasta_file: fasta file
    @param output: output database path
    @param verbose: verbose
    @param title: database title
    @return: None
    """
    if verbose:
        print "Format fasta file for BLAST search ..."

    if OS == "win":
        program_name = "makeblastdb.exe"
    else:
        program_name = "makeblastdb"
    if not prot:
        string = "%s%s -in %s -hash_index -dbtype nucl -parse_seqids -out %s" % (location, program_name, fasta_file, output)
    else:
        string = "%s%s -in %s -hash_index -dbtype prot -parse_seqids -out %s" % (location, program_name, fasta_file, output)
    if title:
        string += " -title %s" % title
    if verbose:
        print string
    os.system(string)
    if verbose:
        print "Format fasta file for BLAST search complete!"


def alias_tool(dblist, output, title):
    """ Create alias database
    Command example:
        blastdb_aliastool -dblist "nematode_mrna nematode_genomic" -dbtype nucl -out nematode_all -title "Nematode RefSeq mRNA + Genomic"
    @param dblist: list of database paths
    @param output: output database path
    @param title: alias database title
    @return: None
    """
    if OS == "win":
        program_name = "blastdb_aliastool.exe"
    else:
        program_name = "blastdb_aliastool"
    string = "%s%s -dblist %s -dbtype nucl -out %s -title %s" % (location, program_name, dblist, output, title)
    print string
    os.system(string)


def bl2seq(input1, input2, output):
    """
    Compare two sequences with blast.
    @todo: do something with parameters
    @param input1: file with fasta
    @param input2: file with fasta
    @param output: output file
    @return: None
    """
    if OS == "win":
        _location = settings["blast_settings"]["blast_location_WIN"]
    else:
        _location = settings["blast_settings"]["blast_location_NIX"]

    string = "%s -p blastn -i %s -j %s -o %s -F F -W 5 -D 1 -e %s" % (
                                                        _location,
                                                        input1,
                                                        input2,
                                                        output,
                                                        e)
    print string
    os.system(string)


def create_db_for_genome(file_pattern=None, chromosome_list=None, output=None, title=None):
    """
    Create BLAST database for genome.
    Example:
        create_db_for_genome(file_pattern="E:/eukaryota_genomes/mus_musculus/fasta/%s.fsa_nt",
                     chromosome_list=["wgs.AAHY.1", "wgs.AAHY.2", "wgs.AAHY.3", "wgs.AAHY.4", "wgs.AAHY.5", "wgs.AAHY.6", "wgs.CAAA"],
                     output="mus_musculus_wgsa",
                     title="Mouse WGSAs")
    @param file_pattern: path pattern
    @param chromosome_list: list of chromosomes for file pattern
    @param output: final database path
    @param title: final database title
    @return: None
    """
    files = []
    for x in chromosome_list:
        fasta_file = file_pattern % x
        print "Database construction: ", fasta_file
        create_db(fasta_file, fasta_file, verbose=True, title=fasta_file.split("/")[-1])
        files.append(fasta_file)
    dblist = '"' + " ".join(files) + '"'
    title = '"%s"' % title
    alias_tool(dblist, output, title)


def filter_raw_blast(blast_file, fitlers):
    """ Clear blast file.
    """
    pass


def get_gi_list(gi_score_file, score_limit=90):
    """
    Get GI list from given score file.
    @param gi_score_file: GI score file
    @param score_limit: get only GI with score greater than score limit (default 90)
    @return: list of GI
    """
    result = []
    with open(gi_score_file, "rb") as fh:
        for line in fh:
            data = line.strip().split("\t")
            if float(data[1]) > score_limit:
                result.append(data[0])
    return result


def get_all_gi_from_blast(blast_file, mode="gi"):
    """
    Get all GI from blast results file.
    @param blast_file: blast file
    @param mode: parsing mode GI or ref based
    @raise: Unknown mode
    @return: a dictionary of gi or ref to number of hits
    """
    result = {}
    for blast_obj in sc_iter_tab_file(blast_file, BlastResultModel, remove_starts_with="#"):
        if mode == "gi":
            _id = blast_obj.subject_gi
        elif mode == "ref":
            _id = blast_obj.subject_ref
        else:
            raise Exception("Unknown mode")
        result.setdefault(_id, 0)
        result[_id] += 1

    result = [(_id, result[_id]) for _id in result.keys()]
    result.sort(reverse=True, key=lambda x: x[1])
    return result


def get_all_blast_obj_from_blast(blast_file, mode="ref"):
    """
    Get all GI to blast_obj list dictionary from blast output.
    @param blast_file: blast file
    @param mode: parsing mode (gi or ref)
    @return: dictionary GI to blast_obj
    """
    result = {}
    for blast_obj in sc_iter_tab_file(blast_file, BlastResultModel, remove_starts_with="#"):
        if mode == "gi":
            _id = blast_obj.subject_gi
        elif mode == "ref":
            _id = blast_obj.subject_ref
        else:
            raise Exception
        result.setdefault(_id, [])
        result[_id].append(blast_obj)
    return result


def get_all_blast_obj_from_multiblast(blast_file):
    """
    Get a dictionary query ref to a list of blast objects with hits
    @param blast_file: blast file
    @return: dictionary GI to blast_obj
    """
    result = {}
    for blast_obj in sc_iter_tab_file(blast_file, BlastResultModel, remove_starts_with="#"):
        result.setdefault(blast_obj.query_ref, [])
        result[blast_obj.query_ref].append(blast_obj)
    return result


def bl2seq_search_for_trs(trf_large_file, annotation_bl2seq_folder, temp_file):
    """
    Pairwise search of TRs with bl2seq.
    @param trf_large_file: file with TRs
    @param annotation_bl2seq_folder: folder for output files
    @param temp_file: temp file name
    @return: None
    """
    print "Read TRs..."
    trf_objs = get_all_trf_objs(trf_large_file)
    N = len(trf_objs)
    for i in xrange(N):
        for j in xrange(i, N):
            print "Searching:", i, j, "from", N
            output_file = "%s.%s.blast" % (trf_objs[i].trf_id, trf_objs[j].trf_id)
            blast_output_file = os.path.join(annotation_bl2seq_folder, output_file)
            input1 = trf_objs[i].trf_array
            input2 = trf_objs[j].trf_array
            input_file1 = temp_file + ".1.fasta"
            input_file2 = temp_file + ".2.fasta"
            with open(input_file1, "w") as fh:
                fh.write(">%s\n%s" % (trf_objs[i].trf_id, input1))
            with open(input_file2, "w") as fh:
                fh.write(">%s\n%s" % (trf_objs[j].trf_id, input2))
            bl2seq(input_file1, input_file2, blast_output_file)


def blastn_search_for_trs(trf_large_file, db, annotation_self_folder, temp_file, skip_by_family=None, is_huge_alpha=False, skip_short=None):
    """
    Search TRs in given DB.
    @param trf_large_file: TRs file
    @param db: database
    @param annotation_self_folder: folder for output files
    @param temp_file: temp file
    @param skip_by_family: skip TRs by family name
    @param is_huge_alpha: skip ALPHA families
    @return: None
    """
    for n, u in enumerate(sc_iter_tab_file(trf_large_file, TRModel)):
        pass
    
    alpha_sets = {}
    for i, u in enumerate(sc_iter_tab_file(trf_large_file, TRModel)):
        
        if skip_short and u.trf_array_length < skip_short:
            print "Skipped", i
            continue

        print "%s/%s" % (i, n)
        if skip_by_family:
            if u.trf_family in skip_by_family:
                continue
        if u.trf_family == "ALPHA":
            continue
        blast_output_file = os.path.join(annotation_self_folder, "%s.blast" % u.trf_id)
        if os.path.isfile(blast_output_file):
            with open(blast_output_file) as fh:
                    data = fh.read()
            if data.startswith("ALPHA"):
                continue
        seen = False
        if is_huge_alpha:
            for key, tr_set in alpha_sets.items():
                if u.trf_id in tr_set:
                    with open(blast_output_file, "w") as fh:
                        fh.write("ALPHA\t%s" % key)
                    seen = True
        if seen:
            continue

        if os.path.isfile(blast_output_file) and os.stat(blast_output_file).st_size != 0:
            trids = _get_trids_from_blast_file(blast_output_file)
            alpha_sets[u.trf_id] = trids
            print "ADDED ALPHA by %s with length %s" % (u.trf_id, len(trids)) 
            continue

        with open(temp_file, "w") as fh:
            fh.write(">%s\n%s" % (u.trf_id, u.trf_array))
        blastn(db, temp_file, blast_output_file)
        if os.path.isfile(temp_file):
            os.remove(temp_file)

        if is_huge_alpha:
            trids = _get_trids_from_blast_file(blast_output_file)
            alpha_sets[u.trf_id] = trids
            print "ADDED ALPHA by %s with length %s" % (u.trf_id, len(trids))


def _get_trids_from_blast_file(blast_output_file):
    """
    Get trids from blast output file.
    @param blast_output_file: blast output file
    @return: return list of trids
    """
    result = set()
    gi2obj = get_all_blast_obj_from_blast(blast_output_file)
    for gi, hits in gi2obj.items():
        for hit in hits:
            if hit.score >= 90:
                result.add(int(hit.subject_ref))
    return list(result)



# if __name__ == '__main__':
    

#     assembly = "/Users/akomissarov/Dropbox/ARGolden.scafSeq.kgf.gap1_10.filter"
#     assembly_len_file = "/Users/akomissarov/Dropbox/ARGolden.scafSeq.kgf.gap1_10.filter.length"
#     markers_file = "/Users/akomissarov/Dropbox/workspace/PySatDNA/arowana_lg.tcv"
#     markers_fasta_file = "/Users/akomissarov/Dropbox/workspace/PySatDNA/arowana_lg.fa"
#     markers_length_file = "/Users/akomissarov/Dropbox/workspace/PySatDNA/arowana_lg.legnth"
#     query_len_file = "/Users/akomissarov/Dropbox/workspace/PySatDNA/arowana_lg.tcv.obj"
#     blast_file = "/Users/akomissarov/Dropbox/workspace/PySatDNA/ARGolden.scafSeq.kgf.gap1_10.filter.blast.lg.blast"
#     parsed_blast_file = "/Users/akomissarov/Dropbox/workspace/PySatDNA/ARGolden.scafSeq.kgf.gap1_10.filter.blast.lg.blast.pickle"
    
#     from trseeker.seqio.fasta_file import sc_iter_fasta
#     from trseeker.models.sequence_model import SequenceModel, MarkerModel
#     import pickle
#     from collections import Counter

#     print "Load makers"
#     if not os.path.isfile(query_len_file):
#         markers = {}      
#         for i, marker_obj in enumerate(sc_iter_tab_file(markers_file, MarkerModel)):
#             markers[marker_obj.genbank_id] = marker_obj
#             print "Read ", i, marker_obj.genbank_id
#         with open(query_len_file, "wb") as fh:
#             pickle.dump(markers, fh)    
#     else:
#         with open(query_len_file) as fh:
#             markers = pickle.load(fh)

            

#     print "Load markers lengths"
#     if not os.path.isfile(markers_length_file):        
#         markers_lengths = {}
#         for seq_obj in sc_iter_fasta(markers_fasta_file):
#             ref = seq_obj.header.split("|")[-2].split(".")[0]
#             markers_lengths[ref] = seq_obj.length
#             print ref, seq_obj.length
#         with open(markers_length_file, "wb") as fh:
#             pickle.dump(markers_lengths, fh)
#     else:
#         with open(markers_length_file) as fh:
#             markers_lengths = pickle.load(fh)

#     print "Load lengths"
#     if not os.path.isfile(assembly_len_file):        
#         lengths = {}
#         for seq_obj in sc_iter_fasta(assembly):
#             lengths[seq_obj.header] = seq_obj.length
#             print seq_obj.header, seq_obj.length
#         with open(assembly_len_file, "wb") as fh:
#             pickle.dump(lengths, fh)
#     else:
#         with open(assembly_len_file) as fh:
#             lengths = pickle.load(fh)

#     print "Read blast db"
#     if not os.path.isfile(parsed_blast_file):
#         blast_objs = get_all_blast_obj_from_multiblast(blast_file)
#         with open(parsed_blast_file, "wb") as fh:
#             pickle.dump(blast_objs, fh)
#     else:
#         with open(parsed_blast_file) as fh:
#             blast_objs = pickle.load(fh)

#     print "Show results"
#     result = []
#     scaffold2lg = {}
#     for key, items in blast_objs.items():
#         for item in items:
#             ref = key.split("|")[-2].split(".")[0]
#             item.fraction_of_query = abs(item.query_start - item.query_end) / float(markers_lengths[ref])
#             print ref, item.fraction_of_query

#     #         scaffold2lg.setdefault(item.subject_ref, Counter())
#     #         scaffold2lg[item.subject_ref].update([markers[ref].lg])


            
#     #         if len(items) == 1:
#     #             result.append((key, ref, len(items), item.subject_ref, lengths[item.subject_ref], item.subject_start, item.subject_end, item.alignment_length, item.proc_identity, markers[ref].lg, markers[ref].cm_average.replace(",",".")))
#     # result.sort(key=lambda x: (int(x[-2]), float(x[-1])))
#     # # for r in result:
#     # #     print r
#     # for s in scaffold2lg:
#     #     print s, scaffold2lg[s]


if __name__ == '__main__':

    fasta_markers = "/home/akomissarov/Dropbox/PySatDNA/OtherFiles/arowana_lg.fa"

    ref2length = {}
    for seq_obj in sc_iter_fasta(fasta_markers):
        ref = seq_obj.header.split()[0]
        ref2length[ref] = seq_obj.length

    marker2hits = {}

    makers_lgs = "/home/akomissarov/Dropbox/PySatDNA/OtherFiles/arowana_lg.tcv"
    ref2lg = {}
    for data in sc_iter_simple_tab_file(makers_lgs):
        ref2lg[data[1]] = (data[2], data[3], data[4], data[5])
        #print ref2lg

        marker2hits[data[1]] = 0



    # blast_file = "/storage1/akomissarov/genome_arowana/assemblies/golden/ARGolden.scafSeq.kgf.gap1_10.filter.blast.lg.blast.hq90"
    # output_file = "/home/akomissarov/Dropbox/PySatDNA/OtherFiles/golden_hits.data"
    # md_file = "/home/akomissarov/Dropbox/PySatDNA/gold.microd.tsv"
    # md_lgs_file = "/home/akomissarov/Dropbox/PySatDNA/gold.microd.lgs.tsv"
    # lgs_file = "/home/akomissarov/Dropbox/PySatDNA/gold.microd.lgs.plus.tsv"

    # blast_file = "/storage1/akomissarov/genome_arowana/assemblies/red/ARRed.scafSeq.kgf.gap1_10.filter.blast.lg.blast"
    # output_file = "/home/akomissarov/Dropbox/PySatDNA/OtherFiles/red_hits.data"
    # md_file = "/home/akomissarov/Dropbox/PySatDNA/red.microd.tsv"
    # md_lgs_file = "/home/akomissarov/Dropbox/PySatDNA/red.microd.lgs.tsv"
    # lgs_file = "/home/akomissarov/Dropbox/PySatDNA/red.microd.lgs.plus.tsv"

    # # blast_file = "/storage1/akomissarov/genome_arowana/assemblies/green/ARGreen.scafSeq.kgf.gap1_10.filter.blast.lg.blast"
    # # output_file = "/home/akomissarov/Dropbox/PySatDNA/OtherFiles/green_hits.data"
    # # md_file = "/home/akomissarov/Dropbox/PySatDNA/green.microd.tsv"
    # # md_lgs_file = "/home/akomissarov/Dropbox/PySatDNA/green.microd.lgs.tsv"
    # # lgs_file = "/home/akomissarov/Dropbox/PySatDNA/green.microd.lgs.plus.tsv"

    # scaf2lgs = defaultdict(list)
    # for data in sc_iter_simple_tab_file(output_file):
    #     scaf = data[3]
    #     lg = data[7]
    #     scaf2lgs[scaf].append(lg)

    # from trseeker.tools.statistics import get_simple_statistics

    # result = {}

    # with open(md_lgs_file, "w") as fh:
    #     for data in sc_iter_simple_tab_file(md_file):
    #         if data[0] == "-":
    #             header =  "%s\n" % "\t".join(data)
    #             fh.write(header)
    #             continue
    #         scaffold, aroA2, aro1, aro3, aro4, aroA1 = data
    #         d = [(int(x), i) for i,x in enumerate([aroA2, aro1, aro3, aro4, aroA1])]
    #         lg = ",".join(set(scaf2lgs[data[0]]))
    #         max_cov = max(d)

    #         stats = get_simple_statistics(map(int,[aroA2, aro1, aro3, aro4, aroA1]))

    #         print scaffold, d, max_cov, "LG:", lg, "std:", stats["standard_deviation"]

    #         if not lg:
    #             lg = "unk"

    #         result.setdefault(lg, [0,0,0,0,0])
    #         result[lg][0] += int(aroA2)
    #         result[lg][1] += int(aro1)
    #         result[lg][2] += int(aro3)
    #         result[lg][3] += int(aro4)
    #         result[lg][4] += int(aroA1)

    #         data.append(lg)
    #         fh.write("%s\n" % "\t".join(map(str,data)))

    # with open(lgs_file, "w") as fh:
    #     lgs = result.keys()
    #     lgs.sort()

    #     fh.write(header)

    #     for lg in lgs:
    #         s = "%s\t%s\n" % (lg, "\t".join(map(str, result[lg])))
    #         fh.write(s)

    files = [
        "/home/akomissarov/Dropbox/PySatDNA/gold.microd.lgs.plus.tsv",
        "/home/akomissarov/Dropbox/PySatDNA/red.microd.lgs.plus.tsv",
        "/home/akomissarov/Dropbox/PySatDNA/green.microd.lgs.plus.tsv",
    ]
    

    header = "-\tgold_aroA2\tgold_aro1\tgold_aro3\tgold_aro4\tgold_aroA1\tred_aroA2\tred_aro1\tred_aro3\tred_aro4\tred_aroA1\tgreen_aroA2\tgreen_aro1\tgreen_aro3\tgreen_aro4\tgreen_aroA1\n"
    fh = open("/home/akomissarov/Dropbox/PySatDNA/all.microd.lgs.plus.tsv", "w")
    results = {}
    for i, file_name in enumerate(files):
        for data in sc_iter_simple_tab_file(file_name):
            if data[0] == "-":
                continue
            scaffold, aroA2, aro1, aro3, aro4, aroA1 = data
            results.setdefault(scaffold, [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
            results[scaffold][i*5+0] = aroA2
            results[scaffold][i*5+1] = aro1
            results[scaffold][i*5+2] = aro3
            results[scaffold][i*5+3] = aro4
            results[scaffold][i*5+4] = aroA1
    fh.write(header)

    lgs = results.keys()
    lgs.sort()
    for lg in lgs:
        s = "%s\t%s\n" % (lg, "\t".join(map(str, results[lg])))
        fh.write(s)


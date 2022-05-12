#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 26.12.2010
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com 
'''
Functions related to Repbase datasets.
'''
import os, re
from trseeker.tools.sequence_tools import clear_sequence
from PyExp import sc_iter_filepath_folder
from trseeker.seqio.fasta_file import sc_iter_fasta
from trseeker.seqio.tab_file import sc_iter_tab_file
from trseeker.models.repbase_model import RepbaseModel
from collections import Counter, defaultdict
from trseeker.tools.statistics import get_simple_statistics
from trseeker.tools.sequence_tools import get_revcomp

def _join_repbase_fastas(fasta_file, fasta_output, index_file_name=None, file_name=None, start_id=None, increment=False, repbase_names_fix=False):
    ''' Function join all Repbase fasta files in one huge fasta; reformat headers for compatibility with NCBI tools
    '''
    if not increment:
        if os.path.isfile(fasta_output):
            os.unlink(fasta_output)

    with open(fasta_output, "a") as fa_fw:
        if start_id:
            i = start_id
        else:
            i = 0
        for seq_obj in sc_iter_fasta(fasta_file):
            i += 1
            head = seq_obj.seq_head
            if file_name and repbase_names_fix:
                file_name = file_name.split(".")[0]
                head_info = "%s:%s" % (file_name, head.replace("  ", " ").replace(" ", "_").replace(">", "").replace("(", "").replace(")", "").replace("\t", ":").strip())
            else:
                head_info = head.strip().replace(">", "")

            fa_fw.write(">%s\n%s\n" % (head_info, seq_obj.sequence))
    return i

def join_repbase_files(input_folder, output_file):
    ''' Function join all Repbase fasta files in one huge fasta; reformat headers for compatibility with NCBI tools.'''
    start_id = 0
    for input_file in sc_iter_filepath_folder(input_folder, mask="."):
        if not input_file.endswith("ref"):
            continue
        start_id = _join_repbase_fastas(input_file,
                                   output_file,
                                   file_name=os.path.split(input_file)[-1],
                                   start_id=start_id,
                                   increment=True,
                                   repbase_names_fix=True)

def fix_repeatmasker_out_file(file_name):
    ''' Fix names in RepeatMasker file.
    '''
    def preprocess_function(x):
        x = re.sub("\s+","\t",x)
        x = x.strip()
        x = re.sub("[)()]","",x)
        x += "\n"
        return x
    with open(file_name, "r") as fh:
        data = fh.readlines()
    data = [preprocess_function(x) for x in data]
    with open(file_name, "w") as fh:
        fh.writelines(data)


def compute_repbase_stats(file_name, output_file, input_fasta, output_folder, length_cutoff=0.9):
    ''' Compute statistics for RepeatMasker output file.
    '''
    result = {}
    family2n = defaultdict(int)
    family2length = defaultdict(int)


    print "Read genome sequences"
    genome = {}
    for seq_obj in sc_iter_fasta(input_fasta):
        genome[seq_obj.header] = seq_obj.sequence


    fix_repeatmasker_out_file(file_name)

    for i, rm_object in enumerate(sc_iter_tab_file(file_name, RepbaseModel)):
        if i % 100000 == 0:
            print i
        rm_object.set_element_coverage()
        if not rm_object.name in result:
            result[rm_object.name] = {
                'nall': 0,
                'n0': 0,
                'n50': 0,
                'n80': 0,
                'n98': 0,
                'dall': [],
                'd0': [],
                'd50': [],
                'd80': [],
                'd98': [],
            }
        result[rm_object.name]['nall'] += 1
        if rm_object.pcov <= 0.5:
            result[rm_object.name]['n0'] += 1
            result[rm_object.name]['d0'].append(rm_object.pdivergence)
        if rm_object.pcov > 0.5:
            result[rm_object.name]['n50'] += 1
            result[rm_object.name]['d50'].append(rm_object.pdivergence)
        if rm_object.pcov > 0.8:
            result[rm_object.name]['n80'] += 1
            result[rm_object.name]['d80'].append(rm_object.pdivergence)
        if rm_object.pcov > 0.98:
            result[rm_object.name]['n98'] += 1
            result[rm_object.name]['d98'].append(rm_object.pdivergence)


        if rm_object.pcov > length_cutoff:
            sequence = rm_object.get_sequence(genome)
            if not sequence:
                print "Error. Not found.", rm_object
            else:
                if rm_object.query.endswith("|_1"):
                    continue
                name = rm_object.name.replace("/","__")
                family_file_name = os.path.join(output_folder, name)
                family_file_name_fa = os.path.join(output_folder, name+".fa")
                with open(family_file_name, "a") as fh:
                    s = str(rm_object).strip()
                    fh.write("%s\t%s\n" % (s, sequence))
                with open(family_file_name_fa, "a") as fh:
                    s = str(rm_object).strip()
                    if rm_object.strand == "C":
                        sequence = get_revcomp(sequence)
                    head = "%s_%s_%s_%s_%s" % (rm_object.query, rm_object.qstart, rm_object.qend, rm_object.name, rm_object.pdivergence)
                    fh.write(">%s\n%s\n" % (head, sequence))

                family2n[rm_object.name] += 1
                family2length[rm_object.name] += len(sequence)


    print
    print "Save dataset to", output_file
    with open(output_file, "w") as fh:
        s = "name\tn\tmean_length\tnall\tn0\td0\td0std\tn50\td50\td50std\tn80\td80std\tn98\td98\td98std\n"
        fh.write(s)
        for te in result:
            stat = get_simple_statistics(result[te]['d0'])
            result[te]['d0'] = stat["mean"]
            result[te]['d0std'] = stat["standard_deviation"]
            stat = get_simple_statistics(result[te]['d50'])
            result[te]['d50'] = stat["mean"]
            result[te]['d50std'] = stat["standard_deviation"]
            stat = get_simple_statistics(result[te]['d80'])
            result[te]['d80'] = stat["mean"]
            result[te]['d80std'] = stat["standard_deviation"]
            stat = get_simple_statistics(result[te]['d98'])
            result[te]['d98'] = stat["mean"]
            result[te]['d98std'] = stat["standard_deviation"]
            result[te]['name'] = te
            result[te]['n'] = family2n[te]
            result[te]['mean_length'] = family2length[te]
            if family2n[te]:
                result[te]['mean_length'] /= family2n[te]
            s = "%(name)s\t%(n)s\t%(mean_length)s\t%(nall)s\t%(n0)s\t%(d0)s\t%(d0std)s\t%(n50)s\t%(d50)s\t%(d50std)s\t%(n80)s\t%(d80std)s\t%(n98)s\t%(d98)s\t%(d98std)s\n" % result[te]
            fh.write(s)


if __name__ == '__main__':
    
    file_name = "/storage1/akomissarov/genome_leo/repeatmasker_all/ref_assembly/pantheraLeoScaffolds.fa.ori.out"
    output_file = "/home/akomissarov/Dropbox/PySatDNA/leo.result.tsv"
    file_name = "/storageSM2/akomissarov/genomes_mammals/pid_gmam_felis_catus/repeatmasker_all/ref_assembly/Felis_catus-6_2.fasta.ori.out"
    output_file = "/home/akomissarov/Dropbox/PySatDNA/felis.result.tsv"

    file_name = "/storage1/akomissarov/genome_nerpa/repeatmasker_all/ref_assembly/out_gapClosed.fa.ori.out"
    output_file = "/home/akomissarov/Dropbox/PySatDNA/nerpa.result.tsv"

    file_name = "/storage1/akomissarov/genome_pangolin/repeatmasker_all/ref_assembly/sga_scaff_min1k.fa.ori.out"
    output_file = "/home/akomissarov/Dropbox/PySatDNA/pangolin.result.tsv"    


    file_name = "/storage1/akomissarov/pacbio/repeatmasker/pacbio2_rna/L_RNA_scaffolder.fasta.ori.out"
    output_file = "/home/akomissarov/Dropbox/PySatDNA/seabass.result.tsv"    

    file_name = "/storageSM2/akomissarov/genomes_mammals/pid_gmam_canis_lupus_familiaris/repeatmasker_all/ref_assembly/CanFam3_1.fasta.ori.out"
    output_file = "/home/akomissarov/Dropbox/PySatDNA/canis.result.tsv"   

    file_name = "/storageSM2/akomissarov/genomes_mammals/pid_gmam_panthera_tigris_altaica/repeatmasker_all/ref_assembly/PanTig1_0.fasta.ori.out"
    output_file = "/home/akomissarov/Dropbox/PySatDNA/tiger.result.tsv"            


    file_name = "/stripe/akomissarov/cheetah/repeatmasker_all/ref_assembly/cheetah.v3.scaffold.fa.ori.out"
    genome_file = "/stripe/akomissarov/cheetah/ref_data/v3/cheetah.v3.scaffold.fa"
    output_file = "/home/akomissarov/Dropbox/PySatDNA/cheetah.result.tsv"
    output_folder = "/stripe/akomissarov/cheetah/repeatmasker_all/ref_assembly/families"

    # file_name = "/storageSM2/akomissarov/genomes_mammals/pid_gmam_ailuropoda_melanoleuca/repeatmasker_all/ref_assembly/AilMel_1_0.fasta.ori.out"
    # output_file = "/home/akomissarov/Dropbox/PySatDNA/panda.result.tsv" 
    # genome_file = "/storageSM2/akomissarov/genomes_mammals/pid_gmam_ailuropoda_melanoleuca/ref_data/AilMel_1_0/AilMel_1_0.fasta"
    # output_folder = "/storageSM2/akomissarov/genomes_mammals/pid_gmam_ailuropoda_melanoleuca/repeatmasker_all/ref_assembly/families"

    # file_name = "/storage1/akomissarov/pacbio/repeatmasker/pacbio2_rna/L_RNA_scaffolder.fasta.ori.out"
    # output_file = "/home/akomissarov/Dropbox/PySatDNA/seabass.result.tsv"    
    # genome_file = "/storage1/akomissarov/pacbio/ref_data/pacbio2_trans/L_RNA_scaffolder.fasta"
    # output_folder = "/storage1/akomissarov/pacbio/repeatmasker/pacbio2_rna/families"


    # file_name = "/storage1/akomissarov/pacbio/repeatscout/pacbio2_rna/repeat_masker/L_RNA_scaffolder.fasta.ori.out"
    # genome_file = "/storage1/akomissarov/pacbio/ref_data/pacbio2_trans/L_RNA_scaffolder.fasta"
    # output_file = "/home/akomissarov/Dropbox/PySatDNA/seabass.denovo.result.tsv"
    # output_folder = "/storage1/akomissarov/pacbio/repeatscout/pacbio2_rna/repeat_masker/families"

    # fix_repeatmasker_out_file(file_name)

    # genome = {}
    # for i, seq_obj in enumerate(sc_iter_fasta(genome_file)):
    #     head = seq_obj.header.split()[0]
    #     print "Read:", i, head, "\r",
    #     genome[head] = seq_obj.sequence
    # print

    # compute_repbase_stats(file_name, output_file, genome, output_folder)


    from trseeker.tools.ngrams_tools import process_list_to_kmer_index
    from trseeker.tools.kmer_distance import compute_distances_for_index, compute_distances_for_index_by_raw_kmers, _process_index_data_to_file


    fasta_file = "/stripe/akomissarov/cheetah/repeatmasker_all/ref_assembly/SINEs.all.fa"
    dist_file = "/stripe/akomissarov/cheetah/repeatmasker_all/ref_assembly/SINEs.all.dist"
    ngram_index_file = "/stripe/akomissarov/cheetah/repeatmasker_all/ref_assembly/SINEs.all.index"
    k = 23



    sequences = []
    id2length = {}
    id2seq = {}
    id2head = {}
    print "Read lengths"
    for i, seq_obj in enumerate(sc_iter_fasta(fasta_file)):
        print "Read fasta", i, "\r",
        id2length[i] = seq_obj.length
        sequences.append(seq_obj.sequence)
        id2seq[i] = seq_obj.sequence
        id2head[i] = seq_obj.header
        # if i > 100:
        #     break
    print
    print "Loaded", len(id2length)

    from trseeker.models.ngrams_model import read_kmer_index
    from trseeker.tools.kmer_distance import compute_distances

    
    print "Read index.."
    index_objs = read_kmer_index(ngram_index_file, None)

    # with open(families_profile_file, "w") as fh:
    #     for i, consensus in index_data.items():
    #         fh.write("%s\t%s\n" % (i, consensus))
    print "Compute distances..."
    D = compute_distances(index_objs, id2length, index_function=None)


    # print "Compute index"
    # index_data = process_list_to_kmer_index(sequences, k, docids=True, cutoff=None, verbose=True)
    # print "Compute distances"
    # print "Sort data..."
    # result = []
    # skipped_by_dust = 0
    # for i, (key, revkey, tf, df, docids, freqs) in enumerate(index_data):
    #     new_doc_ids = docids
    #     # sort docsis by freqs
    #     items = []
    #     for pos in xrange(len(new_doc_ids)):
    #         all_items = id2length[pos] - k + 1
    #         if all_items <= 0:
    #             continue
    #         items.append((
    #                 freqs[pos] * 1. / all_items,
    #                 new_doc_ids[pos]
    #             ))
    #     items.sort()
    #     new_doc_ids = [str(x[1]) for x in items]
    #     freqs = [str(x[0]) for x in items]

    #     data = [key, revkey, tf, df, ",".join(new_doc_ids), ",".join(freqs)]
    #     data = "%s\n" % "\t".join(map(str, data))
    #     result.append(data)
    # print "Save data to %s..." % ngram_index_file
    # with open(ngram_index_file, "w") as fh:
    #     fh.writelines(result)



    # path_to_distance_exe = "/home/akomissarov/libs/compute_distance"
    # ## CALL C++ fuction
    # command = "%s %s %s" % (path_to_distance_exe, ngram_index_file, dist_file)
    # print command
    # os.system(command)



    # D = compute_distances_for_index(index, id2length, index_function=None)

    # data = [key, revkey, tf, df, ",".join(new_doc_ids), ",".join(freqs)]

    # r = []
    # print "Sort distances"
    # for key,d in D.items():
    #     a,b = map(int,key.split("\t"))
    #     r.append((d,a,b))
    # r.sort(reverse=True)
    # print "Save distances"
    # with open(dist_file, "w") as fh:
    #     for d,a,b in r:
    #         print d,a,b, id2head[a], id2head[b]
    #         s = "%s\t%s\t%s\n" % (id2head[a], id2head[b], d)
    #         fh.write(s)


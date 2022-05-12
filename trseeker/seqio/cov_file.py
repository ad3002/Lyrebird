#!/usr/bin/env python# -*- coding: utf-8 -*-
#
#@created: 05.06.2011
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com 
"""
Cov file assembly stack:

1. Check if is there any assmebly flank:

    NNN AAA
        or
    AAA
        AAA is last one

2. Take next and check if it is possible to extend:

    a) not possible
        AAA NNN BBB
            BBB is last one
    b) possible only left
        AAA aaN BBB
            BBB is last one
    c) possible only right
        AAA Nbb BBB
            BBB is last one
    d) possible both
        AAA ccc BBB
            BBB is last one
3. Fill the last one
    
    BBB NNN
        or
    BBB
        THE END.
    
"""
import os



class CovData(object):

    def __init__(self, data):
        self.type = data[0]
        self.pos = int(data[1])
        self.cov = int(data[2])
        self.kmer = data[3]
        self.pA = int(data[5])
        self.pC = int(data[6])
        self.pG = int(data[7])
        self.pT = int(data[8])
        self.nA = int(data[10])
        self.nC = int(data[11])
        self.nG = int(data[12])
        self.nT = int(data[13])
        self.pS = sum([self.pA, self.pC, self.pG, self.pT])
        self.nS = sum([self.nA, self.nC, self.nG, self.nT])

    def __str__(self):
        return "%s: %s-%s" % (self.pos, self.kmer, self.cov)


def sc_iter_cov_file(file_name):
    """ Iter over cov file."""
    with open(file_name) as fh:
        for line in fh:
            d = line.strip().split("\t")
            if not d:
                continue
            if d[2] == "-":
                continue    
            yield CovData(d)


def get_seq_dashes(file_name, verbose=False):

    seq = None
    last_pos = None
    result = []
    start = None
    end = None
    start_ext_cov = None
    end_ext_cov = None
    for x in sc_iter_cov_file(file_name):
        if seq is None:
            start = x.pos
            start_ext_cov = x.pS
            seq = x.kmer
            last_pos = x.pos
            end_ext_cov = x.nS
            continue
        if x.pos-1 == last_pos:
            last_pos = x.pos
            seq += x.kmer[-1]
            end_ext_cov = x.nS
        else:
            end = last_pos
            result.append((start, end, start_ext_cov, end_ext_cov, seq))
            if verbose:
                print len(result), start, end, len(seq)
            last_pos = None
            last_pos = x.pos
            start_ext_cov = x.pS
            start = x.pos
            seq = x.kmer
            end_ext_cov = x.nS
    if last_pos is not None:
        end = last_pos
        result.append((start, end, start_ext_cov, end_ext_cov, seq))
        if verbose:
            print len(result), start, end, len(seq)
    return result


if __name__ == '__main__':
    
    import jellyfish
    from trseeker.tools.jellyfish_tools import Kmer2tfAPI
    from AriadnaPy.tools.patches import get_patches_from_debrijin_simple
    from AriadnaPy.tools.unitigs import extend_to_right, extend_to_right_from_seed, extend_to_left_from_seed

    total_file = "/home/akomissarov/Dropbox/PySatDNA/all.seq"
    with open(total_file, "w") as fh:
        pass

    for sid in xrange(1,47):

        file_name = "/mnt/GNDRbackup/akomissarov/antilopa/ref.%s.cov" % sid
        jf_path = "/mnt/GNDRbackup/akomissarov/antilopa/trim.%s" % sid
        final_file = "/mnt/GNDRbackup/akomissarov/antilopa/fin.%s.seq" % sid

        result = get_seq_dashes(file_name, verbose=True)
        k = 23

        jf_api = jellyfish.QueryMerFile(jf_path)
        kmer2tf = Kmer2tfAPI(jf_api)

        sequences = ""
        pE1 = None
        end1 = None

        for x in result:
            print x[:4]

        

        for i in xrange(len(result)):

            start2, end2, pS2, pE2, seq2 = result[i]

            print sequences, end1, pE1
            print start2, end2, pS2, pE2, seq2
            # raw_input("?")

            ### Case 1. The first round.
            if not sequences:
                if start2 > 1:
                    if pS2 == 0:
                        sequences += "N"*start2
                        sequences += seq2
                        # raw_input("1a?")
                    else:
                        kmer = seq2[:k]
                        seq, cov, R = extend_to_left_from_seed(kmer, kmer2tf, error_cutoff=0, return_r=True)
                        sequences = seq[-start2:] + seq2
                        # raw_input("1b?")
                    pE1 = pE2
                    end1 = end2
                    
                    continue

            ### Case 2a. The first round.
            if pE1 == 0 and pS2 == 0:
                sequences += "N"*(start2-end1)
                sequences += seq2
                pE1 = pE2
                end1 = end2
                # raw_input("2a?")
                continue

            ### Case 2b. The first round.
            if pE1 > 0 and pS2 == 0:
                
                kmer = sequences[-k:]
                print sequences
                seq, cov, R = extend_to_right_from_seed(kmer, kmer2tf, error_cutoff=0, return_r=True)
                patch = seq[k:]
                sequences += patch
                sequences += "N"*(start2-end1-len(patch))
                sequences += seq2
                pE1 = pE2
                end1 = end2
                # raw_input("ext to right?")
                # raw_input("2b?")
                continue

            ### Case 2c. The first round.
            if pE1 == 0 and pS2 > 0:
                kmer = seq2[:k]
                print seq2
                seq, cov, R = extend_to_left_from_seed(kmer, kmer2tf, error_cutoff=0, return_r=True)
                patch = seq[:-k]
                
                sequences += "N"*(start2-end1-len(patch))
                sequences += patch
                sequences += seq2
                pE1 = pE2
                end1 = end2
                print seq
                # raw_input("ext to left?")
                # raw_input("2c?")
                continue

            ### Case 2d. The first round.
            if pE1 > 0 and pS2 > 0:

                left = sequences[-k:]
                right = seq2[:k]
                error_rate = 0
                max_length = 10000
                max_steps = 10000
                max_queue = 10000
                ext, seq2status, status = get_patches_from_debrijin_simple(kmer2tf, left, right, error_rate, max_length, max_steps, max_queue, k=23)

                if ext:
                    print pE1, pS2, ext, seq2status, status
                    patch = ext[0][23:]
                    sequences += patch
                    sequences += seq2
                    pE1 = pE2
                    end1 = end2
                    # raw_input("solid?")
                else:
                    ### A. ext left
                    kmer = sequences[-k:]
                    seq, cov, R = extend_to_right_from_seed(kmer, kmer2tf, error_cutoff=0, return_r=True)
                    left_patch = seq[k:]
                    
                    kmer = seq2[:k]
                    seq, cov, R = extend_to_left_from_seed(kmer, kmer2tf, error_cutoff=0, return_r=True)
                    right_patch = seq[:-k]
                    
                    n_patch = start2 - end1 - len(left_patch) - len(right_patch)

                    if n_patch <= 0:

                        print n_patch,start2 ,end1, start2-end1
                        print left_patch
                        print right_patch

                        sequences += "N"*(start2 - end1)
                        sequences += seq2

                        print "Error with insert length"
                        # raw_input("???")

                    else:

                        assert n_patch > 0

                        sequences += left_patch
                        sequences += "N"*n_patch
                        sequences += right_patch
                        sequences += seq2

                    pE1 = pE2
                    end1 = end2

                    # raw_input("not solid?")
                # raw_input("2d?")
                continue

        ### Case 3. The first round.
        # if pE1 > 0 and pS2 == 0:
        #     sequences += "N"*(start2-end1)
        #     sequences += seq2
        #     pE1 = pE2
        #     end1 = end2
            # raw_input("3?")

        with open(final_file, "w") as fh:
            fh.write(sequences)

        with open(total_file, "a") as fh:
            fh.write(">%s\n%s\n" % (sid, sequences))
            
        
        
        
            




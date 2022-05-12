#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 23.02.2011
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com 
'''
Functions related to ngram (k-mer).
'''
from collections import defaultdict
from trseeker.tools.sequence_tools import get_revcomp
from trseeker.seqio.tab_file import sc_iter_tab_file
from trseeker.models.trf_model import TRModel
from trseeker.seqio.fasta_file import sc_iter_fasta
from trseeker.tools.complexity import get_zlib_complexity


def generate_ngrams(text, n=12):
    ''' Yields all ngrams of length k from given text. 
    
    - n: ngram length
    '''
    for i in xrange(0, len(text) - n + 1):
        yield i, text[i:i + n]

def generate_kmers(seq, k):
    ''' Yields all kmers of length k from given seq. 
    '''
    for i in xrange(len(seq) - k + 1):
        yield seq[i:i + k]


def generate_window(text, n=12, step=None):
    ''' Yields all ngrams of length n from given text. 
    - n: ngram length
    - step: window step
    '''
    if not step:
        step = n / 2
    for i in xrange(0, len(text) - n + 1, step):
        yield i, text[i:i + n]

def get_ngrams(text, m=None, n=23, k=None, skip_n=False):
    ''' Returns m most frequent (ngram of length n, tf) tuples for given text.
    
    - m: number of returned ngrams
    - n: ngram length
    '''
    if k:
        n = k
    ngrams = {}
    for pos, ngram in generate_ngrams(text, n=n):
        if skip_n and 'n' in ngram:
            continue
        ngrams.setdefault(ngram, 0)
        ngrams[ngram] += 1
    ngrams = [(key, value) for key, value in ngrams.items()]
    ngrams.sort(reverse=True, key=lambda x: x[1])
    if m:
        return ngrams[:m]
    else:
        return ngrams

def get_kmer2tf(sequence, k):
    ''' Return kmer to tf dictionary.
    '''
    kmers = defaultdict(int)
    for i in xrange(0, len(sequence) - k + 1):
        kmer = sequence[i:i+k].lower()
        if 'n' in kmer:
            continue
        kmers[kmer] += 1
    return kmers


def get_kmers(sequence, k):
    ''' Return kmers.
    '''
    kmers = set()
    for i in xrange(0, len(sequence) - k + 1):
        kmer = sequence[i:i+k].lower()
        if 'n' in kmer:
            continue
        kmers.add(kmer)
    return list(kmers)

def get_kmers_both_strands(sequence, k):
    ''' Return kmers.
    '''
    kmers = set()
    for i in xrange(0, len(sequence) - k + 1):
        kmer = sequence[i:i+k].lower()
        if 'n' in kmer:
            continue
        kmers.add(kmer)
        kmers.add(get_revcomp(kmer))
    return list(kmers)

def get_ngrams_freq(text, m=None, n=23, k=None):
    ''' Returns m most frequent (ngram of length n, fraction of possible ngrams) tuples for given text.
    
    - m: number of returned ngrams
    - n: ngram length
    '''
    if k:
        n = k
    ngrams = defaultdict(int)
    for pos, ngram in generate_ngrams(text, n=n):
        ngrams[ngram] += 1
    text_length = float(len(text) - n + 1)
    ngrams = [(key, value, value / text_length) for key, value in ngrams.items()]
    ngrams.sort(reverse=True, key=lambda x: x[1])
    if m:
        return ngrams[:m]
    return ngrams

def get_most_freq_kmer(text, k=None):
    ''' Returns the most frequent kmer
    '''
    ngrams = defaultdict(int)
    for pos, ngram in generate_ngrams(text, n=k):
        ngrams[ngram] += 1
    ngrams = [(key, value, k*value) for key, value in ngrams.items()]
    ngrams.sort(reverse=True, key=lambda x: x[1])
    return ngrams[0]

def get_ngrams_feature_set(text, m=5, n=12):
    ''' Returns a feature set {'ngram':'ngram',...}  of m most frequent ngram of length n for given text.
    
    - m: number of returned ngrams
    - n: ngram length
    '''

    ngrams = {}
    for pos, ngram in generate_ngrams(text, n=n):
        ngrams.setdefault(ngram, 0)
        ngrams[ngram] += 1
    result = [(key, value) for key, value in ngrams.items()]
    ngrams.sort(reverse=True, key=lambda x: x[1])
    data = {}
    for key, value in enumerate(result[:m]):
        data[key] = key
    return data

def get_ngram_freq_distance(ngrams_a, ngrams_b):
    ''' Returns a distance between two ngram sets where distance is a sum(min(ngram_a, ngram_b) for each common ngram)
    
    - ngrams_a: dictionary {ngram:n, ...}
    - ngrams_b: dictionary {ngram:n, ...}
    '''
    distance = 0
    common_ngrams = [ngram for ngram, fr in ngrams_a.items() if ngram in ngrams_b]
    for ngram in common_ngrams:
        distance += min(ngrams_a[ngram], ngrams_b[ngram])
    return distance

def get_ngram_common_distance(ngrams_a, ngrams_b):
    ''' Returns a distance between two ngram sets where distance is a len()
    
    - ngrams_a: dictionary {ngram:n, ...}
    - ngrams_b: dictionary {ngram:n, ...}
    '''
    common_ngrams = [ngram for ngram, fr in ngrams_a.items() if ngram in ngrams_b]
    return len(common_ngrams)

def get_repeatness_coefficent(length, k, kmern):
    ''' Return repeatness coefficient. From 0 (e.g. polyA) to 1 (unique sequence).
    '''
    N = length - k + 1.
    return kmern * (N-1) / (N**2)

def get_expessiveness_coefficent(kmern, k):
    ''' Return expressivenesss coefficient.
    '''
    return kmern * 1. / 4^k

def count_kmer_tfdf(sequence, tf_dict, df_dict, k):
    ''' Update tf and df data with k-mers from given sequence.
    '''
    seen = set()
    local_tf = defaultdict(int)
    local_df = defaultdict(int)
    sequence = sequence.lower()
    for (ngram, tf, nf) in get_ngrams_freq(sequence, n=k):
        if 'n' in ngram:
            continue
        seen.add(ngram)
        tf_dict[ngram] += tf
        local_tf[ngram] += tf
    for ngram in seen:
        df_dict[ngram] += 1
        local_df[ngram] += 1
    return tf_dict, df_dict, local_tf, local_df

def get_kmer_tf_df_for_data(data, k, docids=False, verbose=True):
    '''
    '''
    df = defaultdict(int)
    tf = defaultdict(int)
    kmer2ids = defaultdict(list)
    kmer2freq = defaultdict(list)
    N = len(data)
    if N>100:
        verbose = True
    for i, sequence in enumerate(data):
        if verbose:
            print("Process tf/df: ", i, N, sep=" ")
        tf, df, local_tf, local_df = count_kmer_tfdf(sequence, tf, df, k)
        if docids:
            for key in local_tf:
                kmer2ids[key].append(i)
                kmer2freq[key].append(local_tf[key])
    if verbose:
        print()
    if docids:
        return tf, df, kmer2ids, kmer2freq
    return tf, df

def get_kmer_tf_df_pos_for_list(data, k):
    '''
    '''
    df = defaultdict(int)
    tf = defaultdict(int)
    kmer2pos = defaultdict(list)
    for docid, sequence in enumerate(data):
        df_set = set()
        for i in xrange(len(sequence)-k+1):
            kmer = sequence[i:i+k]
            kmer2pos[kmer].append((docid,i))
            tf[kmer] += 1
            if not kmer in df_set:
                df_set.add(kmer)
                df[kmer] += 1
    return (tf, df, kmer2pos)

    
def get_df_stats_for_list(data, k, kmer2df):
    ''' Compute max df, number and percentage of sequence with given ngram.
    Return (maxdf, nmaxdf, pmaxdf)
    '''
    df = defaultdict(int)
    tf = defaultdict(int)
    n = len(data)
    ngram_seqs = []
    for sequence in data:
        tf, df, local_tf, local_df = count_kmer_tfdf(sequence, tf, df, k)
    result = [(v,k) for (k,v) in df.items()]
    result.sort()
    maxdf = result[-1][0]
    ngram_seqs = [(k, tf[k]) for v,k in result if v == maxdf]
    ngram_seqs.sort(key=lambda x: x[1], reverse=True)
    nmaxdf = len(ngram_seqs)
    pmaxdf = round(float(maxdf)/n, 3)
    ngram_seqs = [":".join((k,str(f))) for k,f in ngram_seqs[:10]]
    return (maxdf, nmaxdf, pmaxdf, ngram_seqs)

def process_list_to_kmer_index(data, k, docids=True, cutoff=None, verbose=True):
    ''' Get list of string.
    Return list of (kmer, revkmer, tf, df, None, None)
    OR
    Return list of (kmer, revkmer, tf, df, docids, freqs)
    '''
    if docids:
        (tf_dict, df_dict, doc_data, freq_data) = get_kmer_tf_df_for_data(data, k, docids=docids, verbose=verbose)
    else:
        (tf_dict, df_dict) = get_kmer_tf_df_for_data(data, k, docids=docids, verbose=verbose)
    result = []
    seen = set()
    skipped = 0
    added = 0
    if verbose:
        print("Join data...")
    for key in df_dict:
        if verbose:
            print(skipped, added, sep=" ")
        if key in seen:
            continue
        revkey = get_revcomp(key)
        if revkey in seen:
            continue
        if revkey in df_dict:
            df = df_dict[key] + df_dict[revkey]
            tf = tf_dict[key] + tf_dict[revkey]
            if docids:
                ids = doc_data[key] + doc_data[revkey]
                freqs = freq_data[key] + freq_data[revkey]
        else:
            df = df_dict[key]
            tf = tf_dict[key]
            if docids:
                ids = doc_data[key]
                freqs = freq_data[key]
        # skip by df
        if cutoff and df <= cutoff:
            skipped += 1
            continue
        added += 1
        if revkey < key:
            key, revkey = revkey, key
        if docids:
            result.append((key, revkey, tf, df, ids, freqs))
        else:
            result.append((key, revkey, tf, df, None, None))
        seen.add(key)
        seen.add(revkey)
    if verbose:        
        print()
    result.sort(key=lambda x: x[-3], reverse=True)
    return result

def compute_kmer_index_for_fasta_file(file_name, index_file, k=23):
    """ 
    """
    data = []
    print("Read arrays...")
    for i, seq_obj in enumerate(sc_iter_fasta(file_name)):
        data.append(seq_obj.sequence)
    print("Readed %s arrays." % i)
    print("Compute k-mers...")
    result = process_list_to_kmer_index(data, k, docids=False)
    print("Save index...")
    with open(index_file, "w") as fh:
        for item in result:
            s  = "%s\n" % "\t".join(map(str, item))
            fh.write(s)
    return result

def compute_kmer_index_for_trf_file(file_name, index_file, k=23, max_complexity=None, min_complexity=None):
    """
    TODO: replace gzip complexity with DUST filter
    """
    data = []
    print("Read arrays...")
    for i, trf_obj in enumerate(sc_iter_tab_file(file_name, TRModel)):
        data.append(trf_obj.trf_array)
    print("Readed %s arrays." % i)
    print("Compute k-mers...")
    result = process_list_to_kmer_index(data, k, docids=False)
    print("Save index...")
    with open(index_file, "w") as fh:
        for item in result:
            if max_complexity:
                if get_zlib_complexity(item[0]) > max_complexity:
                    continue
            if min_complexity:
                if get_zlib_complexity(item[0]) < min_complexity:
                    continue
            s  = "%s\n" % "\t".join(map(str, item))
            fh.write(s)
    return result

def get_sequence_kmer_coverage(sequence, kmers, k):
    ''' 
    '''
    n = len(sequence)
    match = 0.
    mismatch = 0.
    variability = set()
    for i, kmer in generate_ngrams(sequence, n=k):
        if kmer in kmers:
            variability.add(kmer)
            match += 1
        else:
            mismatch += 1
    return match/(n-k+1), variability

def compute_kmers_libraries_from_fasta(fasta_file, k_diaposon):
    '''
    '''
    index2name = {}
    arrays = []
    for i, seq_obj in enumerate(sc_iter_fasta(fasta_file)):
        index2name[i] = seq_obj.seq_head[1:]
        arrays.append(seq_obj.sequence)
    libraries = {}
    for k in k_diaposon:
        library = {}
        index = process_list_to_kmer_index(arrays, k, docids=True)
        for kmer, revkmer, tf, df, docids, freqs in index:
            items = []
            for i, docid in enumerate(docids):
                name = index2name[docid]
                freq = freqs[i]
                items.append("%s:%s" % (name, freq))
            items = ",".join(items)
            library[kmer] = (tf, df, items)
        libraries[k] = library
    return libraries

def get_for_and_rev_kmers_from_fastq(fasta_file, k):
    '''
    '''
    index2name = {}
    arrays = []
    for i, seq_obj in enumerate(sc_iter_fasta(fasta_file)):
        index2name[i] = seq_obj.seq_head[1:]
        arrays.append(seq_obj.sequence)
    library = {}
    index = process_list_to_kmer_index(arrays, k, docids=True)
    for kmer, revkmer, tf, df, docids, freqs in index:
        items = []
        for i, docid in enumerate(docids):
            name = index2name[docid]
            freq = freqs[i]
            items.append("%s:%s" % (name, freq))
        items = ",".join(items)
        library[kmer] = (tf, df, items)
        library[revkmer] = (tf, df, items)
    return library

def print_prev_cutoff(kmer, kmer2freq, cutoff=0):
    ''' Returns:
    - R as a dictionary nucleotide to tf
    - n as a number of pathes
    - nucleotides as a list with nucleotides
    - max_hits as a sorted list of (tf, nucleotide)
    '''
    R = {}
    n = 0
    nucleotides = []
    max_hits = []
    a = kmer2freq['A'+kmer[:-1]]
    if a and a > cutoff:
        R["A"] = a
        n += 1
        nucleotides.append('A')
        max_hits.append((a,'A'))
    c = kmer2freq['C'+kmer[:-1]]
    if c and c > cutoff:
        R["C"] = c
        n += 1
        nucleotides.append('C')
        max_hits.append((c,'C'))
    t = kmer2freq['T'+kmer[:-1]]
    if t and t > cutoff:
        R["T"] = t
        n += 1
        nucleotides.append('T')
        max_hits.append((t,'T'))
    g = kmer2freq['G'+kmer[:-1]]
    if g and g > cutoff:
        R["G"] = g
        n += 1
        nucleotides.append('G')
        max_hits.append((g,'G'))

    max_hits.sort(reverse=True)
    return R, n, nucleotides, max_hits

def print_next_cutoff(kmer, kmer2freq, cutoff=0):
    ''' Returns
    - R as a dictionary nucleotide to tf
    - n as a number of pathes
    - nucleotides as a list with nucleotides
    - max_hits as a sorted list of (tf, nucleotide)
    '''
    R = {}
    n = 0
    nucleotides = []
    max_hits = []
    a = kmer2freq[kmer[1:]+'A']
    if a and a > cutoff:
        R["A"] = a
        nucleotides.append('A')
        n += 1
        max_hits.append((a,'A'))
    c = kmer2freq[kmer[1:]+'C']
    if c and c > cutoff:
        R["C"] = c
        nucleotides.append('C')
        n += 1
        max_hits.append((c,'C'))
    t = kmer2freq[kmer[1:]+'T']
    if t and t > cutoff:
        R["T"] = t
        nucleotides.append('T')
        max_hits.append((t,'T'))
        n += 1
    g = kmer2freq[kmer[1:]+'G']
    if g and g > cutoff:
        R["G"] = g
        nucleotides.append('G')
        n += 1
        max_hits.append((g,'G'))
        
    max_hits.sort(reverse=True)
    return R, n, nucleotides, max_hits


def print_next_R(kmer, kmer2freq, cutoff=0):
    ''' Returns
    - R as a dictionary nucleotide to tf
    '''
    R = {}
    a = kmer2freq[kmer[1:]+'A']
    if a and a > cutoff:
        R["A"] = a
    c = kmer2freq[kmer[1:]+'C']
    if c and c > cutoff:
        R["C"] = c
    t = kmer2freq[kmer[1:]+'T']
    if t and t > cutoff:
        R["T"] = t
    g = kmer2freq[kmer[1:]+'G']
    if g and g > cutoff:
        R["G"] = g
    return R


def print_prev_R(kmer, kmer2freq, cutoff=0):
    ''' Returns:
    - R as a dictionary nucleotide to tf
    '''
    R = {}
    a = kmer2freq['A'+kmer[:-1]]
    if a and a > cutoff:
        R["A"] = a
    c = kmer2freq['C'+kmer[:-1]]
    if c and c > cutoff:
        R["C"] = c
    t = kmer2freq['T'+kmer[:-1]]
    if t and t > cutoff:
        R["T"] = t
    g = kmer2freq['G'+kmer[:-1]]
    if g and g > cutoff:
        R["G"] = g
    return R
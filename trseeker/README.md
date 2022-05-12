# Trseeker: a Python library and framework for working with various genomics data

## Contents

- [Introduction](#_toolkit_intro)
- [Settings](#_toolkit_settings)
- [Models](#_models)
    - [DNA Sequence](#_models_seq)
    - [Tandem repeat](#_models_trf)
    - [Organism](#_models_organism)
    - [Dataset](#_models_dataset)
    - [Blast output](#_models_blast)
    - [Chromosome](#_models_chr)
    - [WGS assembly](#_models_wgs)
    - [Genome](#_models_genome)
    - [Kmer](#_models_kmer)
    - [Kmers](#_models_kmers)
- [Readers](#_io)
    - [Tab file](#_io_tab)
    - [Block file](#_io_block)
    - [Fasta file](#_io_fasta)
    - [GBFF file](#_io_gbff)
    - [FTP](#_io_ftp)
    - [NCBI FTP](#_io_ncbi_ftp)
    - [MONGO DB](#_io_mongo)
    - [TRF output](#_io_trf)
    - [TRs file](#_io_trs)
    - [Ngram file](#_io_ngram)
    - [Blast output](#_io_blast)
    - [Fastq file](#_io_sra)
- [Toolkit](#_tools)
    - [Annotation tools](#_annotation)
    - [Assembly tools](#_assembly_tools)
    - [Tulip](#_tools_tulip)
    - [TRs nomenclature](#_tools_trs_group)
    - [Sequence](#_tools_sequence)
    - [Various for files](#_tools_var)
    - [Various](#_tools_other)
    - [Kmers](#_tools_kmers)
    - [Edit distance](#_tools_ed)
    - [Repbase](#_tools_repbase)
    - [Trace](#_tools_trace)
    - [Statistics](#_tools_stat)
    - [Parsing](#_tools_parsers)
    - [TRF](#_tools_trf)
    - [TRs datasets](#_tools_trs)
    - [SRA datasets](#_tools_sra)
    - [Sequence patterns](#_tools_patterns)
    - [Blast](#_tools_blast)
    - [ESA](#_tools_sa)
    - [NCBI genomes](#_tools_ncbi_genome)
    - [Networks](#_tools_networks)
    - [Networks](#_tools_ncbi_ann)
    - [Jellyfish](#_tools_jellyfish)
    - [Trees](#_tools_tree)
    - [Classifiction TRs in types](#_trs_types)
    - [Working with Entrez NCBI databases](#_tools_entrez)



<a name="_toolkit_intro"/>

## Introduction

Framework состоит из следующих частей:

- модели данных
- чтение/запись данных согласно моделям
- инструменты работы с этими данными

<a name="_toolkit_settings"/>

## Настройки фреймворка

В файле settings.py:

```python
SETTINGS_FILENAME = "settings.yaml"
NGRAM_LENGTH = 23
NGRAM_N = 100000000
```

Настройки можно прочитать и записать:

```python
from trseeker.settings import load_settings
from trseeker.settings import save_settings

settings_dict = load_settings()
save_settings(settings_dict)
```

Настройки содержат следующие параметры:

```yaml	
trseeker:
    os: linux64
    root_dir: /root/Dropbox/workspace/trseeker
    work_dir: /home
blast_settings:
    blast_location:
    blast_location_NIX: /home/ncbi-blast-2.2.26+/bin/legacy_blast.pl bl2seq
    blast_location_WIN: c:\Program Files\NCBI\blast-2.2.26+\bin\legacy_blast.pl bl2seq.exe
    jellyfish_location: jellyfish
    repbase_db_folder: /home/rebase_blast_db/repbase
    blast_e: 1e-20
    blast_b: 20000000
    blast_v: 20000000
trf_settings:
    trf_location: /root/trf404.linux64
    trf_match: 2
    trf_mismatch: 5
    trf_indel: 7
    trf_p: 80
    trf_q: 10
    trf_threshold: 50
    trf_length: 2000
    trf_param_postfix: 2.5.7.80.10.50.2000
    trf_masked_file: False
    trf_flanked_data: True
    trf_data_file: True
    trf_nohtml: True
    overlapping_cutoff: 10
    overlapping_cutoff_proc: 30
    overlapping_gc_diff: 0.05
ngrams_settings:
    ngram_length: 23
    ngram_m: 10000000
```
<a name="_models"/>

## Avaliable Models


<a name="_models_seq"/>

### DNA Sequence

```python
from trseeker.models.sequence_model import SequenceModel
```

Attribites:

- seq_gi (int)
- seq_ref
- seq_description
- seq_sequence
- seq_length (int)
- seq_gc (float)
- seq_revcom, reverce complement
- seq_gapped (int)
- seq_chr
- seq_head
- seq_start_position (int)
- seq_end_position (int)

Properties

- length (self.seq_length)
- sequence (self.seq_sequence)
- fasta
- header
- contige_coverage if **_cov_** key in name

```python
print seq_obj.fasta
>>> ">[seq_ref]\n[seq_sequence]\n"
```

- sa_input

```python
print seq_obj.sa_input
>>> "[seq_sequence]$"
```

- ncbi_fasta

```python
print seq_obj.ncbi_fasta
>>> ">gi|[seq_gi]|ref|[seq_ref]|[seq_description]\n[seq_sequence]\n"
```

Methods:

- add_sequence_revcom()
- set_dna_sequence(self, title, sequence, description=None)

```python
self.seq_ref = title
```

- set_ncbi_sequence(self, head, sequence)

Chromosome name is **?** or setted with parse_chromosome_name(head).

```python	
(self.seq_gi, self.seq_ref, self.seq_description) = parse_fasta_head(head)
```

- set_gbff_sequence(self, head, sequence)

Head is a dictionary with gi, ref, description keys.

Chromosome name is **?** or setted with parse_chromosome_name(head["description"]).

Sequence is cleared with clear_sequence(s) function. Lowercase and all non-DNA characters replacing with **n**. If the sequence has **n** then it is gapped.

<a name="_models_trf"/>

### TRF results

```python
from trseeker.models.trf_model import TRModel
```

Attributes:

- project, project name
- id (float)
- trf_id (int)
- trf_type,
- trf_family,
- trf_family_prob (float),
- trf_l_ind (int)
- trf_r_ind (int)
- trf_period (int)
- trf_n_copy (float)
- trf_pmatch (float)
- trf_pvar (float)
- trf_entropy (float)
- trf_consensus
- trf_array
- trf_array_gc (float)
- trf_consensus_gc (float)
- trf_gi
- trf_head
- trf_param
- trf_array_length (int)
- trf_chr
- trf_joined (int)
- trf_superfamily
- trf_superfamily_self
- trF_superfamily_ref
- trf_family
- trf_subfamily
- trf_subsubfamily
- trf_family_network
- trf_family_self
- trf_family_ref
- trf_hor (int)
- trf_n_chrun (int)
- trf_chr_refgenome
- trf_bands_refgenome
- trf_repbase
- trf_strand

Methods

- set_project_data(project), set self.project to given project
- set_raw_trf(head, body, line), head, body and line from TRF parser
- get_index_repr()

```python	
print trf_obj.get_index_repr()
'''
Tab delimted string with \n-symbol:
trf_id
trf_period
trf_array_length
trf_pvar
trf_gi
trf_l_ind
trf_r_ind
trf_chr
'''
```

- get_numerical_repr()

```python
print trf_obj.get_numerical_repr()
>>> [trf_period]\t[trf_array_length]\t[trf_array_gc]\n
```

- get_fasta_repr(), where head is trf_obj.trf_id and sequence is trf_obj.trf_array
- or **fasta** property
- get_monomer_fasta_repr(), where head is trf_obj.trf_id and sequence is trf_obj.trf_consensus
- get_family_repr()
- get_gff3_string(self, chromosome=True, trs_type="complex_tandem_repeat", probability=1000, tool="PySatDNA", prefix=None, properties={"id":"trf_id", "family": "trf_family_self"})

```python
print trf_obj.get_family_repr()
'''
Tab delimted string with \n-symbol:
trf_id
trf_period
trf_array_length
trf_array_gc
trf_pvar
trf_gi
trf_l_ind
trf_r_ind
trf_chr
trf_repbase
trf_superfamily
trf_family
trf_subfamily
'''
```

For network slice added one more index - gid (group id)

```python
from trseeker.models.trf_model import NetworkSliceModel	

slice_obj = NetworkSliceModel()
```

#### TRsClassificationModel(AbstractModel):

Model for keeping classification data.

Attributes:

- project
- id (int)
- trf_id (int)
- trf_period (int)
- trf_array_length (int)
- trf_array_gc
- trf_type
- trf_family
- trf_subfamily
- trf_family_prob (float)
- trf_family_kmer
- trf_subfamily_kmer
- trf_family_self
- class_ssr
- class_tssr
- class_sl
- class_good
- class_micro
- class_100bp
- class_perfect
- class_x4
- class_entropy
- class_gc
- trf_consensus

```python

class_obj = TRsClassificationModel()
class_obj.set_with_trs(trf_obj)

print class_obj.network_head()

```

<a name="_models_organism"/>

### Organism model

```python
from trseeker.models.organism_model import OrganismModel		
```

Attributes:

- organism_taxon
- organism_common_name
- organism_acronym
- organism_description
- organism_wgs_projects
- organism_genome_assemblies

<a name="_models_dataset"/>

### Dataset model

```python
from trseeker.models.dataset_model import DatasetModel		
```

Attributes:

- dataset_taxon
- dataset_id
- dataset_sources
- dataset_description
- dataset_gc (float)
- dataset_length (int)
- dataset_trs_n (int)
- dataset_trs_length (int)
- dataset_trs_mean_gc (float)
- dataset_trs_fraq (float)

<a name="_models_blast"/>

### Blast Results Model

```python
from trseeker.models.blast_model import BlastResultModel		
```

Attributes:

- query_id (int)
- query_gi (int)
- query_ref
- subject_id
- subject_gi(int)
- subject_ref
- query_start (int)
- query_end (int)
- subject_start (int)
- subject_end (int)
- evalue  (float)
- bit_score (flaot)
- score (int)
- alignment_length (int)
- proc_identity (float)
- identical (int)
- mismatches (int)
- positives (int)
- gap_opens (int)
- gaps (int)
- proc_positives (float)
- frames
- query_frame (int)
- subject_frame (int)
- fraction_of_query (float)  

Additional functions:

- read_blast_file(blast_file, length), return subject_ref -> list of matches (BlastResultModel models).

```python
from trseeker.models.blast_model import read_blast_file

ref_to_blast_obj = read_blast_file(file_name)
```

<a name="_models_chr"/>

### Chromosome model

```python
from trseeker.models.chromosome_model import ChromosomeModel
```

Attributes:

- chr_genome
- chr_number
- chr_taxon
- chr_prefix
- chr_gpid
- chr_acronym
- chr_contigs
- chr_length
- chr_mean_gc
- chr_trs_all
- chr_trs_3000
- chr_trs_all_proc
- chr_trs_3000_proc
- chr_trs_all_length
- chr_trs_3000_length
- genome_gaps
- chr_sum_gc

<a name="_models_wgs"/>

### WGS assembly model

```python
from trseeker.models.wgs_model import WGSModel
```

Attributes:

- wgs_prefix
- wgs_taxon
- wgs_gpid
- wgs_acronym
- wgs_contigs (int)
- wgs_length (int)
- wgs_mean_gc (float)
- wgs_trs_all (int)
- wgs_trs_3000 (int)
- wgs_trs_1000 (int)
- wgs_trs_500 (int)
- wgs_trs_all_proc (float)
- wgs_trs_3000_proc (float)
- wgs_trs_1000_proc (float)
- wgs_trs_500_proc (float)
- wgs_trs_all_length (int)
- wgs_trs_3000_length (int)
- wgs_trs_1000_length (int)
- wgs_trs_500_length (int)
- wgs_sum_gc (float)

Methods:

Clear trf information (set to 0):

```python
wgs_obj.clear_trf()
```

<a name="_models_genome"/>

### Genome model

```python
from trseeker.models.genome_model import GenomeModel
```

- genome_taxon
- genome_prefix
- genome_gpid
- genome_acronym
- genome_chromosomes
- genome_contigs
- genome_length
- genome_mean_gc
- genome_trs_all
- genome_trs_3000
- genome_trs_all_proc
- genome_trs_3000_proc
- genome_trs_all_length
- genome_trs_3000_length
- genome_gaps
- genome_sum_gc

### RepeatMasker track data

```python
from trseeker.models.genome_model import RepeatMaskerAnnotation
```

Container for RepeatMasker track. For example from http://hgdownload.cse.ucsc.edu/goldenPath/felCat5/database/rmsk.txt.gz

```
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
```

<a name="_models_kmer"/>

### Ngram/kmer model

```python
from trseeker.models.ngram_model import NgramModel
```

```python
ngram_obj = NgramModel(seq_f, seq_r)
ngram_obj.add_tr(trf_obj, tf)

print ngram_obj
>>> '[fseq]\t[rseq]\t[tf]\t[df]\t[len taxons]\t[len fams]\n'

print ngram_obj.get_families()
>>> ???
```

Attributes

- seq_r
- seq_f
- tf (float)
- df (int)
- taxons (set)
- trs (set)
- families (dict)

<a name="_models_kmers"/>

### Kmer index model

```python
from trseeker.models.ngrams_model import KmerIndexModel
```

Attributes:

- kmer
- rkmer
- tf (int)
- df (int)
- docs (list of int)
- freqs (list of float)

### Kmer Slice model

```python
from trseeker.models.ngrams_model import KmerSliceModel
```

#### Attributes

- kmer
- local_tf (int)
- df (int)
- f_trs
- f_wgs
- f_mwgs
- f_trace
- f_sra
- f_ref1
- f_ref2
- f_caroli

### SliceTreeModel

```python
from trseeker.models.ngrams_model import SliceTreeModel
```

#### Attributes

- deep (int)
- size (int)
- blast_fams (list of str)
- maxdf (int)
- nmaxdf (int)
- pmaxdf (float)
- gc_var (float)
- units (list of int)
- trs (list of int)
- kmers (list of SliceTreeModel)

#### Additional function

Yield KmerIndexModel:

```python
sc_ngram_reader(file_name)
```

Yield (ngram id, [(seq id, tf), ...]):

```python
sc_ngram_trid_reader(file_name)
```

Read kmer index data as list:

Parameters:

- index file
- microsatellite kmer dictionary, kmer->index


```python
read_kmer_index(ngram_index_file, micro_kmers, cutoff=1)
```

<a name="_io"/>

## IO functions


<a name="_io_tab"/>

### Tab file

```python
from trseeker.seqio.tab_file import TabDelimitedFileIO

reader = TabDelimitedFileIO(skip_first=False, format_func=None, delimiter='\t', skip_startswith='#')

reader.read_from_file(file_name)
```

Avaliable all functions from parent class AbstractFileIO from PyExp package.

Useful functions:

- sc_iter_tab_file(input_file, data_type, remove_starts_with=None)
- sc_iter_simple_tab_file(input_file)
- sc_read_dictionary(dict_file, value_func=None)
- sc_read_simple_tab_file(input_file)
- sc_write_model_to_tab_file(output_file, objs)

```python	
from trseeker.seqio.tab_file import sc_iter_tab_file

for wgs_obj in sc_iter_tab_file(file_name, WGSModel):
	print wgs_obj.wgs_taxon
```

```python
from trseeker.seqio.tab_file import sc_iter_simple_tab_file

for item in sc_iter_simple_tab_file(file_name):
	print item[0]
```

```python
from trseeker.seqio.tab_file import sc_read_dictionary

for data in sc_read_dictionary(file_name):
	for key, value in data.items():
		print key, value
```

```python
from trseeker.seqio.tab_file import sc_read_simple_tab_file

for data in sc_read_simple_tab_file(file_name):
	for item in data:
		print "id = ", item[0]	
```

<a name="_io_block"/>

### Block file

```python
from trseeker.seqio.block_file import AbstractBlockFileIO

reader = AbstractBlockFileIO(token, **args)
for (head, body, start, next) in reader.read_online(file_name):
	pirnt head
```

Avaliable all functions from parent class AbstractFileIO from PyExp package.

<a name="_io_fasta"/>

### Fasta file

```python
from trseeker.seqio.fasta_file import FastaFileIO

reader = FastaFileIO()
```

#### Useful functions:

- sc_iter_fasta(file_name)
- fasta_reader(file_name), same as sc_iter_fasta
- sc_iter_fasta_simple(file_name)
- save_all_seq_with_exact_substring(fasta_file, substring, output_file)
- sort_fasta_file_by_length(file_name)

```python
from trseeker.seqio.fasta_file import sc_iter_fasta

for seq_obj in sc_iter_fasta(file_name):
	print seq_obj.sequence
```

```python
from trseeker.seqio.fasta_file import sc_iter_fasta_simple

for (gi, sequence) in sc_iter_fasta_simple(file_name):
	print gi
```

```python
from trseeker.seqio.fasta_file import save_all_seq_with_exact_substring

save_all_seq_with_exact_substring(fasta_file, substring, output_file)
```

Sort by length and save fasta file:

```python
sort_fasta_file_by_length(file_name)
```

<a name="_io_gbff"/>

### GBFF file

```python
from trseeker.seqio.gbff_file import GbffFileIO

reader = GbffFileIO()
```

#### Useful functions:

- sc_iter_gbff(file_name)
- sc_iter_gbff_simple(file_name)
- sc_parse_gbff_in_folder(gbbf_folder, fasta_folder, fasta_postfix, mask)

```python
from trseeker.seqio.gbff_file import sc_iter_gbff

for seq_obj in sc_iter_gbff(file_name):
	print seq_obj.length
```

```python
from trseeker.seqio.gbff_file import sc_iter_gbff_simple

for (gi, sequence) in sc_iter_gbff_simple(file_name):
	print gi
```

```python
from trseeker.seqio.gbff_file import sc_parse_gbff_in_folder

sc_parse_gbff_in_folder("/home/user/", "/home/user/fasta", "fa", "mouse")
```


<a name="_io_ftp"/>

### FTP IO

```python	
from trseeker.seqio.frt_io import AbstractFtpIO

reader = AbstractFtpIO(ftp_address=address)

reader.connect()

reader.cd(['/home', 'user', 'data'])

reader.ls()
>>> ['readme.txt', 'data.fa']

file_size = reader.size(file_name)

reader.get(file, output_file)

reader.unzip(file_name)
```

Download from NCBI ftp with aspera:

```python
download_with_aspera_from_ncbi(source, destination)
```

<a name="_io_ncbi_ftp"/>

### NCBI ftp

```python
from trseeker.seqio.ncbi_ftp import NCBIFtpIO

reader = NCBIFtpIO()

reader.download_wgs_fasta(wgs_list, file_suffix, output_folder, unzip=False)

reader.download_all_wgs_in_fasta(output_folder, aspera=False, get_size=False, unzip=False)

reader.download_all_wgs_in_gbff(output_folder)

reader.download_chromosomes_in_fasta(ftp_address, name, output_folder)

reader.download_trace_data(taxon, output_folder)

reader.download_sra_from_ddbj(ftp_address, output_folder)

reader.download_with_aspera(local_path, remove_server, remote_path)
```

<a name="_io_mongo"/>

### Mongo db reader

Not implemented yet.

```python	
from trseeker.seqio.mongodb_reader import MongoDBReader

reader = MongoDBReader()
db = reader.get_trdb_conn()
```

<a name="_io_trf"/>

### TRF file

#### Useful functions:

- sc_parse_raw_trf_folder(trf_raw_folder, output_trf_file, project=None)

```python
from trseeker.seqio.trf_file import TRFFileIO

from trseeker.seqio.trf_file import sc_parse_raw_trf_folder

for trf_obj in TRFFileIO(file_name, filter=True):
	print trf_obj.trf_id
```

```python
sc_parse_raw_trf_folder(trf_raw_folder, output_trf_file, project="mouse_genome")
```

При чтение данных TRF происходит их фильтрация по следующим параметрам:

1. Убираются все вложенные поля меньшей длины.
2. Если поля overlapping, то сейчас ничего не делается, а раньше если 

```python
overlap_proc_diff >= settings["trf_settings"]["overlapping_cutoff_proc"] 
gc_dif <= settings["trf_settings"]["overlapping_gc_diff"]
``` 
поля объяединяются в одно поле. Иначе считаем, что это отдельные поля.
3. Если поля совпадают, то выбирается то которое с большим trf_pmatch.

TODO: убрать пересечение, так как это видимо в том числе и ошибки ассмеблера и это будет мешать корректной классификации. Убрать в отдельную часть.

<a name="_io_trs"/>
### TRs file

```python
from trseeker.seqio.tr_file import *
```

Functions:

- save_trs_dataset(trs_dataset, output_file)
- save_trs_class_dataset(tr_class_data, output_file)
- get_trfid_obj_dict(trf_large_file)
- get_classification_dict(fam_kmer_index_file)

Save trf file as fasta file:

```python
# don't save trf_obj if trf_family is ALPHA
save_trs_as_fasta(trf_file, fasta_file, add_project=False, skip_alpha=False)
```

```python
from trseeker.seqio.tr_file import get_trfid_obj_dict

trid2trfobj = get_trfid_obj_dict(trf_large_file)
```

- get_all_trf_objs(trf_large_file)
- get_all_class_objs(trf_class_file)
- get_class_objs_dict(trf_class_file)

```python
from trseeker.seqio.tr_file import get_all_trf_objs
from trseeker.seqio.tr_file import get_all_class_objs

trf_obj_list = get_all_trf_objs(trf_large_file)
class_obj_list = get_all_class_objs(trf_class_file)
```

- get_trf_objs_dict(trf_large_file)

```python
trf_obj_dics = get_trf_objs_dict(trf_large_file)
```

- read_trid2meta(file_name)

Load trid to full index dictionary as string:

```python
read_trid2ngrams(annotation_ngram_folder, trf_large_file)
```

TODO: rewrite this with AbstractReaders

<a name="_io_ngram"/>

### Ngram file

TODO: rewrite this with AbstractReaders

```python
from trseeker.seqio.ngram_file import save_ngram_index

save_ngram_index(ngram_index_file,
		       hash2id,
		       result_tf,
                     	       result_df, 
                     	       result_rtf, 
                     	       result_rdf,
	                      seen_rev, 
	                      hash2rev,
                                    ngram2kmerann)

save_ngram_pos_index(ngram_trids_file, id2trids, id2trid2tf)    

save_distance_data(dist_file, distances)
```

<a name="_io_blast"/>

### Blast results file

```python	
from trseeker.seqio.blast_file import get_blast_result

get_blast_result(blast_file, length, gap_size=1000, min_align=500, min_length=2400, format_function=None, min_score=90)
```

Для фильтров используются following parameters:

1. **cutoff** - TRs array length.
2. **match** and **min_match** from **blast_match_proc** settings
3. **min_align** is minimal length of alignment (array_length / .min_alig)
4. **blast_gap_size** is size of the gap between two match regions (array_length * min_match)
5. **min_length** is minimal total length (array_length * match)

Этапы анализа:

1. From blast output files extracted ref_gi to BlastResultModel dictionary.
2. Removed all nested matches.
3. Overlapping matches joined.
4. Joined matches with length less than **min_align** are dropped.
5. Matches with gap length less than **gap_size** are joined
6. Joined gapped matches with length less than **min_align** are dropped.
7. Results formatted with given format_function


Dataset updating functions:

```python
update_with_repbase_blast_result(trs_dataset, annotation_self_folder, filters)

# inside function
# result is semicolon-delimited
trs_dataset[i].trf_repbase = result

# filters example
filters = {
	"blast_gap_size": 300,
	"min_align": 1000,
	"min_length": 3000,
}
```

```python
update_with_self_blast_result(trs_dataset, annotation_self_folder, filters_obj)
	
# inside function
# result comma-delimited
# inside get_filters(array_length) generate parameters by array_length:
# 	filters = filters_obj.get_filters(trf_obj.trf_array_length)
trs_dataset[i].trf_family_self = result
```

```python
update_with_ref_blast_result(trs_dataset, annotation_self_folder, filters)

# inside function
# comma-delimited data
trs_dataset[i].trf_family_ref = result
```

TODO: refractor to more generic form with setattr(...)

Self and ref annotation comma-delimited. Repbase annotation semicolon-delimited. Avalibale format functions:

- format_function_repbase(blast_dataset)
- format_function_self(blast_dataset) [DEFAULT]

While blast results parsing:
- skipped all inner matches
- joined overlapped
- skipped all alignments with length less thatn min_align paramter
- joined gapped if gap less than gap_size paramter

<a name="_io_sra"/>

### SRA file

TODO: move FastqObj to models

```python
from trseeker.seqio.sra_file import FastqObj

fastq_obj = FastqObj(head, seq, strain, qual_str, phred33=False)
print fastq_obj.fastq
print fastq_obj.fastq_with_error(error)
print fastq_obj.sequence
print fastq_obj.fasta
print fastq_obj.trimmed_fastq
print fastq_obj.gc
```

Properties:

- head
- seq
- qual
- strain
- id
- trimmed_seq
- trimmed_qual
- qv: list of Q - phredX
- adaptor_positions: positions of adapter
- adapter_contamination
- qual_junk
- status
- parts

Methods related to trimming:

- trim(), NotImplemented
- trim_by_quality_cutoff(cutoff=30)
- trim_exact_adaptor(adaptor, verbose=False), NotImplemented
- trim_inexact_adaptor(adaptor, verbose=False, num_errors=2), NotImplemented

Additional functions:

- fastq_reader

```python
for fastq_obj in fastq_reader(fastq_file, phred33=False):
	print fastq_obj.seq
```


<a name="_tools"/>

## Toolkit

<a name="_annotation"/>

### Assembly annotation

Compute intersection between two datasets.
Dataset format: list of tuples (chrName, startPos, endPos, feature_id).

```python
from trseeker.tools.annotation import *

compute_intersection_intervals(data_a, data_b)
```

Usage example:

```python
data_b = open("cat_good_trf.gff3").readlines()
data_b = [x.strip().split("\t") for x in data_b]
data_b = [(x[0], int(x[3]), int(x[4]), x[-1]) for x in data_b]

for x in compute_intersection_intervals(data_a, data_b):
    print x
```

<a name="_assembly_tools"/>

### Assembly tools

```python
from trseeker.tools.assembly_tools import *

N50_contig_length, N50, shortest_seq, longest_seq = get_N50(lengths)
```

<a name="_tools_tulip"/>

### Tulip files

Format and write network data to tulip format from distance_file with given cutoff (defaul=90).
    
- distance_file: (node_a, node_b, distance)
- tulip_file_output: network in tulip format
- cutoff: distance cutoff

```python
from trseeker.tools.tulip_tools import *

format_distances_to_tulip(distance_file, tulip_file_output, cutoff=90)
```

<a name="_tools_trs_nomen"/>

### TRs nomenclature

TODO: rename package
```python
from trseeker.tools.trs_groups import *
```
Get next family index. Index limited to latin alphabet characters. Otherwise will return X1,X2 and so on:
```python
letter = get_index(i)
```
Get most frequent element of list:
```python
element = get_popular(s)
```
Get unit and letter for family. Unit picked by most common pmatch value. A seen_units is a dictionary unit_size to frequence, it is needed for letter choosing:
```python
unit, letter, seen_units = get_family_name(ids, seen_units)	
```
Join families with common members, sort by family size:
```python
join_families_with_common(families)
```

<a name="_tools_sequence"/>

### Sequence tools

```python
from trseeker.tools.sequence_tools import *
```
Return complementary sequence:
```python
revcom = get_revcomp(sequence)
```

Get shift variants for TR monomer:

```python
sequences = get_revcomp(sequence)
```

Return normalized sequence with following rules:

1. if T > A then return reverse complement
2. if A == T and G > C  then return reverse complement
3. else return sequence

```python
sequence = fix_strand(sequence)
```
Count GC content:
```python
gc = get_gc(sequence)
```
Check for N|n in sequence. Return n(with N) or w(whole):
```python
bool_value = check_gapped(sequence)
```
Return subsequence:
```python
get_subseq(seq, start, end)
```
Clear sequence (ACTGN alphabet):

1. lower case
2. remove any gaps
3. remove all letters except atgcn

```python
sequence = clear_sequence(sequence)
```
Return list of sequences splited by pattern with added end. Pattern is regexp:
```python
restriction(sequence, pattern, end=None)
```
Return number of fragments after restriction of sequence with given regexp pattern:
```python
restriction_fragments_n(sequence, pattern)
```
Check tandem repeat synonims  between two sequence with same length:
```python
check_cyclic_repeats(seq_a, seq_b)
```
Return sequence with n mutations:
```python
random_mutation(seq, n, alphabet="actgn +")
```
Return consensus string for given list of strings:
```python
get_consensus(strs)
```

Take a minimal sequence from lexicographically sorted rotations of sequence and its reverse complement.
Example: ACT, underlined - reverse complement seqeunces
ACT, AGT, CTA, GTA, TAC, TAG.

Find all possible multimers, e.g. replace GTAGTAGTA consensus sequence with ACT.

Return:

- list of sorted TRs
- list of (df, consensus) pairs sorted by array number

```python
remove_consensus_redundancy(trf_objs)
```

<a name="_tools_var"/>

### Various useful functions for working with files

```python
from trseeker.tools.seqfile import *
```

- save_list(file_name, data)
- save_dict(file_name, dictionary)
- save_sorted_dict(d, file, by_value=True, reverse=True, min_value=None, key_type=None)
- count_lines(file)
- sort_file_by_int_field(file_name, field)

<a name="_tools_other"/>

### Other useful functions (other_tools.py)

```python
from trseeker.tools.other_tools import *
```

Sort dictionary by key. Retrun list of (k, v) pairs:
```python
sort_dictionary_by_key(d, reverse=False)
```
Cast list elements to string:
```python
as_strings(alist)
```
Return list from dictionary, return list of (key, value) pairs:
```python
dict2list(adict)
```
Sort dictionary by value. Retrun list of (v, k) pairs:
```python
sort_dictionary_by_value(d, reverse=False)
```
Remove duplicates from list and sort it:
```python
remove_duplicates_and_sort(data)

remove_duplicates(data)
```
Remove redundancy or empty elements in given list:
```python
remove_redundancy(alist)
```
Remove nested fragments, ata format [start, end, ...]:
```python
clear_fragments_redundancy(data, extend=False, same_case_func=None)
```




<a name="_tools_kmers"/>

### Ngrams (kmers) tools

```python
from trseeker.tools.ngrams_tools import *
```

Yields all ngrams of length n from given text:
```python
generate_ngrams(text, n=12)
```
Yields all ngrams of length n from given text:
```python
generate_window(text, n=12, step=None)
```
Returns m most frequent (ngram of length n, tf) tuples for given text:
```python
get_ngrams(text, m=None, n=23, k=None, skip_n=False)
```
Returns m most frequent (ngram of length n, fraction of possible ngrams) tuples for given text:
```python
get_ngrams_freq(text, m=None, n=23, k=None)
```
Returns a feature set {'ngram':'ngram',...}  of m most frequent ngram of length n for given text:
```python
get_ngrams_feature_set(text, m=5, n=12)
```
Returns a distance between two ngram sets where distance is a sum(min(ngram_a, ngram_b) for each common ngram). Format for ngrams_a and ngrams_b is a dictionary {ngram:n, ...}:
```python
get_ngram_freq_distance(ngrams_a, ngrams_b)
```
Returns a distance between two ngram sets where distance is a len(). Format for ngrams_a and ngrams_b is a dictionary {ngram:n, ...}:
```python
get_ngram_common_distance(ngrams_a, ngrams_b)
```
Return repeatness coefficient. From 0 (e.g. polyA) to 1 (unique sequence):
```python
get_repeatness_coefficent(length, k, kmern)
# Inside this function:
# kmern * (N-1) / (N**2)
```
Return expressivenesss coefficient.
```python
get_expessiveness_coefficent(kmern, k)
# Inside this function:
# kmern * 1. / 4^k
```
Update tf and df data with k-mers from given sequence (it will be lowered):
```python	
from trseeker.tools.ngrams_tools import count_kmer_tfdf

k = 23
for sequence in sequences:
	tf_dict, df_dict, local_tf, local_df = count_kmer_tfdf(sequence, tf_dict, df_dict, k)
```

Compute max df, number and procent of sequence with given ngram. Return (maxdf, nmaxdf, pmaxdf, ngram_seqs)

```python
(maxdf, nmaxdf, pmaxdf, ngram_seqs) = get_df_stats_for_list(data, k, kmer2df):
```

Get tf and df frequencies for kmer list:

```python
get_kmer_tf_df_for_data(data, k, docids=False)
```

Get list of (kmer, revkmer, tf, df, docids) for given data:

```python
process_list_to_kmer_index(data, k, docids=True, cutoff=None)
```

Get list of (kmer, revkmer, tf, df, docids) for multifasta file:

```python
compute_kmer_index_for_fasta_file(file_name, index_file, k=23)
```

Get list of (kmer, revkmer, tf, df, docids) for TRs file, filter ny sequence complexity defined by get_zlib_complexity function:

```python
compute_kmer_index_for_trf_file(file_name, index_file, k=23, max_complexity=None, min_complexity=None)
```

Compute kmer coverage and set of kmers:

```python
(m, kmers) = get_sequence_kmer_coverage(sequence, kmers, k)
```

<a name="_tools_ed"/>

### Edit distance functions

```python
from trseeker.tools.edit_distance import *
```

Return list of element (*d*) positions in given list:
```python
get_pos(L, d)
```
Return edit distance between two given strings. Edit distance dynamic programming implementation. Length of second sequence does not change similarity value:

1. monomer mode: double short monomer sequence if len2/len1 < 2
2. full_info: boolean, return full information

Return: percent of ED similarity, or if full_info (distance, number of d, positions of d, S matrix)
```python
get_ed_similarity(s1, s2, monomer_mode=False, full_info=False, verbose=False)
```

Return two sequences edit distance full information.    

1. s1: sequence of tandem repeat monomer
2. s2: sequence of tandem repeat monomer
3. verbose: boolean
4. monomer mode: double short monomer sequence if len2/len1 < 2.

Return: ("max_sim Nmax_sim %sim pos_list", pos_list, all_result, length seq1, length seq2)

```python
get_edit_distance_info(s1, s2, verbose=False, monomer_mode=False)
```

Return ED valie for two sequences:
```python
get_edit_distance(s1, s2)
```

Get last row of ED matrix between two strings:
```python
get_edit_distance_row(s1, s2)
```

Get Hamming distance: the number of corresponding symbols that differs in given strings:
```python
hamming_distance(s1, s2)
```

<a name="_tools_repbase"/>

### Working with Repbase files

TODO: move to Readers

```python
from trseeker.tools.repbase_tools import join_repbase_files
```
Function join all Repbase fasta files in one huge fasta; reformat headers for compatibility with NCBI tools:
```python
join_repbase_files(input_folder, output_file)
```

<a name="_tools_trace"/>

### Working with Trace files

```python
from trseeker.tools.trace_tools import unclip_trace_file
```

Unclip. Remove vector flanks from Trace data.
```python
unclip_trace_file(fasta_file, clip_file, uncliped_file)
```

Read clip dictionary from clip_file_name gz archive:
```python
get_clip_data(clip_file_name)
```


<a name="_tools_stat"/>

### Function related to statistics

```python
from trseeker.tools.statistics import *
```

Calculate sigma, return sum(module(xi - mean):
```python
get_sigma(data)
```
Return dictionary containing simple statistics for given list. Dictionary keys: mean, variance, sigma, sample deviation:
```python
get_mean(data)

get_standard_deviation(variance)

t_test(sample_mean, dist_mean, variance, N)

get_element_frequences(data)

get_simple_statistics(data)
```

<a name="_tools_parsers"/>

### Parsers

```python
from trseeker.tools.parsers import *

parse_fasta_head(fa_head)

parse_chromosome_name(head)

trf_parse_line(line)

trf_parse_param(line)

trf_parse_head(line)

get_wgs_prefix_from_ref(ref)

get_wgs_prefix_from_head(head)
```

<a name="_tools_trf"/>

### Functions related to TRF
```python
from trseeker.tools.trf_tools import *
```

```python
trf_search(file_name)
```
```python
trf_search_in_dir(folder, verbose=True, file_suffix=".fa", output_folder=None)
```

Parallel TRF seeach:

```python
trf_worker(args)
trf_search_in_dir_parallel(folder, verbose=True, file_suffix=".fa", output_folder=None, threads=1)
```

Create output TRF file with tandem repeats with length greater than from input file. Function returns number of tandem repeats in output file:

```python
trf_filter_by_array_length(trf_file, output_file, cutoff)
```
Create output TRF file with tandem repeats with unit length greater than from input file. Function returns number of tandem repeats in output file:

```python
trf_filter_by_monomer_length(trf_file, output_file, cutoff)
```
Create output TRF file with tandem repeats with GI that don't match GI_LIST. List of GI, see TRF and FA specifications, GI is first value in TRF row.

```python
trf_filter_exclude_by_gi_list(trf_file, output_file, gi_list_to_exclude)
```
Write TRF file tab delimited representation.Representation: numerical|index|agc_apm|with_monomer|family
```python
trf_representation(trf_file, trf_output, representation)
```
Write statistics data: field, N:
```python
trf_write_field_n_data(trf_file, file_output, field, field_format="%s")
```
Write statistics data: field_a, field_b:
```python
trf_write_two_field_data(trf_file, file_output, field_a, field_b)
```
Function prints chr, all trs, 3000 trs, 10000 trs:
```python
count_trs_per_chrs(all_trf_file)
```
Function prints number of items with given fasta head fragment:
```python
count_trf_subset_by_head(trf_file, head_value)
```
Several fasta heads impossible to parse, so it is simpler to fix them postfactum:
```python
fix_chr_names(trf_file, temp_file_name=None, case=None)
```


<a name="_tools_trs"/>

### Working with TRs datasets
```python
from trseeker.tools.trs_dataset import *
```

- COLORS dictionary
- COLORS_50 dictionary

```python
COLORS_50 = {
  "denim": ("#1560bd", (21,96,189)),
}
# TODO: finish it
```

- get_colors(family)
- get_colors_rgb(family)
- create_mathematice_dataset_by_family(trs_dataset, path_to_mathematica_folder, min_borders, max_borders)

<a name="_tools_sra"/>

### Working with SRA data
```python
from trseeker.tools.sra_tools import *
```

Iterate over fastaq data.
```python
sra_fastaq_reader(file_name)
```
Iterate over fasta SRA data.

```python
sra_fasta_reader(file_name)
```

- read_fastaq_freq_to_memory(file_name)
- fastaq_to_fasta(file_name, output_file)
- write_fastaq_repeats(input_file, output_file, min_tf=1000)
- seq_to_bin(seq)
- bin_to_seq(bseq)

Read SRA fasta data without read repeats.
```python
write_reduced_fasta(input_file, output_reduced)
```
Write ngrams data from SRA fasta data.
```python
write_ngrams(input_file, output_ngram, NGRAM_N)
```


<a name="_tools_patterns"/>

### Working with sequence patterns
```python
from trseeker.tools.sequence_patterns import *
```

Avaliable patterns:

- MARS1
- MARS2
- CENPB
- PRDB9

Functions:

Translate sequence in DNA15 in reg exp:
```python
re_translate(sequence)
```

Return list of sequences where one letter changed to N:
```python
re_get_plus_minus(sequence)
re_get_mutations(sequence)
```

Return list of patterns one not mutated and second mutated:
```python
get_double_pattern(pattern_static, pattern_dynamic)
```

- get_mutated_pattern_twice(pattern)
- get_mutated_pattern_trice(pattern)
- pattern_search(name, sequence, pattern_function, pattern_function_params)

<a name="_tools_blast"/>

### Working with BLAST
```python
from trseeker.tools.blast_tools import *
```

- blastn(database, query, output, e_value=None)
- create_db(fasta_file, output, verbose=False, title=None)
- alias_tool(dblist, output, title)
- bl2seq(input1, input2, output)
- create_db_for_genome(file_pattern=None, chromosome_list=None, output=None, title=None)
- get_gi_list(gi_score_file, score_limit=90)
- get_all_gi_from_blast(blast_file, mode="gi")
- get_all_blast_obj_from_blast(blast_file, mode="ref")
- blastn_search_for_trs(trf_large_file, db, annotation_self_folder, temp_file, skip_by_family=None, is_huge_alpha=False)
- bl2seq_search_for_trs(trf_large_file, annotation_bl2seq_folder, temp_file)

<a name="_tools_sa"/>

### Working with suffix arrays
```python
from trseeker.tools.sa_tools import *
```

Fasta to SA input. Text file of sequences delimeted by $ symbol:

```python
fasta_to_sa_input(fasta_file, sa_input_file, index_file_name=None, file_name=None, start_id=None, increment=False)
```

SA input file to fasta file:

```python
sa_input_to_fasta(input_file, output_file)
```

Function precompiles doc2id and id2doc pickled dictionary:

```python
pickle_dictionary_for_docid_trid(sa_doc_index, doc_to_trf_file, trf_to_doc_file)
```

Function filters and writes sa data with given min tf and df:

```python
filter_sa_dataset(sa_file, output_file, min_tf, min_df)
```

Yields corpus texts. Corpus is $ delimited text file:

```python
iterate_sa_corpus(corpus)
```

<a name="_tools_ncbi_genome"/>

### Working with NCBI genome information
```python
from trseeker.tools.ncbi_genomes import *
```

Scrap finished WGS genome projects wtih given SubGroup name:
```python
load_genome_projects(subgroup)
```
Scrap finished WGS genome projects wtih given SubGroup name:
```python
load_genome_projects_exclude(notsubgroups)
```
Print data for initiating PySatDNA projects:
```python
print_add_project_data(genome_projects, pid_type, pr_type)
```

<a name="_tools_networks"/>

### Working with graph data
```python
from trseeker.tools.network_tools import *
```

Compute network slices using graph_tool library or networkx.
```python
compute_network(network_file, output_file_pattern, id2nodename)
```

- create_graph_on_fly(network_file)
- load_graph(ml_file)
- analyse_network_graph_tool(G, output_file_pattern, trf_large_index_file)
- write_classification_graph_tool(output_file, components, id2nodename)
- write_classification_neworkx(output_file, components, trid2meta)
- create_graphml(network_file, ml_file, id2nodename)
- load_networkx(network_file)
- init_graph_networkx(network_data, start=0, precise=1, id2nodename=None)
- analyse_networkx(G, network_data, output_file_pattern, id2nodename)

<a name="_tools_ncbi_ann"/>

### Working with NCBI annotations
```python
from trseeker.tools.ncbi_annotation_tools import *
```

Read ideogram file and return return dict chr -> list.

```python
get_ideogram_dict(idiogram_file, mode="NCBI")

get_gene_list(file_gene_list, color='#000000', left_padding=30, gene_group_label='C57BL/6J')
# NB: Not implemented yet.
``` 

<a name="_tools_jellyfish"/>

### Jellyfish wrapper

```python
from trseeker.tools.jellyfish_tools import *

# jellyfish settings:
location = settings["blast_settings"]["jellyfish_location"]
```

 Count kmers with Jellyfish

```python
count_kmers(input_fasta, ouput_prefix, k, mintf=None)

merge_kmers(folder, ouput_prefix, ouput_file)

stats_kmers(db_file, stats_file)

histo_kmers(db_file, histo_file)

dump_kmers(db_file, fasta_file, dumpmintf)

query_kmers(db_file, query_hashes, both_strands=True, verbose=True)

get_kmer_db_and_fasta(folder, input_file, kmers_file, k=23, mintf=None)

query_and_write_coverage_histogram(db_file, query_sequence, output_file, k=23)
# Save coverage histogram into output_file for given query_sequence.
```

Shortcut:

```python
from trseeker.tools.jellyfish_tools import sc_compute_kmer_data

sc_compute_kmer_data(fasta_file, jellyfish_data_folder, jf_db, jf_dat, k, mintf, dumpmintf)
```

<a name="_trs_types"/>

### Classifiction TRs in types

Classification of TRs into types.

```python
from trseeker.tools.trs_types import *
```

Available types:

```python
trs_types = {
    "all": 0,
    "micro": 0,
    "perfect": 0,
    "too_short": 0,
    "100bp": 0,
    "gc": 0,
    "x4": 0,
    "entropy": 0,
    "ok": 0,
}
```

Usage:

```python
trs_types, class_objs, trf_objs = get_trs_types(trf_all_file, trf_all_class_file, settings)
```

Used settings:

```python
settings["other"]["trs_types"]["possible_array_length"]
settings["other"]["trs_types"]["microsat_monomer"]
settings["other"]["trs_types"]["min_array_length"]
settings["other"]["trs_types"]["min_gc"]
settings["other"]["trs_types"]["max_gc"]
settings["other"]["trs_types"]["min_copy_number"]
settings["other"]["trs_types"]["min_entropy"]
```

<a name="_tools_tree"/>

### Tree data

```python
from trseeker.tools.tree_tools import *
```

#### NodeData class

Stores tree-relevant data associated with nodes.

```python
node = NodeData(taxon=None, branchlength=1.0,
                 support=None, comment=None, bootstrap=None,
                 taxon_name=None, distance=None, elements=None)
```

#### Chain class

Stores a list of nodes that are linked together.

```python
chain = Chain()

# Return a list of all node ids
ids = chain.all_ids()

# Attaches node to another: (self, node, prev)
chain.add(node, prev=None)

# Deletes node from chain and relinks successors to predecessor: collapse(self, id).
chain.collapse(id)

# Kills a node from chain without caring to what it is connected
chain.kill(id)

# Disconnects node from his predecessor
chain.unlink(id)

# Connects son to parent
chain.link(id)

# Check if grandchild is a subnode of parent
chain.is_parent_of(parent, grandchild)

# Returns a list of all node_ids between two nodes (excluding start, including end)
chain.trace(start, finish)
```

#### Node datastructure

A single node

```python
# Represents a node with one predecessor and multiple successors
node = Node(data=None)
node.set_id(id)
node.set_data(data)
id = node.get_id()

# Returns a list of the node's successors
id = node.get_succ()
# Returns the id of the node's predecessor
id = node.get_prev()
node.set_prev()

# Adds a node id to the node's successors
node.add_succ(id)
# Removes a node id from the node's successors
node.remove_succ(id)
# Sets the node's successors
node.set_succ(id)
```

#### Tree datastructure

Represents a tree using a chain of nodes with on predecessor (=ancestor)
    and multiple successors (=subclades)

```python

t = Tree(weight=1.0, rooted=False,
                 name="", data=NodeData, max_support=1.0)

t.is_leaf(node)
t.is_subtree(node)
t.is_terminal(node)
t.is_internal(node)
t.is_preterminal(node)
t.has_support(node=None)

l = t.get_node_to_one_d(tree)
t.is_identical(tree)

node = t.node(id)

nodes = t.get_terminals()
n = t.count_terminals(node=None)
node = t.common_ancestor(node1, node2)
d = t.distance(nod1, node2)

dist_dict = t.get_distance_for_all_from_node(id)

bl = t.sum_branchlength(root=None, node=None)

ids = t.split(parent_id=None, n=2, branchlength=1.0)

# Bootstrap tree with given cutoff value. Nodes with suppert below cutoff will be collapsed
t.bootstrap(value)

# Prunes a terminal taxon from the tree
t.prune(taxon)

# Return a list of all otus downwards from a node
t.get_taxa(node_id=None)

ids = t.search_taxon(taxon)

t.to_string(support_as_branchlengths=False, branchlengths_only=False, plain=True, plain_newick=False, ladderize=None)

t.make_info_string(data, terminal=False)

t.ladderize_nodes(nodes, ladderize=None)

t.newickize(node, ladderize=None)
```

Functions:

Get list of (slice_id, file_path):

```python
slice_files = get_slice_files(folder)
```

read_tree_slices(slice_files):

```python
slice_files = read_tree_slices(slice_files)
```


<a name="_tools_entrez"/>

### Function for downloading datasets from NCBI

```python
from trseeker.tools.entrez_datanase import *
```

Download all proteins according to a given query into output file and then return these proteins as fasta text.

```python
fasta_data = download_proteins_from_ncbi(query, output_file, email, batch=500, verbose_step=1000)
```

Download items from given database according to query, rettype, and retmode. Save output to output_file and return as text.

```python
data = download_items_from_ncbi(query, 
                             database, 
                             output_file, 
                             email, 
                             rettype="fasta", 
                             retmode="text", 
                             batch=500, 
                             verbose_step=1000)

```

Get items from given database according to query, rettype, and retmode.

```python
items = get_items_from_ncbi(query, 
                         database, 
                         email, 
                         rettype="fasta", 
                         retmode="text", 
                         batch=500, 
                         verbose_step=1000)
```

Get RNA SRA datasets from NCBI according to taxid. 

Output includes LIBRARY_SOURCE, STUDY_ABSTRACT, DESIGN_DESCRIPTION, PRIMARY_ID, DESCRIPTION, LINKS.

And the second part of the output contains full xml data.

```python
data, _ = get_rna_sra_datasets_by_taxid(taxid, email, batch=500, verbose_step=1)

all_links = {}

for *item, links in datasets.values():
    links = dict([
                   (url.split("/")[-1], url)
                        for 
                            url 
                                in links])
    for sra in links:
        all_links[sra] = links[sra]
```

Download and unpack SRA filrs from NCBI according to taxid.

```python
download_rna_sra_datasets_by_taxid(taxid, email, output_folder, threads=30, batch=500, verbose_step=1)
```

Download genomes and annotation from NCBI according to taxid.
Return refseq and genbank datasets.

```python
download_genome_assemblies_and_annotation_from_ncbi(taxid, output_folder, threads=30, only_refseq=True)
```



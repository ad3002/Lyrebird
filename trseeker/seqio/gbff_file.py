#!/usr/bin/env python# -*- coding: utf-8 -*-
#
#@created: 05.06.2011
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com 
"""
Classes:
    
- FastaFileIO(AbstractBlockFileIO)
    
"""
import os
import re
from trseeker.seqio.block_file import AbstractBlockFileIO
from trseeker.models.sequence_model import SequenceModel
from trseeker.tools.sequence_tools import clear_sequence
from trseeker.models.genbank import GenbankData
from PyExp import sc_iter_filename_folder


class GbffFileIO(AbstractBlockFileIO):
    """ Working with multi fasta files, where each block starts with '>' token.
        
    Overrided public methods:
    
    - __init__(self)
    
    Inherited public properties:
    
    - data  - iterable data, each item is tuple (head, body)
    - N     - a number of items in data
    
    Inherited public methods:
    
    - read_from_file(self, input_file)
    - read_online(self, input_file) ~> item
    - get_block_sequence(self, head_start, next_head, fh)
    - get_blocks(self, token, fh)
    - gen_block_sequences(self, token, fh)
    - read_from_db(self, db_cursor)
    - write_to_file(self, output_file)
    - write_to_db(self, db_cursor)
    - read_as_iter(self, source)
    - iterate(self) ~> item of data 
    - do(self, cf, args) -> result
    - process(self, cf, args)
    - clear(self)
    - do_with_iter(self, cf, args) -> [result,]
    - process_with_iter(self, cf, args)
        
    """

    def __init__(self):
        """ Overrided. Hardcoded start token."""
        token = "LOCUS       "
        super(GbffFileIO, self).__init__(token)


def get_section(fh, line, keyword):
    text = line.replace(keyword, "").strip()
    while True:
        line = fh.readline()
        if line.startswith(" "):
            text += " " + line.strip()
        else:
            break
    return text, line


def sc_iter_gbff(file_name):
    """ Iter over gbff file."""

    HEAD = {}
    META = {}
    FEATURES = []
    SEQ = []
    
    seq_obj = None
    
    with open(file_name) as fh:
        line = fh.readline()
        while True:
            if not line:
                break
            line = line.rstrip()
            if line == "//":
                line = fh.readline()
                continue
            if line.startswith("LOCUS"):
                if seq_obj:
                    yield seq_obj
                seq_obj = GenbankData()
                seq_obj.parse_locus(line)
                line = fh.readline()
                continue
            if line.startswith("DBSOURCE"):
                text, line = get_section(fh, line, "DBSOURCE")
                seq_obj.gb_references.append("DBSOURCE:%s" % text)
                continue
            if line.startswith("DEFINITION"):
                seq_obj.gb_definition, line = get_section(fh, line, "DEFINITION")
                continue
            if line.startswith("ACCESSION"):
                text, line = get_section(fh, line, "ACCESSION")
                seq_obj.gb_accessions = [x for x in text.split() if x]
                continue
            if line.startswith("VERSION"):
                text, line = get_section(fh, line, "VERSION")
                seq_obj.gb_versions = [x for x in text.split() if x]
                continue
            if line.startswith("KEYWORDS"):
                text, line = get_section(fh, line, "KEYWORDS")
                seq_obj.gb_keywords = [x for x in text.split() if x]
                continue
            if line.startswith("SOURCE"):
                seq_obj.gb_source, line = get_section(fh, line, "SOURCE")
                continue
            if line.startswith("ORGANISM"):
                seq_obj.gb_organism, line = get_section(fh, line, "ORGANISM")
                continue
            if line.startswith("DBLINK"):
                text, line = get_section(fh, line, "DBLINK")
                seq_obj.gb_dblink.append(text)
                continue
            if line.startswith("REFERENCE"):
                text, line = get_section(fh, line, "REFERENCE")
                seq_obj.gb_references.append(text)
                continue        
            if line.startswith("COMMENT"):
                text, line = get_section(fh, line, "COMMENT")
                seq_obj.gb_references.append(text)
                continue
            if line.startswith("FEATURES"):
                feature_text = []
                while True:
                    line = fh.readline()
                    if line.strip().startswith("CONTIG"):
                        text, line = get_section(fh, line, "CONTIG")
                        seq_obj.gb_contig = text
                        break
                    if line.startswith("//") or line.startswith("ORIGIN"):
                        break
                    feature_text.append(line)
                features = "".join(feature_text)
                features = features.replace("\n", "^")
                features = re.sub("\^ {21}[^/]", "", features)
                features = re.sub("\^ {21}/", "\t", features)
                features = re.sub("\^\s{5}", "\n", features)
                # ['CDS             complement(240314..241414)', 'gene="dnaN"', 'locus_tag="PI91_RS21880"', 'old_locus_tag="PI91_21630"', 'EC_number="2.7.7.7"', 'inference="COORDINATES: similar to AAequence:RefSeq:WP_006177590.1"', 'note="Derived by automated computational analysis usingene prediction method: Protein Homology."', 'codon_start=1', 'transl_table=11', 'product="DNA polymerase III subunit beta"', 'protein_id="WP_008500170.1"', 'translation="MKFTVEREHLLKPLQQVSGPLGGRPTLPILGNLLLQVADGTLSLGTDLEMEMIARVTLTQPHDAGATTVPARKFFDICRGLPEGAEIAVQLEGDRMLVRSGSRFSLSTLPAADFPNLDDWQSEVEFTLPQATMKRLIEATQFSMAHQDVRYYLNGMLFTEGEELRTVATDGHRLAVCSMPIGDSLPNHSVIVPRKGVIELMRMLDGGDTPLRVQISNNIRAHVGDFVFTSKLVDGRFPDYRRVLPKNPDKTLEAGCDSLKQAFARAAILSNEFRGVRLYVSENQIKITANNPEQEEAEEILDVTYAGAEMEIGFNVSYVLDVLNALKCEVRILLTDSVSSVQIEDAASQSAAYVVMPMRL"']
                for feature in features.split("\n"):
                    F = {}
                    feature = feature.strip()
                    if not feature:
                        continue
                    d = feature.split("\t")
                    name,coord = [x for x in d[0].split() if x]
                    F["ftype"]= name
                    if "complement" in coord:
                        F["fstrand"] = "-"
                    else:
                        F["fstrand"] = "+"
                    F["frawcoord"] = coord
                    if "join" in coord:
                        F["fcomplex"] = 1
                    else:
                        F["fcomplex"] = 0
                    if ">" in coord or "<" in coord:
                        F["fincomplete"] = 1
                    else:
                        F["fincomplete"] = 0
                    poses = re.findall("(\d+)", coord)
                    poses = list(map(int, poses))
                    if len(poses) == 2:
                        F["fstart"] = poses[0]
                        F["fend"] = poses[1]
                    else:
                        F["fstart"] = min(poses)
                        F["fend"] = max(poses)
                    for f in d[1:]:
                        f = f.split("=")
                        if len(f) == 2:
                            name, value = f
                            value = value.replace('"',"")
                        else:
                            name = f[0]
                            value = 1
                        F[name] = value
                    seq_obj.gb_features.append(F)
                continue
            if line.startswith("ORIGIN"):
                while True:
                    line = fh.readline()
                    if line.startswith("//"):
                        break
                continue
            print("Unexpected line:")   
            print(file_name)
            print(line)
            input("?")
    yield seq_obj

def sc_iter_gbff_simple(file_name):
    """ Iter over gbff file."""

    reader = GbffFileIO()
    for seq_obj in  sc_iter_gbff(file_name):
        yield seq_obj.seq_gi, seq_obj.sequence


def sc_covert_gbff_in_folder_to_fasta(gbbf_folder, fasta_folder, fasta_postfix, mask):
    ''' Convert all gbff files in given gbbf_folder according to mask, write result
    to fasta_folder with corresponding fasta_postfix.
    '''
    for file_name in sc_iter_filename_folder(gbbf_folder, mask):

        print("Convert file: %s" % file_name)

        fasta_name = ".".join( file_name.split(".")[:-1] )

        file_name = os.path.join(gbbf_folder, file_name)

        fasta_output = os.path.join(fasta_folder,
                                    "%s.%s" % (fasta_name,fasta_postfix)
                                    )
        print(file_name, fasta_output)
        if os.path.isfile(fasta_output):
            os.remove(fasta_output)
        with open(fasta_output, "ab") as fw:
            for seq_obj in sc_iter_gbff(file_name):
                fw.write(seq_obj.ncbi_fasta)
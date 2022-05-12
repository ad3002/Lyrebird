#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 10.01.2016
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com

from PyExp import AbstractModel

class MafBlock(AbstractModel):

    dumpable_attributes = [
        "full_name",
        "taxon",
        "name",
        "start",
        "end",
        "length",
        "total_length",
        "strand",
        "sequence",
    ]

    int_attributes = [
        "start",
        "end",
        "length",
        "total_length",
    ]

    def __init__(self, line):
        super(MafBlock, self).__init__()
        data = [x for x in line.strip().split() if x]
        self.full_name = data[1]
        self.taxon = data[1].split(".")[0]
        self.name = ".".join(data[1].split(".")[1:])
        self.start = int(data[2]) + 1
        self.length = int(data[3])
        self.end = self.start + self.length - 1
        self.strand = data[4]
        self.total_length = int(data[5])
        self.sequence = data[6]


class MafSection(AbstractModel):

    dumpable_attributes = [
        "n",
    ]
    int_attributes = [
        "n",
    ]

    def __init__(self):
        super(MafSection, self).__init__()
        self.n = 0
        self.blocks = []

    def add(self, block):
        self.n += 1
        self.blocks.append(block)

    def __str__(self):
        s = []
        for i, block in enumerate(self.blocks):
            s.append("\t".join(map(str, (i, block.sequence, block.full_name, block.start, block.end))))
        return "\n".join(s)

    def show(self):
        for i, block in enumerate(self.blocks):
            if "Anc" in block.taxon:
                continue
            print i, block.sequence, block.full_name, block.start, block.end, block.strand


def iter_maf_file(file_name):

    section = None
    with open(file_name) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            if line.startswith("a"):
                if section:
                    yield section
                section = MafSection()
            if line.startswith("s"):
                section.add(MafBlock(line))
        if section:
            yield section

if __name__ == '__main__':
    
    file_name = "/mnt/peru/akomissarov/CHCHSCF_2015-01-14.maf"

    human_names = set()

    last_start = {}
    last_end = {}
    last_name = None

    sequences = {}

    for section in iter_maf_file(file_name):
        # section.show()
        for block in section.blocks:
            if block.taxon == "AcinonyxJubatus":
                if last_name is None:
                    last_name = block.name
                    last_end.setdefault(last_name, 0)
                    last_start.setdefault(last_name, 0)
                    sequences.setdefault(last_name, [])

                if block.name != last_name:
                    last_name = block.name
                    last_end.setdefault(last_name, 0)
                    last_start.setdefault(last_name, 0)
                    sequences.setdefault(last_name, [])

                dstart = block.start - last_end[last_name]
                if dstart < 1:
                    continue
                if last_end[last_name] and dstart > 1:
                    sequences[last_name].append("N"*(dstart-1))
                sequences[last_name].append(block.sequence)

                print block.name, "%s -> %s-%s [%s]" % (last_end[last_name], block.start, block.end, dstart)

                # if dstart < 1 or dstart > 1000:
                #     print block.sequence, "%s -> %s-%s [%s]" % (last_end, block.start, block.end, dstart)
                #     raw_input("m?")

                last_start[last_name] = block.start
                last_end[last_name] = block.end

    for name, seqs in sequences.items():
        with open("/home/akomissarov/Dropbox/PySatDNA/maf_cheetah_%s.seq" % name, "w") as fw:
            fw.write(">%s\n%s" % (name, "".join(seqs)))
        print seqs, name
        raw_input("n?")
    # print block.sequence, "%s -> %s-%s [%s]" % (last_end, block.start, block.end, dstart)

        # raw_input("Next?")
    # print human_names

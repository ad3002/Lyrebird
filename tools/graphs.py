#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 21.05.2014
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com
"""
Repeat can be revealed as:
1) set of kmers forming it;
2) de Brujin graph from these kmers.

Starting from the most frequent (conservative) kmer,
we extend it to left and to right
getting parents and child kmers by shifting on one nucleotide.

"""
import networkx as nx
from collections import defaultdict
from networkx.algorithms import cycles
from PyExp import Timer


def color_function(x):
    if x > 0.9:
        color = "green"
    elif x > 0.3:
        color = "blue"
    else:
        color = "red"
    return color


class Graph(object):

    def __init__(self):
        self.graph = nx.DiGraph()
        self.nodes = set()
        self.edges = set()
        self.kmer2node = {}
        self.node2kmer = {}
        self.i = 0
        self.edge2color = {}
        self.edge2tf = {}
        self.edge2ptf = {}
        self.edge2fos = {}
        self.short2node = {}


    def collapse_node():
        pass


    def get_node_name(self, kmer, fos):
        if kmer in self.kmer2node:
            _kmer = self.kmer2node[kmer]
        else:
            _kmer = "%s_%s_%s" % (self.i, fos, kmer)
            self.kmer2node[kmer] = _kmer
            self.node2kmer[_kmer] = kmer
            self.i += 1
        return _kmer


    def add_node(self, kmer):
        if not kmer in self.nodes:
            self.graph.add_node(kmer)
            self.nodes.add(kmer)

    def add_edge(self, parent, child, color, ptf, tf, fos):
        if not (parent, child) in self.edges:
            self.graph.add_edge(parent, child, color=color, ptf=ptf, tf=tf, fos=fos)
            self.edge2color[(parent, child)] = color
            self.edges.add((parent, child),)
            self.edge2tf[(parent, child)] = tf
            self.edge2ptf[(parent, child)] = ptf
            self.edge2fos[(parent, child)] = fos


    def add_child(self, seed, kmer, ptf, tf, fos):

        color = color_function(ptf)
        
        kmer = self.get_node_name(kmer, fos)
        seed = self.get_node_name(seed, fos)

        self.add_node(kmer)
        self.add_node(seed)
        self.add_edge(seed, kmer, color, ptf, tf, fos)
    

    def check(self):
        pass

    def replace_node(self, node, nnode):

        parents = self.graph.predecessors(node)
        childs = self.graph.successors(node)
        for p in parents:
            
            self.graph.add_edge(p, 
                                nnode, 
                                color=self.edge2color[(p, node)],
                                fos=self.edge2fos[(p, node)],
                                tf=self.edge2tf[(p, node)],
                                ptf=self.edge2ptf[(p, node)],
                                )
            self.edge2color[(p, nnode)] = self.edge2color[(p, node)]
            self.edge2fos[(p, nnode)] = self.edge2fos[(p, node)]
            self.edge2tf[(p, nnode)] = self.edge2tf[(p, node)]
            self.edge2ptf[(p, nnode)] = self.edge2ptf[(p, node)]
            self.graph.remove_edge(p, node)

        for c in childs:

            self.graph.add_edge(
                                nnode, 
                                c, 
                                color=self.edge2color[(node, c)],
                                fos=self.edge2fos[(node, c)],
                                tf=self.edge2tf[(node, c)],
                                ptf=self.edge2ptf[(node, c)],
                                )
            self.edge2color[(nnode, c)] = self.edge2color[(node, c)]
            self.edge2fos[(nnode, c)] = self.edge2fos[(node, c)]
            self.edge2tf[(nnode, c)] = self.edge2tf[(node, c)]
            self.edge2ptf[(nnode, c)] = self.edge2ptf[(node, c)]

            self.graph.remove_edge(node, c)

        self.graph.remove_node(node)


    def simplify(self):

        translation = {}
        changed = set()

        # 1. translate all childs to nucleotids except terminal nodes
        for node in self.graph.nodes():          
            if node in translation:
                node = translation[node]
            nid, fos, seq = node.split("_")
            if len(seq) == 1:
                continue
            nucleotide = seq[-1]
            parents = self.graph.predecessors(node)
            if len(parents) > 0:
                if not nid in changed:
                    nnode = "%s_%s_%s" % (nid, fos, nucleotide)
                    changed.add(nid)
                    self.short2node[node] = nnode

                self.graph.add_node(nnode)
                self.replace_node(node, nnode)
                translation[node] = nnode
                

    def collapse(self):

        translation = {}

        for node in self.graph.nodes():
            while True:
                if node in translation:
                    node = translation[node]
                else:
                    break
            parents = self.graph.predecessors(node)
            childs = self.graph.successors(node)

            if len(parents) != 1 or len(childs) != 1:
                continue

            parent = parents[0]
            # print node, parents, childs
            # print "\t\t", parent,  self.graph.predecessors(parent), self.graph.successors(parent)


            if len(self.graph.successors(parent)) != 1:
                continue

            pnid, pfos, pseq = parent.split("_")
            nparent = parent + node.split("_")[-1]
            
            self.graph.add_node(nparent)
            # print self.graph.nodes()

            try:
                pp = self.graph.predecessors(parent)
            except:
                pp = []

            for p in pp:
                self.graph.add_edge(p, nparent, color=self.edge2color[(p, parent)])

                self.edge2color[(p, nparent)] = self.edge2color[(p, parent)]

                self.graph.remove_edge(p, parent)
            
            for c in childs:
                self.graph.add_edge(nparent, c, color=self.edge2color[(node, c)])

                self.edge2color[(nparent, c)] = self.edge2color[(node, c)]
                try:
                    self.graph.remove_edge(node, c)
                except Exception, e:
                    print e

            self.graph.remove_node(parent)
            self.graph.remove_node(node)
            translation[node] = nparent
            translation[parent] = nparent

            # print "\t", nparent,  self.graph.predecessors(nparent), self.graph.successors(nparent)
            # print self.graph.nodes()
            # print [(x, translation[x]) for x in translation if '37' in x]
            # print

    def remove_simple_bubbles(self):
        '''
        parent has two childred of equal lentgh that have same one child
        grandchild has only two parents

        '''
        for node in self.graph.nodes():
            childs = self.graph.successors(node)
            if len(childs) != 2:
                continue
            a, b = childs
            a_childs = self.graph.successors(a)
            b_childs = self.graph.successors(b)
            if len(a_childs) !=1 or a_childs != b_childs:
                continue
            if self.graph.successors(a_childs[0]) != 1:
                continue
            node = node + a + a_childs[0]



        pass


    def trim(self):

        beard = 1
        iteration = 0
        while beard:
            beard = 0
            iteration += 1
            for node in self.graph.nodes():
                parents = self.graph.predecessors(node)
                childs = self.graph.successors(node)
                if len(childs) == 0:
                    # trim node
                    for p in parents:
                        self.graph.remove_edge(p, node)
                    self.graph.remove_node(node)
                    beard += 1
                elif len(parents) == 0:
                    # trim node
                    for c in childs:
                        self.graph.remove_edge(node, c)
                    self.graph.remove_node(node)
                    beard += 1
            print "Iteration %s trimmed: %s" % (iteration, beard)


    def cycles(self, settings):
        with Timer("Compute cycles..."):
            cs = []
            for c in cycles.simple_cycles(self.graph):
                cs.append(c) 
        print "Cycles:", len(cs)
        raw_input("Continue?")
        for cycle in cs:
            yield cycle

    def draw(self, pref, fam):

        A = nx.to_agraph(self.graph)
        A.layout(prog='dot')
        A.draw('/home/akomissarov/Dropbox/PySatDNA/%s-%s.png' % (pref, fam))
        A.draw('/home/akomissarov/Dropbox/PySatDNA/%s-%s.svg' % (pref, fam))


    def save(self, step, settings):


        while True:
            was = len(self.graph.nodes())
            print "Before trimming %s" % len(self.graph.nodes())
            self.trim()
            print "After trimming %s" % len(self.graph.nodes())
            print "Before simplify %s" % len(self.graph.nodes())
            self.simplify()
            print "After simplify %s" % len(self.graph.nodes())
            

            break

            if was == len(self.graph.nodes()):
                break

        
        mode = "trimmed"

        # dot = write(self.graph)
        # gvv = gv.readstring(dot)
        # gv.layout(gvv,'dot')
        # gv.render(gvv,'png','/home/akomissarov/Dropbox/PySatDNA/europe-%s.png' % step)
        
        with open('/home/akomissarov/Dropbox/PySatDNA/%s-%s.txt' % (settings["prefix"], settings["family"]), "w") as fh:
            for (a,b),tf in self.edge2tf.items():
                fh.write("%s\t%s\t%s\n" % (a,b,tf))


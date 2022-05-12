#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 08.09.2011
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com 
'''
Functions for network construction.

TODO: check it.
'''
from trseeker.seqio.tr_file import get_all_trf_objs
from trseeker.seqio.tab_file import sc_iter_tab_file
from trseeker.models.trf_model import TRModel
import pickle
import sys
from trseeker.seqio.tr_file import read_trid2meta
import os
import array
from collections import defaultdict
try:
    import numpy
except Exception, e:
    print "Numpy import error: %s" % e
from PyExp import Timer

try:
    import graph_tool
    from graph_tool import topology
    from graph_tool.all import *
    GRAPHTOOL = True
except:
    print "WARNING: install graph_tool"
    import networkx
    GRAPHTOOL = False

def compute_network(network_file, output_file_pattern, id2nodename):
    ''' Compute network slices using graph_tool library or networkx.
    '''

    if GRAPHTOOL:
        print "Using graph_tool"
        ml_file = network_file + ".ml"
        with Timer("Create ml"):
            create_graphml(network_file, ml_file, id2nodename)
        with Timer("Load graph"):    
            G = load_graph(ml_file)
        with Timer("ANALYSE NETWORK"):
             analyse_network_graph_tool(G, output_file_pattern, id2nodename)
    else:
        print "Using networkx"
        with Timer("LOAD NETWORK"):
            network_data = load_networkx(network_file)
        with Timer("INIT NETWORK"):
            G = init_graph_networkx(network_data, start=0, precise=5, 
                id2nodename=id2nodename)
        with Timer("ANALYSE NETWORK"):
            analyse_networkx(G, network_data, output_file_pattern, id2nodename)

def load_graph(ml_file):
    ''' Load graph for graph_tools
    '''
    print "Graph loading..."
    G = graph_tool.load_graph(ml_file, file_format="xml")
    return G

def create_graph_on_fly(network_file):
    '''
    '''
    g = Graph()
    v_trid = g.new_vertex_property("int")
    print "Load and save data..."
    id2vid = {}
    weight2edges = defaultdict(list)
    with open(network_file, "r") as fh:
        for i, line in enumerate(fh):
            if not line:
                continue
            a,b,w = line.strip().split("\t")
            a = int(a)
            b = int(b)
            w = round(float(w), 4)
            if not a in id2vid:
                v = g.add_vertex()
                v_trid[v] = a
                id2vid[a] = v
            if not b in id2vid:
                v = g.add_vertex()
                v_trid[v] = b
                id2vid[b] = v
            e = g.add_edge(id2vid[a], id2vid[b])
            weight2edges[w].append(e)
    return g, v_trid, weight2edges

def analyse_network_graph_tool(G, output_file_pattern, id2nodename):
    ''' Analys graph, created by graphtool.'''
    
    print "Get edges..."
    edge2weigth = G.edge_properties["weight"]
    weigth2edges = defaultdict(list)
    for i, e in enumerate(G.edges()):
         weigth2edges[edge2weigth[e]].append(e)
    edge2weigth = None
    weights = weigth2edges.keys()
    weights.sort()
    
    last = -1
    last_component = -1
    ended = False

    node2trf_id = G.vertex_properties["trid"]

    print "Edges:", i

    print "Iterate over weights"
    for w in weights:
        for e in weigth2edges[w]:
            G.remove_edge(e) 
        components, hist = topology.label_components(G)
        N = len(hist)
        if N == last_component:
            continue
        last_component = N

        # generate data comp_id->trf_ids sorted by size
        comp_data = [[] for i in xrange(N)]
        for node_id, comp_id in enumerate(components.a):
            comp_data[comp_id].append(int(node2trf_id[G.vertex(node_id)]))
        comp_data.sort(key=lambda x: len(x), reverse=True)

        print "Distance: %s | Edges: %s | Components: %s | Largest: %s" % (w,
                                                                         G.num_edges(),
                                                                         N,
                                                                         max(hist)
                                                                         )
        ended = True
        for hist_item in hist:
            if hist_item > 1:
                ended = False
                break
        
        output_file = output_file_pattern % (int(w), N)
        write_classification_graph_tool(output_file, comp_data, id2nodename)

        if ended:
            print "Ended since ended=True"
            break

def write_classification_graph_tool(output_file, components, id2nodename):
    ''' Save slice data to file.
    
    - components is a default dictionary: slice_id -> [trid list]
    - trid2meta is a dictionary: trid -> descripion string
    '''
    with open(output_file, "w") as fw:
        for i, comp in enumerate(components):
            for trf_id in comp:
                fw.write("%s\t%s" % (i, id2nodename[int(trf_id)]))

def write_classification_neworkx(output_file, components, trid2meta):
    ''' Save slice data to file.
    
    - components is a default dictionary: slice_id -> [trid list]
    - trid2meta is a dictionary: trid -> descripion string
    '''
    if os.path.isfile(output_file):
        os.unlink(output_file)
    with open(output_file, "w") as fw:
        for i, comp in enumerate(components):
            for trid in comp:
                fw.write("%s\t%s\n" % (i, trid2meta[int(trid)].strip()))

def create_graphml(network_file, ml_file, id2nodename):
    ''' Create graphml xml file from tab delimited network_file.
    '''

    start = """<?xml version="1.0" encoding="UTF-8"?>
<graphml xmlns="http://graphml.graphdrawing.org/xmlns"  
    xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
    xsi:schemaLocation="http://graphml.graphdrawing.org/xmlns
     http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd">
<key id="d0" for="node" attr.name="trid" attr.type="double"/>
<key id="d1" for="edge" attr.name="weight" attr.type="double"/>
<graph id="G" edgedefault="undirected">
            """
    node = '<node id="n%s"><data key="d0">%s</data></node>\n'
    edge = '<edge id="e%s" source="n%s" target="n%s"><data key="d1">%s</data></edge>\n'
    end = '</graph>\n</graphml>'

    print "Load and save data..."
    with open(ml_file, "w") as fw:
        fw.write(start)
        seen = {}
        print "Write node data"
        for k, id in enumerate(id2nodename):
            print k
            fw.write(node % (k, id))
            seen[id] = k
        print "Write edge data"
        with open(network_file, "r") as fh:
            for i, line in enumerate(fh):
                if not line:
                    continue
                a,b,w = line.strip().split("\t")
                a = int(a)
                b = int(b)
                w = round(float(w), 4)
                fw.write(edge % (i, seen[a], seen[b], w))
        fw.write(end)
    
#
# Networkx section
#

def load_networkx(network_file):
    ''' Load network data if it was pickled previously.'''
    
    print "Network data reading..."

    a_nodes = array.array("I")
    b_nodes = array.array("I")
    weights = array.array("f")

    with open(network_file, "r") as fh:
        network_data = fh.readlines()
    print "Network data parsing..."

    for i, line in enumerate(network_data):
        a, b, w = line.strip().split("\t")
        w = float(w)
        a_nodes.append(int(a))
        b_nodes.append(int(b))
        weights.append(w)

    k = len(network_data)

    network_data = None

    print "#edges: ", k
    return (a_nodes, b_nodes, weights, k)

def init_graph_networkx(network_data, start=0, precise=1, id2nodename=None):
    ''' Init graph with data.'''

    a_nodes, b_nodes, weight_vals, n = network_data

    print "Graph creation..."
    G = networkx.Graph()

    print "Populate nodes..."
    for node_id in set(a_nodes).union(set(b_nodes)):
        G.add_node(int(node_id))

    print "Add edges..."
    for i in xrange(n):
        a = a_nodes[i]
        b = b_nodes[i]
        G.add_edge(a, b)
    return G

def analyse_networkx(G, network_data, output_file_pattern, id2nodename):
    ''' Analize network and save slices.
    G is a networkx graph
    network_data - (a_nodes, b_nodes, weight_vals, n) parsed network data
    '''

    a_nodes, b_nodes, weight_vals, n = network_data

    weights = []
    for i, val in enumerate(weight_vals):
        weights.append((val, i))

    weights.sort(reverse=True)

    print "Read meta trid"
    if id2nodename:
        trid2meta = id2nodename

    print "Network analyzing..."
    last = -1
    last_component = -1
    ended = False
    last_singletons = -1

    while weights:
        (val, k) = weights.pop()
        if val == last:
            G.remove_edge(a_nodes[k], b_nodes[k])
            last = val
            continue


        components = networkx.connected_components(G)
        if len(components) == last_component:
            G.remove_edge(a_nodes[k], b_nodes[k])
            last = val
            continue

        last_component = len(components)

        
        singletons = len([x for x in components if len(x) == 1])
        # if (singletons-1) == last_singletons:
        #     print "...singleton", singletons
        #     last_singletons = singletons
        #     continue

        print "Distance: %s | Edges: %s | Nodes: %s | Sing: %s | Components: %s" % (val,
                                                                         G.number_of_edges(),
                                                                         G.number_of_nodes(),
                                                                         singletons,
                                                                         len(components))
        if len(components[0]) == 1:
            ended = True
        output_file = output_file_pattern % (int(val), len(components))

        write_classification_neworkx(output_file, components, trid2meta)

        if ended:
            break

        G.remove_edge(a_nodes[k], b_nodes[k])
        last = val
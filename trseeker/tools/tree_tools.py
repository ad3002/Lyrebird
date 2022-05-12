#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: ???
#@author: ???
#@contact: ???
# Copyright 2005-2008 by Frank Kauff & Cymon J. Cox. All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.
#
# Bug reports welcome: fkauff@biologie.uni-kl.de or on Biopython's bugzilla.
'''
Functions related to tree datastructure.
'''

import types

class ChainException(Exception):
    pass

class NodeException(Exception):
    pass

class TreeException(Exception):
    pass

class NexusParserException(Exception):
    pass

class Node(object):
    ''' A single node.
    '''

    def __init__(self, data=None):
        '''Represents a node with one predecessor and multiple successors 
        '''
        self.id = None
        if data:
            self.data = data
        else:
            self.data = NodeData()
        self.prev = None
        self.succ = []

    def set_id(self, id):
        '''Sets the id of a node, if not set yet: (self,id).'''
        if self.id is not None:
            raise NodeException('Node id cannot be changed.')
        self.id = id

    def get_id(self):
        '''Returns the node's id: (self).'''
        return self.id

    def get_succ(self):
        '''Returns a list of the node's successors: (self).'''
        return self.succ

    def get_prev(self):
        '''Returns the id of the node's predecessor: (self).'''
        return self.prev

    def add_succ(self, id):
        '''Adds a node id to the node's successors: (self,id).'''
        if isinstance(id, type([])):
            self.succ.extend(id)
        else:
            self.succ.append(id)

    def remove_succ(self, id):
        '''Removes a node id from the node's successors: (self,id).'''
        self.succ.remove(id)

    def set_succ(self, new_succ):
        '''Sets the node's successors: (self,new_succ).'''
        if not isinstance(new_succ, type([])):
            raise NodeException('Node successor must be of list type.')
        self.succ = new_succ

    def set_prev(self, id):
        '''Sets the node's predecessor: (self,id).'''
        self.prev = id

    def get_data(self):
        '''Returns a node's data: (self).'''
        return self.data

    def set_data(self, data):
        '''Sets a node's data: (self,data).'''
        self.data = data

class NodeData(object):
    '''Stores tree-relevant data associated with nodes.'''
    def __init__(self, taxon=None, branchlength=1.0,
                 support=None, comment=None, bootstrap=None,
                 taxon_name=None, distance=None, elements=None):
        self.taxon = taxon
        self.taxon_name = taxon_name
        self.branchlength = branchlength
        self.support = support
        self.comment = comment
        self.bootstrap = bootstrap
        self.distance = distance
        if elements is None:
            self.elements = set()
        else:
            self.elements = elements
        self.size = len(self.elements)

class Chain(object):
    '''Stores a list of nodes that are linked together.'''

    def __init__(self):
        '''Initiates a node chain'''
        self.chain = {}
        self.id = -1

    def _get_id(self):
        '''Gets a new id for a node in the chain.'''
        self.id += 1
        return self.id

    def all_ids(self):
        '''Return a list of all node ids.'''
        return self.chain.keys()

    def add(self, node, prev=None):
        '''Attaches node to another: (self, node, prev).'''
        if prev is not None and prev not in self.chain:
            raise ChainException('Unknown predecessor: ' + str(prev))
        else:
            id = self._get_id()
            node.set_id(id)
            node.set_prev(prev)
            if prev is not None:
                self.chain[prev].add_succ(id)
            self.chain[id] = node
        return id

    def collapse(self, id):
        '''Deletes node from chain and relinks successors to predecessor: collapse(self, id).'''
        if id not in self.chain:
            raise ChainException('Unknown ID: ' + str(id))
        prev_id = self.chain[id].get_prev()
        self.chain[prev_id].remove_succ(id)
        succ_ids = self.chain[id].get_succ()
        for i in succ_ids:
            self.chain[i].set_prev(prev_id)
        self.chain[prev_id].add_succ(succ_ids)
        node = self.chain[id]
        self.kill(id)
        return node

    def kill(self, id):
        '''Kills a node from chain without caring to what it is connected: kill(self,id).'''
        if id not in self.chain:
            raise ChainException('Unknown ID: ' + str(id))
        else:
            del self.chain[id]

    def unlink(self, id):
        '''Disconnects node from his predecessor: unlink(self,id).'''
        if id not in self.chain:
            raise ChainException('Unknown ID: ' + str(id))
        else:
            prev_id = self.chain[id].prev
            if prev_id is not None:
                self.chain[prev_id].succ.pop(self.chain[prev_id].succ.index(id))
            self.chain[id].prev = None
            return prev_id

    def link(self, parent, child):
        '''Connects son to parent: link(self,son,parent).'''
        if child not in self.chain:
            raise ChainException('Unknown ID: ' + str(child))
        elif parent not in self.chain:
            raise ChainException('Unknown ID: ' + str(parent))
        else:
            self.unlink(child)
            self.chain[parent].succ.append(child)
            self.chain[child].set_prev(parent)

    def is_parent_of(self, parent, grandchild):
        '''Check if grandchild is a subnode of parent: is_parent_of(self,parent,grandchild).'''
        if grandchild == parent or grandchild in self.chain[parent].get_succ():
            return True
        else:
            for sn in self.chain[parent].get_succ():
                if self.is_parent_of(sn, grandchild):
                    return True
            else:
                return False

    def trace(self, start, finish):
        '''Returns a list of all node_ids between two nodes (excluding start, including end): trace(start,end).'''
        if start not in self.chain or finish not in self.chain:
            raise NodeException('Unknown node.')
        if not self.is_parent_of(start, finish) or start == finish:
            return []
        for sn in self.chain[start].get_succ():
            if self.is_parent_of(sn, finish):
                return [sn] + self.trace(sn, finish)

class Tree(Chain):
    '''Represents a tree using a chain of nodes with on predecessor (=ancestor)
    and multiple successors (=subclades).
    '''

    def __init__(self, weight=1.0, rooted=False,
                 name="", data=NodeData, max_support=1.0):
        '''Ntree(self,tree).'''
        super(Tree, self).__init__()
        self.dataclass = data
        self.raw_node_data = {}
        self.max_support = max_support
        self.weight = weight
        self.rooted = rooted
        self.name = name
        root = Node(data())
        self.root = self.add(root)

    def create_from_tree(self, tree):
        raise NotImplemented

    def is_leaf(self, node):
        return type(node) is not types.ListType

    def is_subtree(self, node):
        return type(node) is types.ListType

    def _add_subtree(self, parent_id=None, tree=None):
        '''Adds leaf or tree (in newick format) to a parent_id. 
        (self,parent_id,tree).
        '''
        if parent_id is None:
            raise TreeException('Need node_id to connect to.')
        for st in tree:
            nd = self.dataclass()
            nd = self._add_nodedata(nd, st)
            if self.is_subtree(st): # it's a subtree
                sn = Node(nd)
                self.add(sn, parent_id)
                self._add_subtree(sn.id, st)
            else: # it's a leaf
                nd.taxon = st
                leaf = Node(nd)
                self.add(leaf, parent_id)

    def _add_nodedata(self, nd, st):
        if self.is_leaf(st):
            nd.taxon_name = self.raw_node_data["id_tax_dict"][st]
        key = frozenset(self.get_node_to_one_d(st))
        if key in self.raw_node_data["p_values"]:
            nd.bootstrap = self.raw_node_data["p_values"][key]
            nd.support = self.raw_node_data["p_values"][key]
        if key in self.raw_node_data["h_values"]:
            nd.distance = self.raw_node_data["h_values"][key]
        return nd

    def get_node_to_one_d(self, tree):
        '''
        Return list view of tree or node
        '''
        result_tree = []
        def tree_traversal (node, result_tree):
            # do_prework(node)
            if not self.is_leaf(node):
                for child in node:
                    tree_traversal(child, result_tree)
            else:
                # do_leafwork(node)
                result_tree.append(node)
            # do_postwork(node)
        tree_traversal(tree, result_tree)
        return result_tree

    def node(self, node_id):
        '''Return the instance of node_id.
        node = node(self,node_id)
        '''
        if node_id not in self.chain:
            raise TreeException('Unknown node_id: %d' % node_id)
        return self.chain[node_id]

    def get_terminals(self):
        '''Return a list of all terminal nodes.'''
        return [i for i in self.all_ids() if self.node(i).succ == []]

    def is_terminal(self, node):
        '''Returns True if node is a terminal node.'''
        return self.node(node).succ == []

    def is_internal(self, node):
        '''Returns True if node is an internal node.'''
        return len(self.node(node).succ) > 0

    def is_preterminal(self, node):
        '''Returns True if all successors of a node are terminal ones.'''
        if self.is_terminal(node):
            return False not in [self.is_terminal(n) for n in self.node(node).succ]
        else:
            return False

    def count_terminals(self, node=None):
        '''Counts the number of terminal nodes that are attached to a node.'''
        if node is None:
            node = self.root
        return len([n for n in self._walk(node) if self.is_terminal(n)])

    def is_identical(self, tree2):
        '''Compare tree and tree2 for identity.

        result = is_identical(self,tree2)
        '''
        return self.set_subtree(self.root) == tree2.set_subtree(tree2.root)

    def common_ancestor(self, node1, node2):
        '''Return the common ancestor that connects two nodes.
        
        node_id = common_ancestor(self,node1,node2)
        '''

        l1 = [self.root] + self.trace(self.root, node1)
        l2 = [self.root] + self.trace(self.root, node2)
        return [n for n in l1 if n in l2][-1]

    def distance(self, node1, node2):
        '''Add and return the sum of the branchlengths between two nodes.
        dist = distance(self,node1,node2)
        '''
        ca = self.common_ancestor(node1, node2)
        return self.sum_branchlength(ca, node1) + self.sum_branchlength(ca, node2)

    def get_distance_for_all_from_node(self, node1_id):
        '''
        Rerurn dictionary taxon_id : distance from given node_id
        '''
        result = {}
        for node2_id in self.get_terminals():
            d = self.distance(node1_id, node2_id)
            node = self.node(node2_id)
            result[node.data.taxon] = d
        return result

    def sum_branchlength(self, root=None, node=None):
        '''Adds up the branchlengths from root (default self.root) to node.
        
        sum = sum_branchlength(self,root=None,node=None)
        '''
        if root is None:
            root = self.root
        if node is None:
            raise TreeException('Missing node id.')
        blen = 0.0
        while node is not None and node is not root:
            blen += self.node(node).data.branchlength
            node = self.node(node).prev
        return blen

    def _walk(self, node=None):
        '''Return all node_ids downwards from a node.'''

        if node is None:
            node = self.root
        for n in self.node(node).succ:
            yield n
            for sn in self._walk(n):
                yield sn

    def split(self, parent_id=None, n=2, branchlength=1.0):
        '''Speciation: generates n (default two) descendants of a node.
        
        [new ids] = split(self,parent_id=None,n=2,branchlength=1.0):
        '''
        if parent_id is None:
            raise TreeException('Missing node_id.')
        ids = []
        parent_data = self.chain[parent_id].data
        for i in range(n):
            node = Node()
            if parent_data:
                node.data = self.dataclass()
                # each node has taxon and branchlength attribute
                if parent_data.taxon:
                    node.data.taxon = parent_data.taxon + str(i)
                node.data.branchlength = branchlength
            ids.append(self.add(node, parent_id))
        return ids

    def bootstrap(self, value):
        '''Bootstrap tree with given cutoff value
        Nodes with suppert below cutoff will be collapsed
        '''
        for id in self.all_ids():
            node = self.node(id)
            if node.data.support:
                if node.data.support <= value:
                    self.collapse(id)

    def has_support(self, node=None):
        '''Returns True if any of the nodes has data.support != None.'''
        for n in self._walk(node):
            if self.node(n).data.support:
                return True
        else:
            return False

    def prune(self, taxon):
        '''Prunes a terminal taxon from the tree.
        
        id_of_previous_node = prune(self,taxon)
        If taxon is from a bifurcation, the connectiong node will be collapsed
        and its branchlength added to remaining terminal node. This might be no
        longer a meaningful value'
        '''

        id = self.search_taxon(taxon)
        if id is None:
            raise TreeException('Taxon not found: %s' % taxon)
        elif id not in self.get_terminals():
            raise TreeException('Not a terminal taxon: %s' % taxon)
        else:
            prev = self.unlink(id)
            self.kill(id)
            if len(self.node(prev).succ) == 1:
                if prev == self.root: # we deleted one branch of a bifurcating root, then we have to move the root upwards
                    self.root = self.node(self.root).succ[0]
                    self.node(self.root).branchlength = 0.0
                    self.kill(prev)
                else:
                    succ = self.node(prev).succ[0]
                    new_bl = self.node(prev).data.branchlength + self.node(succ).data.branchlength
                    self.collapse(prev)
                    self.node(succ).data.branchlength = new_bl
            return prev

    def get_taxa(self, node_id=None):
        '''Return a list of all otus downwards from a node (self, node_id).

        nodes = get_taxa(self,node_id=None)
        '''

        if node_id is None:
            node_id = self.root
        if node_id not in self.chain:
            raise TreeException('Unknown node_id: %d.' % node_id)
        if self.chain[node_id].succ == []:
            if self.chain[node_id].data:
                return [self.chain[node_id].data.taxon]
            else:
                return None
        else:
            list = []
            for succ in self.chain[node_id].succ:
                list.extend(self.get_taxa(succ))
            return list

    def search_taxon(self, taxon):
        '''Returns the first matching taxon in self.data.taxon. Not restricted to terminal nodes.
        
        node_id = search_taxon(self,taxon)
        '''
        for id, node in self.chain.items():
            if node.data.taxon == taxon:
                return id
        return None

    def __str__(self):
        '''Short version of to_string(), gives plain tree'''
        return self.to_string(plain=True)

    def to_string(self, support_as_branchlengths=False, branchlengths_only=False, plain=True, plain_newick=False, ladderize=None):
        '''Return a paup compatible tree line.
       
        to_string(self,support_as_branchlengths=False,branchlengths_only=False,plain=True)
        '''
        # if there's a conflict in the arguments, we override plain=True
        if support_as_branchlengths or branchlengths_only:
            plain = False
        self.support_as_branchlengths = support_as_branchlengths
        self.branchlengths_only = branchlengths_only
        self.plain = plain

        def make_info_string(self, data, terminal=False):
            '''Creates nicely formatted support/branchlengths.'''
            # CHECK FORMATTING
            if self.plain: # plain tree only. That's easy.
                return ''
            elif self.support_as_branchlengths: # support as branchlengths (eg. PAUP), ignore actual branchlengths
                if terminal:    # terminal branches have 100% support
                    return ':%1.2f' % self.max_support
                else:
                    return ':%1.2f' % (data.support)
            elif self.branchlengths_only: # write only branchlengths, ignore support
                return ':%1.5f' % (data.branchlength)
            else:   # write suport and branchlengths (e.g. .con tree of mrbayes)
                if terminal:
                    return ':%1.5f' % (data.branchlength)
                else:
                    if data.branchlength is not None and data.support is not None:  # we have blen and suppport
                        return '%1.2f:%1.5f' % (data.support, data.branchlength)
                    elif data.branchlength is not None:                             # we have only blen
                        return '0.00000:%1.5f' % (data.branchlength)
                    elif data.support is not None:                                  # we have only support
                        return '%1.2f:0.00000' % (data.support)
                    else:
                        return '0.00:0.00000'

        def ladderize_nodes(self, nodes, ladderize=None):
            '''Sorts node numbers according to the number of terminal nodes.'''
            if ladderize in ['left', 'LEFT', 'right', 'RIGHT']:
                succnode_terminals = [(self.count_terminals(node=n), n) for n in nodes]
                succnode_terminals.sort()
                if (ladderize == 'right' or ladderize == 'RIGHT'):
                    succnode_terminals.reverse()
                if succnode_terminals:
                    succnodes = zip(*succnode_terminals)[1]
                else:
                    succnodes = []
            else:
                succnodes = nodes
            return succnodes

        def newickize(self, node, ladderize=None):
            '''Convert a node tree to a newick tree recursively.'''

            if not self.node(node).succ:    #terminal
                return "%s%s" % (self.node(node).data.taxon,
                                 make_info_string(self.node(node).data,
                                 terminal=True))
            else:
                succnodes = ladderize_nodes(self.node(node).succ, ladderize=ladderize)
                subtrees = [newickize(sn, ladderize=ladderize) for sn in succnodes]
                return '(%s)%s' % (','.join(subtrees), make_info_string(self.node(node).data))

            #treeline = ['tree']
            treeline = []
            # if self.name:
            #     treeline.append(self.name)
            # else:
            #     treeline.append('a_tree')
            # treeline.append('=')
            if self.weight != 1:
                treeline.append('[&W%s]' % str(round(float(self.weight), 3)))
            if self.rooted:
                treeline.append('[&R]')
            succnodes = ladderize_nodes(self.node(self.root).succ)
            subtrees = [newickize(sn, ladderize=ladderize) for sn in succnodes]
            treeline.append('(%s)' % ','.join(subtrees))
            if plain_newick:
                return treeline[-1]
            else:
                return ' '.join(treeline) + ';'


def get_slice_files(folder):
    ''' Get list of (slice_id, file_path).'''
    slice_files = []
    for file_path in sc_iter_filepath_folder(folder):
            slice_id = int(file_path.split("_")[-1].replace(".net", ""))
            data = [slice_id, file_path]
            slice_files.append(data)
    slice_files.sort()
    return slice_files

def read_tree_slices(slice_files):
    ''' Read all slice files.'''
    for i, (id, file_path) in enumerate(slice_files):
        print i, "\r",
        result = []
        for slice_obj in sc_iter_tab_file(file_path, NetworkSliceModel):
            result.append((slice_obj.gid, slice_obj.trf_id))
        slice_files[i] = result
    print "Done"
    return slice_files
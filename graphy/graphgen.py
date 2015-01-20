from __future__ import division, print_function, absolute_import
import six
range = six.moves.range

import networkx as nx
import numpy as np

def hierarchical(n, level):
    """ Generates hiearchical graph using method proposed of:
    Ravasz E, Barabasi AL, Hierarchical organization in complex networks, PRE, 2003.

    Parameters
    ----------
    n : int
        Number of nodes
    level : int
        Number of hierarchical levels to create

    """

    if level == 0:
        return nx.complete_graph(n)
    else:
        fullG = nx.Graph()

        #get lower-order graphs
        for i in range(n):
            fullG = nx.union(hierarchical(n, level-1), fullG, rename=(str(i)+'.',''))
        
        edges = []
        suffix = ''
        for l in range(level-1):
            suffix += '.0'
            
        #connect outer nodes to the center
        center = '0.0'+suffix
        for node in fullG.nodes():
            if not '0' in node:
                edges.append((node, center))
        
        fullG.add_edges_from(edges)
        return fullG

def get_hierarchical_pos(net):
    """ Get x,y positions for plotting hierarchical graph.

    For example:

    >>> from graphy import graphgen
    >>> G = graphgen.hierarchical(3, 2)
    >>> pos = graphgen.get_hierarchical_pos(G)
    ...

    Parameters
    ----------
    net : networkx graph
        Graph whose nodes to layout

    """

    pos = {}

    for node in net.nodes():
        x, y = 0, 0
        xl = [0, -1, -1,  1,  1]
        yl = [0, -1,  1, -1,  1]
        for level_ndx, level in enumerate(map(int, node.title().split("."))):
            x += 0.5**level_ndx*(xl[level])
            y += 0.5**level_ndx*(yl[level])

        pos[node] = (x,y)

    return pos

def get_block_matrix(membership, intra_community_p, inter_community_p):
    membership = np.asarray(membership)
    N  = len(membership)
    mx = np.random.rand(N,N) > (1-inter_community_p)
    for comm in set(membership):
        ixs = np.flatnonzero(membership == comm)
        comm_size = len(ixs)
        for r in ixs:
            mx[r,ixs] = np.random.rand(comm_size) > (1-intra_community_p)
    np.fill_diagonal(mx, 0)
    return mx

def get_barbell_matrix(membership, num_conns=1):
    membership = np.asarray(membership)
    comms = list(set(membership))
    N  = len(membership)
    mx = np.zeros(shape=(N,N))
    for ndx1, comm in enumerate(comms):
        ixs = np.flatnonzero(membership == comm)
        for r in ixs:
            mx[r,ixs] = 1
        for ndx2, othercomm in enumerate(comms):
            if ndx1 == ndx2:
                continue
            ixs2 = np.flatnonzero(membership == othercomm)
            for cndx in range(num_conns):
                mx[ixs[cndx], ixs2[cndx]] = 1
                mx[ixs2[cndx], ixs[cndx]] = 1
    np.fill_diagonal(mx, 0)
    return mx


"""Module implementing tools to generate graphs and connectivity matrices.
Generally, `gen_...` functions return generatively-created objects, while
`get_...` return graphs corresponding to passed in membership vectors. 

"""


from __future__ import division, print_function, absolute_import
import six
range = six.moves.range

import networkx as nx
import numpy as np
import itertools

def gen_hierarchical_net(n, level):
    """Generate hiearchical graph using method proposed of:
    Ravasz E, Barabasi AL, Hierarchical organization in complex networks, PRE, 2003.

    Parameters
    ----------
    n : int
        Number of nodes
    level : int
        Number of hierarchical levels to create

    Returns
    -------
    networkx graph
        The hierchically-structured network

    """

    if level == 0:
        return nx.complete_graph(n)
    else:
        fullG = nx.Graph()

        #get lower-order graphs
        for i in range(n):
            fullG = nx.union(gen_hierarchical_net(n, level-1), fullG, rename=(str(i)+'.',''))
        
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

def get_hierarchical_net_pos(net):
    """ Get x,y positions for plotting hierarchical graph.

    For example:

    .. plot::
        :include-source:

        >>> from graphy import graphgen
        >>> G = graphgen.gen_hierarchical_net(5, 2)
        >>> pos = graphgen.get_hierarchical_net_pos(G)
        >>> import networkx as nx
        >>> nx.draw_networkx(G, with_labels=False, node_size=50, pos=pos)

    Parameters
    ----------
    net : networkx graph
        Graph whose nodes to layout

    Returns
    -------
    dict 
        Dictionary containing node:(x,y) entries

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


def sample_connection_matrix(prob_mx):
    """Sample from a matrix of connection probabilities in order to create a 
    binary connection matrix with 0s on the diagonal.

    Parameters
    ----------
    prob_mx : 2-dimensional np.array of floats
        Matrix of connection probabilities

    Returns
    -------
    2-dimensional np.array
        Binary connectivity matrix

    """

    mx = (np.random.rand(*prob_mx.shape) > (1-prob_mx)).astype('int')
    np.fill_diagonal(mx, 0)
    return mx


def get_weighted_block_matrix(membership, intra_community_w, inter_community_w):
    """Get weighted block-structured matrix corresponding to membership vector
    with different weights for intra- versus inter-community connections, and 
    0s on the diagonals.

    For example:

    .. plot::
        :include-source:

        >>> from graphy import graphgen
        >>> cmx = graphgen.get_weighted_block_matrix([0,0,0,0,0,1,1,1,1,1], 0.5, 0.1)
        >>> plt.imshow(cmx, interpolation='none') # doctest: +SKIP

    Parameters
    ----------
    membership : list or np.array of ints
        Array containing assignment of each node to communities
    intra_community_w : float
        Weight for intra-community connections
    inter_community_w : float
        Weight for inter-community connections

    Returns
    -------
    np.array matrix
        The connectivity matrix
        
    """

    membership = np.asarray(membership)
    N  = len(membership)
    mx = np.zeros(shape=(N,N)) + inter_community_w
    for comm in set(membership):
        ixs = np.flatnonzero(membership == comm)
        for r in ixs:
            mx[r,ixs] = intra_community_w
    np.fill_diagonal(mx, 0)
    return mx

def gen_hierarchical_weighted_block_matrix(blocksize, numblocks, numlevels, level_weights):
    """Generate hierarchical weighted block matrix. 

    Parameters
    ----------
    blocksize : int
        Number of nodes to include in each lowest-level block
    numblocks : int
        Number of blocks each level consists of
    numlevels : int
        Number of levels
    level_weights : list of float
        Strength of connection between members for each level (plus one more for 
        the 'top' level). Order is from lowest-level (smallest blocks) to highest-level

    Returns
    -------
    np.array matrix
        The generated connectivity matrix


    For example:

    .. plot::
        :include-source:

        >>> from graphy import graphgen
        >>> cmx = graphgen.gen_hierarchical_weighted_block_matrix(4, 4, 2, [0.3, 0.2, 0.1])
        >>> plt.imshow(cmx, interpolation='none') # doctest: +SKIP

    """

    def get_match_level(i,j):
        K = len(i)
        for ndx, (ival, jval) in enumerate(zip(i,j)):
            if ival != jval:
                return K-ndx
        return 0

    N = blocksize * (numblocks ** numlevels)
    mx = np.zeros( (N,N) )
    blockcoords = list(itertools.product(*[list(range(numblocks)) for i in range(numlevels)]))

    for indx, i in enumerate(blockcoords):
        for jndx, j in enumerate(blockcoords):
            w = level_weights[get_match_level(i,j)]
            mx[indx*blocksize:(indx+1)*blocksize,jndx*blocksize:(jndx+1)*blocksize] = w

    np.fill_diagonal(mx, 0)
    return mx


def get_barbell_matrix(membership, num_conns=1):
    """Get a matrix of completely-connected communities
    connected by paths corresponding to membership vector.

    For example:

    .. plot::
        :include-source:

        >>> from graphy import graphgen
        >>> cmx = graphgen.get_barbell_matrix([0,0,0,0,0,1,1,1,1,1])
        >>> import matplotlib.pylab as plt
        >>> plt.imshow(cmx, interpolation='none') # doctest: +SKIP


    Parameters
    ----------
    membership : list or np.array of ints
        Array containing assignment of each node to communities
    num_conns : int
        Number of links that should run between communities 

    Returns
    -------
    np.array matrix
        The connectivity matrix
        
    """

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


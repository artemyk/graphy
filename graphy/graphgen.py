"""Module implementing tools to generate graphs and connectivity matrices.
Generally, ``gen_...`` functions return generatively-created objects, while
``get_...`` return graphs corresponding to passed in membership vectors. 

Use ``sample_connection_matrix`` to generate binary connectivity matrices
sampled from matrices of connection probabilities.

"""


from __future__ import division, print_function, absolute_import
import six
range = six.moves.range

import networkx as nx
import numpy as np
import itertools
import scipy

def gen_ring_matrix(N, neighs_per_side=1):
    """Generate a ring-lattice matrix.

    For example:

    .. plot::
      :include-source:

      >>> import graphy
      >>> mx = graphy.graphgen.gen_ring_matrix(30, 5)
      >>> graphy.plotting.plot_graph(mx) # doctest: +SKIP


    Parameters
    ----------
    N : int
        How many nodes.
    neighs_per_side : int
        How many neighbors on each side should be connected to each node.

    Returns
    -------
    np.array matrix
        The connectivity matrix
        
    """

    if (2*neighs_per_side + 1) > N:
        raise ValueError('For a ring of size %d, neighs_per_side can ' 
            'have maximum value of %d' % (N, int((N-1)/2)))

    v = np.zeros(N)
    v[1:(neighs_per_side+1)] = 1
    v[-(neighs_per_side):] = 1
    return scipy.linalg.toeplitz(v)


def get_clique_of_rings_net_and_pos(sizes, neighs_per_side=1, ring_weight=1.0, clique_weight=1.0):
    """Generate several ring-lattice graphs of different sizes.  Then 
    choose a single node from each ring and interconnect those into 
    a clique.

    .. plot::
      :include-source:

      >>> import graphy
      >>> net, pos = graphy.graphgen.get_clique_of_rings_net_and_pos([5, 10, 20])
      >>> graphy.plotting.plot_graph(net, pos=pos) # doctest: +SKIP    


    Parameters
    ----------
    sizes : list of int
        How many nodes in each ring.
    neighs_per_side : int (default 1)
        For nodes in ring lattices, how many neighbors on each side to connect to.
    ring_weight : float (default 1.0)
        Strength of connections in each ring lattice.
    clique_weight : float (default 1.0)
        Strength of connections in clique connecting single nodes in each ring.

    Returns
    -------
    networkx graph
        The clique of rings network
    dict
        { node : xyposition } dictionary of positions for good node layout

    """

    mx = np.zeros((sum(sizes),sum(sizes)))
    offset=0
    N = sum(sizes)
    pos = {}
    interconnect_nodes = []

    groundtruth = np.hstack([np.ones(s,dtype='int')*ndx for ndx, s in enumerate(sizes)])
    for ndx, s in enumerate(sizes):
        centerrad = 2 * np.pi * np.sum(np.sqrt(sizes[:ndx])) / np.sum(np.sqrt(sizes))
        base_pos = np.sqrt(N*1)*np.array([np.cos(centerrad),np.sin(centerrad)])
        rads = 2 * np.pi * np.linspace(0, 1, s, endpoint=False)
        radius = np.sqrt(s)
        nodepos = base_pos + radius * np.vstack([np.cos(rads), np.sin(rads)]).T

        interconnect_nodes.append(offset+np.argsort(np.linalg.norm(nodepos, axis=1))[0])
        
        for i in range(s):
            pos[offset+i] = nodepos[i]

        mx[offset:offset+s,offset:offset+s] = ring_weight * gen_ring_matrix(s, neighs_per_side)
        offset+=s

    for i in interconnect_nodes:
        mx[i,interconnect_nodes] = clique_weight

    np.fill_diagonal(mx, 0)

    return nx.from_numpy_matrix(mx), pos


def gen_hierarchical_net(n, level):
    """Generate hiearchical graph using method proposed of:
    Ravasz E, Barabasi AL, Hierarchical organization in complex networks, PRE, 2003.

    Parameters
    ----------
    n : int
        Number of nodes in the lowest level.
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

        >>> import graphy
        >>> G = graphy.graphgen.gen_hierarchical_net(5, 2)
        >>> pos = graphy.graphgen.get_hierarchical_net_pos(G)
        >>> graphy.plotting.plot_graph(G, pos=pos)

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


    For example:

    .. plot::
        :include-source:

        >>> from graphy import graphgen
        >>> cmx = graphgen.gen_hierarchical_weighted_block_matrix(4, 4, 2, [0.3, 0.2, 0.1])
        >>> plt.imshow(cmx, interpolation='none') # doctest: +SKIP


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
            mlevel = get_match_level(i,j)
            if mlevel >= len(level_weights):
                raise ValueError('Not enough level_weights specified')
            w = level_weights[mlevel]
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


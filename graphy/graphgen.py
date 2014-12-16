import networkx as nx
import numpy as np

def hierarchical(n, level):
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

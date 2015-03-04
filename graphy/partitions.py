"""Module implementing tools to work with partitions of sets as well as 
perform heuristic optimization over set partitions.
"""

from __future__ import division, print_function, absolute_import
import six
range = six.moves.range

import numpy as np
import networkx as nx
import matplotlib.pylab as plt
import copy

def to_str(membership):
    """Convert membership array to pretty string.

    Example:

    >>> from graphy import partitions
    >>> print(partitions.to_str([0,0,0,1,1,1]))
    [0 0 0 1 1 1]

    Parameters
    ----------
    membership : np.array or list
        Membership array to convert

    Returns
    -------
    str
        Pretty string

    """

    return "[" + " ".join(map(str, membership)) + "]"

def to_alphanum_str(membership):
    """Convert membership array to short string.

    Example:

    >>> from graphy import partitions
    >>> print(partitions.to_alphanum_str([0,0,0,1,1,1,2,3,4,20,20,20]))
    [000111234KKK]

    Parameters
    ----------
    membership : np.array or list
        Membership array to convert

    Returns
    -------
    str
        Short pretty string

    """
    membership = np.asarray(membership)
    names = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz@#$%^&*_=+-/?><"
    if membership.max() >= len(names):
        return to_str(membership)
    return "[" + "".join([names[i] for i in membership]) + "]"


def get_minsize_assignment(N, min_comm_size):
    """Create membership vector where each community contains at least
    as a certain number of nodes.

    Parameters
    ----------
    N : int
        Desired length of membership vector
    min_comm_size : int
        Minimum number of nodes each community should have.

    Returns
    -------
    np.array
        Membership vector

    """
    num_comms = int(N / min_comm_size)
    membership = -np.ones(N, dtype='int')  # -1 means non-assigned
    for c in range(num_comms):
        left_to_assign = np.flatnonzero(membership == -1)
        assign = np.random.choice(left_to_assign, min_comm_size, replace=False)
        membership[assign] = c

    membership[membership == -1] = np.random.randint(num_comms, size=np.sum(membership == -1))
    
    return membership




def renumber_membership(membership):
    """Renumber membership vector so that community numbers begin with 0 and increase
    consecutively.

    Example:

    >>> from graphy import partitions
    >>> print(partitions.renumber_membership([10,10,10,31,31,31,32,33,4,20]))
    [3 3 3 5 5 5 0 1 2 4]

    Parameters
    ----------
    membership : np.array
        Membership array to renumber

    Returns
    -------
    np.array
        Renumber membershp array

    """
    membership = np.asarray(membership)
    r = membership.copy()
    remap = { old_comm:new_comm for new_comm, old_comm in enumerate(set(membership)) }
    for i in range(len(membership)):
        r[i] = remap[membership[i]]
    return r



def remap2match(partition1, partition2):
    """Renumber membership assignments so that community identities have maximal 
    overlap with another membership assignment.

    For example:

    >>> from graphy import partitions
    >>> print(remap2match([3,3,1,1,0],[2,2,3,3,3]))
    [2 2 3 3 4]

    Parameters
    ----------
    partition1 : np.array or list
        Membership assignment to remap
    partition2 : np.array or list
        Membership assignment to match

    Returns
    -------
    np.array
        Remapped assignment

    """
    partition1 = np.asarray(partition1)
    partition2 = np.asarray(partition2)

    nmap = {}
    to_remap = set(partition1)
    gtlist = partition2.tolist() if isinstance(partition2, np.ndarray) else partition2
    allowed_matches = set(gtlist + list(range(partition2.max()+1,partition2.max()+len(partition1))))
    while len(to_remap):
        max_overlap, saved_pair = None, None
        for c1 in to_remap:
            for c2 in allowed_matches:
                overlap = np.logical_and(partition1 == c1, partition2 == c2).sum()
                if max_overlap is None or overlap > max_overlap:
                    max_overlap = overlap
                    saved_pair = (c1, c2)
        old_c, new_c = saved_pair
        nmap[old_c] = new_c
        to_remap        = to_remap        - set([old_c,])
        allowed_matches = allowed_matches - set([new_c,])
    return np.array([nmap[c] for c in partition1], dtype='int')




def find_optimal_across_time(qualityObj, timepoints, num_runs=1):
    import igraph

    saved_best = []
    best_membership = None
    last_time = 0
    last_best_membership = None
    last_best_membership_q = 0
    for t in sorted(timepoints):
        qualityObj.end_state_advance(t - last_time)
        last_time = t

        #remap_target = groundtruth if groundtruth is not None else best_membership

        # Run search
        best_membership, best_membership_q = qualityObj.find_optimal(debug_level=0, 
            initial_membership=last_best_membership,
            num_runs=num_runs)
        #best_membership2 = graphy.partitions.find_optimal(qualityObj, debug_level=0, num_runs=num_runs) # random start
        #best_membership = best_membership1 if qualityObj.quality(best_membership1) >= qualityObj.quality(best_membership2) else best_membership2

        if last_best_membership is not None:
            if np.abs(last_best_membership_q - best_membership_q) < 1e-4:
                best_membership = last_best_membership
            else:
                best_membership = renumber_membership(best_membership)
            vi  = igraph.compare_communities(best_membership, last_best_membership, method='vi')
            nmi = igraph.compare_communities(best_membership, last_best_membership, method='nmi')
        else:
            vi = 0.0
            nmi = 0.0

        #if groundtruth is not None:
        #    #Relabel to best match
        #    best_membership = graphy.partitions.remap2match(best_membership, groundtruth)

        print('t=%2d vi=%0.4f nmi=%0.4f #=%2d q=%0.4f %s' % (t, vi, nmi, len(set(best_membership)), best_membership_q, to_alphanum_str(best_membership)))
        #print 'Ground truth   : q=%0.4f %s' % (qualityObj.quality(groundtruth)    , graphy.partitions.to_str(groundtruth))

        saved_best.append( (t, best_membership, copy.deepcopy(qualityObj)) )
        last_best_membership, last_best_membership_q = best_membership, best_membership_q

    return saved_best


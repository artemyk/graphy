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
import random


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
        if max_overlap == 0:
            new_c = max(list(nmap.values()) + [0,]) + 1
        nmap[old_c] = new_c
        to_remap        = to_remap        - set([old_c,])
        allowed_matches = allowed_matches - set([new_c,])
    return np.array([nmap[c] for c in partition1], dtype='int')


def greedy_search(qualityfunc, N, initial_membership=None, num_runs=1, debug_level=0):
    """Find optimal decomposition using a greedy search.

    Parameters
    ----------
    qualityfunc : function
        Function which returns the quality measure to be maximized.
    N : int
        Number of nodes to partition.
    initial_membership : np.array, optional
        Initial membership assignment.  If None specified, each component is
        assigned to separate subsystem.
    num_runs : int, optional
        Number of runs to try, can improve quality of decompositions. Default
        is 1.
    debug_level : int, optional
        Amount of debugging information to display, from 0 (no debugging 
        information) to 3 (maximal debugging information)

    Returns
    -------
    np.array
        Optimal membership array
    float
        Q value corresponding to optimal membership

    """

    class CommMerger(object):
        @staticmethod
        def get_elements(membership):
            return list(set(membership))

        @staticmethod
        def prop_memberships(el, membership):
            for diff_comm in list(set(membership)):
                if el == diff_comm:
                    continue

                prop_membership = membership.copy()
                prop_membership[prop_membership == el] = diff_comm 

                yield prop_membership

    class CommSpliter(CommMerger):
        @staticmethod
        def prop_memberships(el, membership):
            c_nodes = np.flatnonzero(membership == el)
            if len(c_nodes) <= 1:
                return

            about_half = (len(c_nodes)+1)/2
            new_comm = max(membership)+1
            for _ in range(10):
                random.shuffle(c_nodes)
                prop_membership = membership.copy()
                prop_membership[c_nodes[:about_half]] = new_comm
                yield prop_membership

    class NodeMover(object):
        @staticmethod
        def get_elements(membership):
            return list(range(len(membership)))

        @staticmethod
        def prop_memberships(el, membership):
            for diff_comm in list(set(membership)):
                if membership[el] == diff_comm:
                    continue

                prop_membership = membership.copy()
                prop_membership[el] = diff_comm 

                yield prop_membership

    class NodeSwapper(NodeMover):
        @staticmethod
        def prop_memberships(el, membership):
            for diff_el_ndx in range(len(membership)):
                if membership[el] == membership[diff_el_ndx]:
                    continue

                prop_membership = membership.copy()
                prop_membership[el], prop_membership[diff_el_ndx] = prop_membership[diff_el_ndx], prop_membership[el] 

                yield prop_membership

    _done = set()
    def get_quality(membership):
        _done.add(tuple(membership.tolist()))
        return qualityfunc(membership)

    def greedy_moves(membership, mover_class):

        if debug_level >= 1:
            classname = mover_class.__name__.ljust(15)

        old_quality = None
        cur_quality = get_quality(membership)

        iter_num = 0   
        while old_quality is None or cur_quality > (old_quality + 1e-5):
            old_quality = cur_quality
            elements = mover_class.get_elements(membership)
            random.shuffle(elements)

            for v in elements:

                all_proposed = [m for m in mover_class.prop_memberships(v, membership) if tuple(m.tolist()) not in _done]

                if not len(all_proposed):
                    continue

                random.shuffle(all_proposed)

                best_move_quality, best_move_membership = cur_quality, None
                for c in all_proposed:
                    q = get_quality(c)
                    if debug_level >= 4:
                        print(classname, 
                              "Trying: %s -> %s [q=%0.3f vs. old q=%0.3f]"
                              % (to_str(membership), to_str(c), q, cur_quality)
                             )
                    if q >= best_move_quality:
                        best_move_quality = q
                        best_move_membership = c

                if best_move_quality > cur_quality: 
                    cur_quality = best_move_quality
                    if debug_level >= 3:
                        print(classname, 
                              "Accepted move: %s -> %s [q=%0.3f]"
                              % (to_str(membership), to_str(best_move_membership), best_move_quality)
                             )

                    membership = best_move_membership

            membership = renumber_membership(membership)

            if debug_level >= 2:
                print(classname, 
                      "Iteration %d, #=%d quality=%5.3f (improvement=%5.3f), m=%s" %
                      (iter_num, len(set(membership)), cur_quality, cur_quality - old_quality, to_str(membership))
                     )
                
        return membership, cur_quality

    # ***************************************************
    # Main function body
    # ***************************************************

    best_membership, best_quality = None, None
    for i in range(num_runs):

        if initial_membership is None:
            membership = np.arange(N, dtype='int')
        else:
            if len(initial_membership) != N:
                raise ValueError(
                  'Length of initial_membership (%d) is different from expected (%d)' % 
                  (len(initial_membership), N) )
            membership = initial_membership.copy()

        if debug_level >= 1:
            print("*** Run %d ***" % i)

        old_quality, cur_quality = None, None
        while old_quality is None or cur_quality >= (old_quality + 1e-5):
            old_quality = cur_quality
            membership, cur_quality = greedy_moves(membership, mover_class=NodeMover)
            #membership, cur_quality = greedy_moves(membership, mover_class=NodeSwapper)
            membership, cur_quality = greedy_moves(membership, mover_class=CommMerger)
            membership, cur_quality = greedy_moves(membership, mover_class=NodeMover)
            membership, cur_quality = greedy_moves(membership, mover_class=CommSpliter)
            
        if best_quality is None or best_quality < cur_quality:
            best_membership = membership
            best_quality = cur_quality
            
    return best_membership, best_quality



def find_optimal_across_time(qualityObj, timepoints, num_runs=1, debug_level=0):
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

        vi = 0.0
        nmi = 0.0

        if last_best_membership is not None:
            if np.abs(last_best_membership_q - best_membership_q) < 1e-4:
                best_membership = last_best_membership
            else:
                best_membership = renumber_membership(best_membership)
            nmi = igraph.compare_communities(best_membership, last_best_membership, method='nmi')

        if debug_level > 0:
            print('t=%2d nmi=%0.4f #=%2d q=%0.4f %s' % (t, nmi, len(set(best_membership)), best_membership_q, to_alphanum_str(best_membership)))

        saved_best.append( (t, best_membership, copy.deepcopy(qualityObj)) )
        last_best_membership, last_best_membership_q = best_membership, best_membership_q

    return saved_best


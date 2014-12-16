import numpy as np
import random
import networkx as nx
import matplotlib.pylab as plt

def to_str(membership):
    return "[" + " ".join(map(str, membership)) + "]"

class FindOptimal(object):

    def __init__(self, N, quality_func_obj, debug_level=0):

        self.N = N
        self.quality_func_obj = quality_func_obj
        self.debug_level = debug_level

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
            return range(len(membership))

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

    def greedy_moves(self, membership, mover_class):
        old_quality = None
        cur_quality = self.quality_func_obj.quality(membership)

        iter_num = 0   
        while old_quality is None or cur_quality > (old_quality + 1e-5):
            old_quality = cur_quality
            elements = mover_class.get_elements(membership)
            random.shuffle(elements)

            for v in elements:

                all_proposed = list(mover_class.prop_memberships(v, membership))

                if not len(all_proposed):
                    continue

                random.shuffle(all_proposed)

                memb_qualities = []
                for c in all_proposed:
                    q = self.quality_func_obj.quality(c)
                    memb_qualities.append((c, q))

                best_move_membership, best_move_quality = sorted(memb_qualities, reverse=True, key=lambda x: x[1])[0] 

                if best_move_quality > cur_quality: 
                    cur_quality = best_move_quality
                    if self.debug_level >= 2:
                        print mover_class.__name__, "Accepted move: %s -> %s [q=%0.3f]" \
                            % (to_str(membership), to_str(best_move_membership), best_move_quality)

                    membership = best_move_membership

            remap = { old_comm:new_comm for new_comm, old_comm in enumerate(set(membership)) }
            for i in range(len(membership)):
                membership[i] = remap[membership[i]]

            if self.debug_level >= 1:
                print "Iteration %d, #=%d quality=%5.3f (improvement=%5.3f), m=%s, cls=%s" % \
                        (iter_num, len(set(membership)), cur_quality, cur_quality - old_quality, to_str(membership), mover_class.__name__)
                
        return membership, cur_quality
    
    def run(self, initial_membership = None):
        if initial_membership is None:
            membership = np.arange(self.N, dtype='int32')
        else:
            membership = initial_membership.copy()

        for i in range(2):
            if self.debug_level >= 1:
                print "*** Run through %d ***" % i

            old_quality, cur_quality = None, None

            while old_quality is None or cur_quality >= (old_quality + 1e-5):
                old_quality = cur_quality
                membership, cur_quality = self.greedy_moves(membership, mover_class=self.NodeMover)
                membership, cur_quality = self.greedy_moves(membership, mover_class=self.NodeSwapper)
                membership, cur_quality = self.greedy_moves(membership, mover_class=self.CommMerger)
                membership, cur_quality = self.greedy_moves(membership, mover_class=self.CommSpliter)

        return membership


def get_minsize_assignment(N, min_comm_size):
    # Min-size community assignment
    num_comms = int(N / min_comm_size)
    membership = -np.ones(N, dtype='int32')  # -1 means non-assigned
    for c in range(num_comms):
        left_to_assign = np.flatnonzero(membership == -1)
        assign = np.random.choice(left_to_assign, min_comm_size, replace=False)
        membership[assign] = c

    membership[membership == -1] = np.random.randint(num_comms, size=np.sum(membership == -1))
    
    return membership


def remap2match(partition1, partition2):
    nmap = {}
    to_remap = set(partition1)
    gtlist = partition2.tolist() if isinstance(partition2, np.ndarray) else partition2
    allowed_matches = set(gtlist + range(partition2.max()+1,partition2.max()+len(partition1)))
    while len(to_remap):
        max_overlap, saved_pair = None, None
        for c1 in to_remap:
            for c2 in allowed_matches:
                overlap = np.logical_and(partition1 == c1, partition2 == c2).sum()
                if overlap > max_overlap:
                    max_overlap = overlap
                    saved_pair = (c1, c2)
        old_c, new_c = saved_pair
        nmap[old_c] = new_c
        to_remap        = to_remap        - set([old_c,])
        allowed_matches = allowed_matches - set([new_c,])
    return np.array([nmap[c] for c in partition1], dtype='int32')


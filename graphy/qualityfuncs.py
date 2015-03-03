"""Module implements quality functions for graph decompositions.
"""

from __future__ import division, print_function, absolute_import
import six
range = six.moves.range

import numpy as np
import scipy

import random

from . import partitions


class QualityFunction(object):
  # TODO: Document that should implement N attribute
  def __init__(self):
    """Implements a quality function.
    """
    # precompute quantities here
    pass
  
  def quality(self, membership):
    """Returns a quality score corresponding to membership vector.

    Parameters
    ----------
    membership : np.array
        Membership vector

    Returns
    -------
    float
        Quality

    """
    raise NotImplementedError

  def find_optimal(self, initial_membership=None, num_runs=1, debug_level=0):
    """Find optimal decomposition.

    Parameters
    ----------
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
        return self.quality(membership)

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

                #memb_qualities = []
                best_move_quality, best_move_membership = cur_quality, None
                for c in all_proposed:
                    q = get_quality(c)
                    if debug_level >= 4:
                        print(classname, 
                              "Trying: %s -> %s [q=%0.3f vs. old q=%0.3f]"
                              % (partitions.to_str(membership), partitions.to_str(c), q, cur_quality)
                             )
                    #memb_qualities.append((c, q))
                    if q >= best_move_quality:
                        best_move_quality = q
                        best_move_membership = c

                #best_move_membership, best_move_quality = sorted(memb_qualities, reverse=True, key=lambda x: x[1])[0] 

                if best_move_quality > cur_quality: 
                    cur_quality = best_move_quality
                    if debug_level >= 3:
                        print(classname, 
                              "Accepted move: %s -> %s [q=%0.3f]"
                              % (partitions.to_str(membership), partitions.to_str(best_move_membership), best_move_quality)
                             )

                    membership = best_move_membership

            membership = partitions.renumber_membership(membership)

            if debug_level >= 2:
                print(classname, 
                      "Iteration %d, #=%d quality=%5.3f (improvement=%5.3f), m=%s" %
                      (iter_num, len(set(membership)), cur_quality, cur_quality - old_quality, partitions.to_str(membership))
                     )
                
        return membership, cur_quality
    
    # ***************************************************
    # Main function body
    # ***************************************************

    best_membership, best_quality = None, None
    for i in range(num_runs):

        if initial_membership is None:
            membership = np.arange(self.N, dtype='int')
        else:
            if len(initial_membership) != self.N:
                raise ValueError(
                  'Length of initial_membership (%d) is different from expected (%d)' % 
                  (len(initial_membership), self.N) )
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
    
class Modularity(QualityFunction):
  def __init__(self, mx):
    """Class that implements Newman's modularity quality function.

    Parameters 
    ----------
    mx : 2-dimensional np.array
      Connectivity matrix of graph to partition into communities.

    """
    self.mx = np.asarray(mx).astype('float')
    self.ks = self.mx.sum(axis=1)
    self.total_stubs = self.ks.sum()
    self.N = mx.shape[0]

  def quality(self, membership):
    q = 0
    for comm in set(membership):
        ixs = membership == comm
        q += self.mx[ixs,:][:,ixs].sum() - np.outer(self.ks[ixs], self.ks[ixs]/self.total_stubs).sum()
    return q

class DirectedModularity(Modularity):
  def __init__(self, mx):
    """Class that implements Newman's directed modularity quality function.

    EA Leicht, MEJ Newman, "Community structure in directed networks", PRL, 2007.
    http://arxiv.org/abs/0709.4500

    Parameters 
    ----------
    mx : 2-dimensional np.array
      Connectivity matrix of graph to partition into communities.

    """
    super(DirectedModularity, self).__init__(mx)
    self.k_in  = self.mx.sum(axis=0)
    self.k_out = self.mx.sum(axis=1)

  def quality(self, membership):
    q = 0.0
    for comm in set(membership):
        ndxs = membership == comm
        q += (self.mx[ndxs,:][:,ndxs].sum() - np.outer(self.k_in[ndxs], self.k_out[ndxs]).sum() / self.total_stubs) / self.total_stubs
    return q


class InfoMapCodeLength(QualityFunction):
  def __init__(self, mx, teleportation_prob=0.0):
    """Class that implements Rosvall's `map equation' quality function.


    Parameters 
    ----------
    mx : 2-dimensional np.array
      Connectivity matrix of graph to partition into communities.
    teleportation_prob : float (default 0.0)
      Probability of teleporting to random node.

    """
    mx = np.asarray(mx)
    self.N = mx.shape[0]
    mx = mx.astype('float')
    
    # Create transition matrix
    self.trans = mx/mx.sum(axis=0)[None,:]
    
    if teleportation_prob > 0:
      teleportation_trans = np.ones(mx.shape)
      np.fill_diagonal(teleportation_trans, 0)
      teleportation_trans /= teleportation_trans.sum(axis=0)[None,:]
      
      self.trans = (1-teleportation_prob)*self.trans + teleportation_prob*teleportation_trans
    
    # Get equilibrium probs and entropy
    evals, evecs = scipy.linalg.eig(self.trans)
    equil = np.real(evecs[:,np.isclose(evals,1)])
    self.equil = np.ravel(equil / equil.sum())
    self.equil_entropy = self.normed_entropy(self.equil)
        
  @staticmethod
  def normed_entropy(probs):
    #probs = probs[probs != 0.0]
    psum = sum(probs)
    probs = probs / (psum if psum != 0 else 1)
    return -np.sum(probs * np.log2(probs))

  def quality(self, membership, extrainfo=False):
    comms = set(membership)
    exit_probs = []
    total_code_len = 0.0
    for c in comms:
      ixs             = membership == c
      comm_node_probs = self.equil[ixs]
      comm_p          = comm_node_probs.sum()
      stay_prob       = (self.trans[:,ixs][ixs,:].dot(comm_node_probs)).sum()
      comm_exit_prob  = (comm_p - stay_prob)

      if comm_exit_prob > 0:
        plist = np.array(comm_node_probs.tolist() + [comm_exit_prob,])
      else:
        plist = comm_node_probs

      comm_code_len   = (comm_p + comm_exit_prob) * self.normed_entropy(plist)
      total_code_len += comm_code_len
      exit_probs.append(comm_exit_prob)
        
    total_code_len += sum(exit_probs) * self.normed_entropy(np.array(exit_probs))

    # Return negative because we want to minimize, not maximize
    return -total_code_len

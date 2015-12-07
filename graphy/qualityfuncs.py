"""Module implements quality functions for graph decompositions.
"""

from __future__ import division, print_function, absolute_import
import six
range = six.moves.range

import numpy as np
import scipy
import scipy.sparse as sp

from . import partitions


class QualityFunction(object):
  def __init__(self):
    """Base class for implementing a quality function.
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
    return partitions.greedy_search(self.quality, self.N, 
      initial_membership=initial_membership,
      num_runs=num_runs,
      debug_level=debug_level)

    
class Modularity(QualityFunction):
  def __init__(self, mx):
    """Class that implements Newman's modularity quality function.

    Parameters 
    ----------
    mx : 2-dimensional np.array
      Connectivity matrix of graph to partition into communities.

    """
    mx = mx if sp.isspmatrix(mx) else np.asarray(mx)
    self.mx = mx.astype('float')
    self.ks = self.mx.sum(axis=1)
    self.total_stubs = self.ks.sum()
    self.N = mx.shape[0]

  def quality(self, membership):
    if len(membership) != self.N:
        raise ValueError('Membership should be array of length %d, not %d' % 
                         (self.N, len(membership)))

    q = 0
    for comm in set(membership):
        ixs = membership == comm
        q += (self.mx[ixs,:][:,ixs].sum() - self.ks[ixs].sum()**2/self.total_stubs)
    return q / self.total_stubs

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
    self.k_in  = self.mx.sum(axis=0).flat
    self.k_out = self.mx.sum(axis=1).flat

  def quality(self, membership):
    q = 0.0
    for comm in set(membership):
        ndxs = membership == comm
        q += (self.mx[ndxs,:][:,ndxs].sum() - (self.k_in[ndxs].sum()*self.k_out[ndxs].sum()) / self.total_stubs)
    return q  / self.total_stubs


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
      # np.fill_diagonal(teleportation_trans, 0)  
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
    K = len(comms)
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
    
    if K > 1:
      total_code_len += sum(exit_probs) * self.normed_entropy(np.array(exit_probs))

    # Return negative because we want to minimize, not maximize
    return -total_code_len

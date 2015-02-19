"""Module implements quality functions for graph decompositions.
"""

from __future__ import division, print_function, absolute_import
import six
range = six.moves.range

import numpy as np


class QualityFunction(object):
  def __init__(self):
    """Implements a quality function.
    """
    # precompute quantities here
    pass

  @property
  def N(self):
    """Returns the number of nodes (length of membership vector)
    considered by this QualityFunction.
    """ 
    raise NotImplementedError

  
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
    
class Modularity(QualityFunction):
  def __init__(self, mx):
    """Class that implements Newman's modularity quality function.
    """
    self.mx = np.asarray(mx).astype('float')
    self.ks = self.mx.sum(axis=1)
    self.total_stubs = self.ks.sum()
    self._N = mx.shape[0]

  @property
  def N(self):
    return self._N

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

"""Module implementing useful functions for working with graphs.
"""


from __future__ import division, print_function, absolute_import
import six
range = six.moves.range

import numpy as np


def get_laplacian(mx):
  """Return laplacian of a graph represented by a matrix.

  Parameters
  ----------
  mx : 2-dimensional np.array
      Matrix representing graph

  Returns
  -------
  2-dimensional np.array
      Laplacian of graph

  """

  degs = np.sum(mx,axis=1)
  return -(np.diag(degs) - mx)



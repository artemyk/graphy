from __future__ import division, print_function, absolute_import
import six
range = six.moves.range

import numpy as np
import networkx as nx
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm

def plot_membership(N, membership, pos=None, ax=None, colormap_name='Paired'):
  """Plot a circle visualizing community membership of nodes. 

  For example:

  .. plot::
      :include-source:

      >>> from graphy import plotting
      >>> import numpy as np
      >>> plotting.plot_membership(12, np.arange(12)/3)
      ...


  Parameters
  ----------
  N : int
      Number of nodes.
  membership : list
      Community membership vector.
  pos : list of tuples
      List of (x,y) positions of nodes.  If not specified, nodes
      are arranged along a circle.
  ax : matplotlib axis object (default is current axis)
      Matplotlib axis to plot onto.
  colormap_name : string (default 'Paired')
      Name of colormap to use to visualize community assignments.

  """

  class MplColorHelper:
    def __init__(self, cmap_name, start_val, stop_val):
      self.cmap_name = cmap_name
      self.cmap = plt.get_cmap(cmap_name)
      self.norm = mpl.colors.Normalize(vmin=start_val, vmax=stop_val)
      self.scalarMap = cm.ScalarMappable(norm=self.norm, cmap=self.cmap)

    def get_rgb(self, val):
      return self.scalarMap.to_rgba(val)

  membership = np.asarray(membership)

  if pos is None:
      pos = np.asarray([[np.cos(angle), np.sin(angle)] 
                        for angle in np.arange(0,2*np.pi, 2*np.pi/N)])
      
  g = nx.Graph()
  for i in range(N):
      g.add_node(i)

  if pos is None:
      pos = nx.spring_layout(g, weight='weight')

  COL = MplColorHelper(colormap_name, 0, N)

  membership_colors = []
  for m in membership:
    new_col = np.sum([v*2**-(p+1) 
                      for p, v in enumerate(np.array(map(int, bin(m)[2:]))[::-1])])
    membership_colors.append(int(N*new_col))
  
  nx.draw(g, 
          pos=pos, 
          font_size=10,
          width=2,
          node_size=200,
          with_labels=False, 
          node_color=COL.get_rgb(membership_colors), 
          ax=ax)

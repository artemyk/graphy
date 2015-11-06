from __future__ import division, print_function, absolute_import
import six
range = six.moves.range

import collections

import numpy as np
import networkx as nx
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib import cm
from matplotlib.patches import Ellipse
from scipy.spatial.distance import cdist


def plot_graph(G, pos=None, colors=None, node_labels=None, node_size=0.04, 
  edgescale=1.0, nodeopts={}, labelopts={}, arrowopts={}, cmap='Paired'):
  """Plot a graphs.  Supports both directed and undirected graphs.
  Undirected edges are drawn as lines while directed edges are drawn
  as arrows.

  For example:

  .. plot::
      :include-source:

      >>> from graphy import plotting
      >>> import networkx as nx
      >>> G=nx.karate_club_graph()
      >>> plotting.plot_graph(G, pos=nx.spring_layout(G), colors=range(G.number_of_nodes())) # doctest: +SKIP


  Parameters
  ----------
  G : networkx Graph object or 2-d np.array
      Graph to plot, either instance of networkx Graph or a 2-d connectivity 
      matrix.
  pos : dict
      Dict specifying positions of nodes, as in {node: (x,y).  If not provided, 
      nodes are arranged along a circle.
  colors : list of ints (default None)
      Color(s) to use for node faces, if desired.
  node_labels : list of strings (default None)
      Labels to use for node labels, if desired.
  node_size : float (default 0.05)
      Size of nodes.
  edgescale : float (default 1.0)
      Controls thickness of edges between nodes.
  nodeopts : dict (default {})
      Extra options to pass into plt.Circle call for plotting nodes.
  labelopts : dict (default {})
      Extra options to pass into plt.text call for plotting labels. 
      Could be used to specify fontsize, for example.
  arrowopts : dict (default {})
      Extra options to pass into plt.arrow call for plotting edges.
  cmap : string (default 'Paired')
      Name of colormap to use.
  """


  class MplColorHelper:
    """Class that helps pick colors from specified colormaps.
    """
    def __init__(self, cmap_name, start_val, stop_val):
      self.cmap_name = cmap_name
      self.cmap = plt.get_cmap(cmap_name)
      self.norm = mpl.colors.Normalize(vmin=start_val, vmax=stop_val)
      self.scalarMap = cm.ScalarMappable(norm=self.norm, cmap=self.cmap)

    def get_rgb(self, val):
      return self.scalarMap.to_rgba(val)

  def intersect_two_circles(xy1, r1, xy2, r2):
    d = np.linalg.norm(xy2-xy1)

    a = (r1**2 - r2**2 + d**2) / (2 * d)
    b = d - a
    h = np.sqrt(r1**2 - a**2)
    xy3 = xy1 + a * (xy2 - xy1) / d

    ix = xy3[0] + h * (xy2[1] - xy1[1]) / d
    iy = xy3[1] - h * (xy2[0] - xy1[0]) / d

    return np.array([ix, iy])

  if isinstance(G, (np.ndarray, np.generic) ):
    G = nx.from_numpy_matrix(G, create_using=nx.DiGraph())
  elif not isinstance(G, nx.Graph):
    raise ValueError('Unknown type of graph: %s' % str(type(G)))

  if pos is None:
    pos = nx.circular_layout(G)

  if colors is None:
    colors = 1

  if not isinstance(colors, collections.Iterable):
    colors = [colors,] * G.number_of_nodes()

  colors = np.asarray(colors)

  cmap_helper = MplColorHelper(cmap, colors.min(), colors.max())
  
  bbox = plt.gca().get_window_extent() 
  asp_ratio = bbox.width/bbox.height

  node_map = { n:ndx for ndx, n in enumerate(G.nodes())}
  xys = np.array([pos[n] for n in G.nodes()])
  xys -= xys.min(axis=0)
  maxvalues = xys.max(axis=0)
  maxvalues[maxvalues == 0] = 1
  xys /= maxvalues
  xys[:,0] *= asp_ratio
  
  plt.xlim([-(3*node_size)*asp_ratio,(1+(3*node_size))*asp_ratio])
  plt.ylim([-(3*node_size),1+(3*node_size)])
  plt.axis('off')
  
  try:
    edge_weights = nx.get_edge_attributes(G, 'weight')
  except KeyError:
    edge_weights = {}
      
  headscale = node_size*0.5 if G.is_directed() else 0
  arrowdict = dict(ec='k', fc='k', 
      length_includes_head=True, 
      shape='full',
      head_length=headscale,
      head_width =headscale)
  for k, v in arrowopts.items():
      arrowdict[k] = v

  done_edges = []
  for edge in G.edges():
    startxy = xys[node_map[edge[0]],:].copy()
    endxy   = xys[node_map[edge[1]],:].copy()

    edgesize = edge_weights.get(edge,1.0)*edgescale
    arrowdict['lw'] = edgesize

    if edge[0] == edge[1]:
        loopoffset = np.sign(startxy - xys.mean(axis=0)) * node_size * 1.05
        cloop = plt.Circle(startxy+loopoffset, radius=node_size*0.65, 
                           ec='k', fill=False, lw=edgesize)
        plt.gca().add_artist(cloop)
        arrowloc = intersect_two_circles(startxy, node_size, 
                                         startxy+loopoffset, node_size*0.7)
        arrowlocstart = arrowloc + (arrowloc - startxy)*1e-5
        plt.arrow(arrowlocstart[0], arrowlocstart[1],  
                  arrowloc[0]-arrowlocstart[0], arrowloc[1]-arrowlocstart[1], 
                  **arrowdict)

    else:
        angle   = np.arctan2(endxy[1]-startxy[1], endxy[0]-startxy[0])
        offset  = np.array([np.cos(angle),np.sin(angle)])*node_size
        startxy += offset
        endxy   -= offset            
        
        if not nx.is_directed(G) or (edge[1], edge[0]) in done_edges:
            if (edge[1], edge[0]) in done_edges: # will be drawn in the other direction
                continue
            else:
                midxy = (startxy + endxy) / 2.0
                plt.arrow(midxy[0], midxy[1], 
                          endxy[0]-midxy[0], endxy[1]-midxy[1], **arrowdict)
                plt.arrow(midxy[0], midxy[1], 
                          startxy[0]-midxy[0], startxy[1]-midxy[1], **arrowdict)
        else:
            plt.arrow(startxy[0], startxy[1], 
                      endxy[0]-startxy[0], endxy[1]-startxy[1], **arrowdict)

    done_edges.append((edge[0], edge[1]))
        

  for ndx, xy in enumerate(xys):  # Plot nodes
      cnode = plt.Circle((xy[0],xy[1]), radius=node_size, ec='none',
                         color=cmap_helper.get_rgb(colors[ndx]), **nodeopts)
      plt.gca().add_artist(cnode)
      if node_labels is not None:
          plt.text(xy[0],xy[1], node_labels[ndx], ha='center', va='center', 
                   **labelopts)
        

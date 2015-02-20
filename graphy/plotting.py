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
      >>> plotting.plot_graph(G, pos=nx.spring_layout(G), colors=range(G.number_of_nodes()))
      ...


  Parameters
  ----------
  G : networkx Graph object
      Graph to plot.
  pos : list of tuples
      List of (x,y) positions of nodes.  If not specified, nodes
      are arranged along a circle.
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

  xys = np.array([pos[ndx] for ndx in range(len(pos))])
  xys -= xys.min(axis=0)
  xys /= xys.max(axis=0)
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

  for edge in G.edges():
    startxy = xys[edge[0],:].copy()
    endxy   = xys[edge[1],:].copy()

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
        
        if not nx.is_directed(G) or G.has_edge(edge[1], edge[0]):
            if edge[0] > edge[1]: # will be drawn in the other direction
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
        

  for ndx, xy in enumerate(xys):  # Plot nodes
      cnode = plt.Circle((xy[0],xy[1]), radius=node_size, ec='none',
                         color=cmap_helper.get_rgb(colors[ndx]), **nodeopts)
      plt.gca().add_artist(cnode)
      if node_labels is not None:
          plt.text(xy[0],xy[1], node_labels[ndx], ha='center', va='center', 
                   **labelopts)
        


def DrawModularityFigure(mod_ts, optmod_ts=None, data_ts=None, time=None, 
                         change_points=None, vis_change_points=None,
                         node_size=0.05,
                         y1label='Perturbation \n Modularity', 
                         y2label='Change in Phase', x1label='Time',
                         state_linewidth=2,
                         changepoint_linewidth=1,
                         network_axis_size=0.2
                         ):
  """Make the Perturbation Modularity plots. 

  Parameters
  ----------
  mod_ts : array, num_time x num_partitions
      Modularity over time for partitions of interest (blue dashed).
  optmod_ts : array, num_time x 1
      Modularity over time for optimal partiton (red solid).
  data_ts : array, num_time x SystemSize
      Data timeseries for second plot (phase or change in phase, ect...)
  time : array
      Time points for plots
  change_points : dict ({ time : membership_vector })
      Dictionary of partitions keyed by the first time they become optimal
  vis_change_points : list of time (default None)
      Which change_points to visualize using network plots.  By default, all
      change points will be visualized.
  node_size : int (default 0.05)
      node size for partition graphs
  state_linewidth : int (default 2)
      Linewidth to use for the state plots.
  changepoint_linewidth : int (default 1)
      Linewidth to use for the changepoint lines.
  network_axis_size : float (default 0.2)
      Size of the axis object for visualizing network plots.

  """
    
  black_color = '#262626'
  blue_color = '#779ECB'
  red_color = '#C23B22'
  
  if time is None:
      time = range(len(mod_ts))
  
  fig, (ax1, ax2) = plt.subplots(2, sharex=True, sharey=False, figsize=(6,6))
  #plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.1, hspace=0.1)
  fig.subplots_adjust(top=0.75, hspace=0.0, right=0.95, left=0.15)
  
  # play with the axis
  spine_list = ['top', 'right', 'left', 'bottom']
  for spine in spine_list:
      ax1.spines[spine].set_color(black_color)
      ax2.spines[spine].set_color(black_color)

  ax1.xaxis.label.set_color(black_color)
  ax1.yaxis.label.set_color(black_color)

  ax2.xaxis.label.set_color(black_color)
  ax2.yaxis.label.set_color(black_color)

  for xlabel in ax2.axes.get_xticklabels():
      xlabel.set_fontsize(12)

  for ylabel in ax1.axes.get_yticklabels():
      ylabel.set_fontsize(12)

  for ylabel in ax2.axes.get_yticklabels():
      ylabel.set_fontsize(12)
  
  #make last y2 label hidden
  ax2.yaxis.set_major_locator(MaxNLocator(nbins=len(ax2.get_yticklabels()), prune='upper')) # added 
  
  ax1.set_ylabel(y1label, fontsize=14)
  ax2.set_xlabel(x1label, fontsize=14)
  ax2.set_ylabel(y2label, fontsize=14)
  
  # plot the data timeseries in the bottom axis
  color_idx = np.linspace(0, 1, data_ts.shape[1])
  for i, y in zip(color_idx, data_ts.T):
    plt.plot(time, y, color=plt.cm.Paired(i), lw=state_linewidth)
  plt.xlim([time.min(), time.max()])

  # plot the background modularity timeseries
  print(mod_ts.shape)
  print(time.shape)
  ax1.plot(time, mod_ts, color=blue_color, ls = '--', lw = 1)

  transFigure = fig.transFigure.inverted()
  
  cplist = sorted(change_points.keys())
  # get the mean position for a partition interval
  mean_cp = np.convolve(cplist + [time[-1]], 0.5*np.ones(2), 'same')
 
  if vis_change_points is None:
    vis_change_points = change_points.keys()

  vis_change_points = sorted(vis_change_points)
  num_vis_change_points = float(len(vis_change_points))

  for cp, mp in zip(cplist, mean_cp[1::]):
    coord3 = transFigure.transform(ax1.transData.transform([cp, ax1.get_ylim()[1]]))
    coord4 = transFigure.transform(ax2.transData.transform([cp, ax2.get_ylim()[0]]))
    fig.lines.append(mpl.lines.Line2D((coord3[0],coord4[0]),(coord3[1],coord4[1]),color='gray',
                               transform=fig.transFigure, lw=changepoint_linewidth))

  blankGraph = nx.Graph()
  for i in range(len(change_points.values()[0])):
    blankGraph.add_node(i)

  for cp, mp in zip(cplist, mean_cp[1::]):
    if cp not in vis_change_points:
      continue

    # this is an inset axes over the main axes
    a = plt.axes([0.1 + 0.85*vis_change_points.index(cp)/num_vis_change_points, 
                  .8, network_axis_size, network_axis_size])
    plot_graph(blankGraph, colors=change_points[cp], node_size=0.05)
    a.text(0.5, 0.5, str(len(set(change_points[cp]))), fontsize=20, ha='center', va='center', color = black_color)
    plt.setp(a, xticks=[], yticks=[])

    coord_graph = transFigure.transform(a.transData.transform([0.5, 0.5]))
    coord_fig   = transFigure.transform(ax1.transData.transform([mp, ax1.get_ylim()[1]]))
    
    llen = np.linalg.norm(coord_graph - coord_fig)
    lang = np.arctan2(coord_graph - coord_fig)
    coord_graph = (llen - 1) * np.array([np.cos(lang), np.sin(lang)])
    fig.lines.append(mpl.lines.Line2D(*zip(coord_graph, coord_fig),color=black_color,
                               transform=fig.transFigure, lw=2))

  ax1.plot(time, optmod_ts, color=red_color, lw=3)
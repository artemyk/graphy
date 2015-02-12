from __future__ import division, print_function, absolute_import
import six
range = six.moves.range

import numpy as np
import networkx as nx
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib import cm
from matplotlib.patches import Ellipse
from scipy.spatial.distance import cdist


class MplColorHelper:
  def __init__(self, cmap_name, start_val, stop_val):
    self.cmap_name = cmap_name
    self.cmap = plt.get_cmap(cmap_name)
    self.norm = mpl.colors.Normalize(vmin=start_val, vmax=stop_val)
    self.scalarMap = cm.ScalarMappable(norm=self.norm, cmap=self.cmap)

  def get_rgb(self, val):
    return self.scalarMap.to_rgba(val)


def plot_membership(membership, pos=None, ax=None, colormap_name='Paired',
  node_size=200, labels=None, font_size=10):
  """Plot a circle visualizing community membership of nodes. 

  For example:

  .. plot::
      :include-source:

      >>> from graphy import plotting
      >>> import numpy as np
      >>> plotting.plot_membership(np.arange(12)/3)
      ...


  Parameters
  ----------
  membership : list of int
      Community membership vector.
  pos : list of tuples
      List of (x,y) positions of nodes.  If not specified, nodes
      are arranged along a circle.
  ax : matplotlib axis object (default is current axis)
      Matplotlib axis to plot onto.
  colormap_name : string (default 'Paired')
      Name of colormap to use to visualize community assignments.
  node_size : int (default 200)
      Size of nodes to draw.
  labels : list of str (default None)
      Node labels, if desired.
  font_size : int (default 10)
      Font size with which to plot the labels, if provided.

  """

  membership = np.asarray(membership).astype('int')
  N = len(membership)

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
    bin_repr = np.array(list(map(int, bin(m)[2:])))
    new_col = np.sum([v*2**-(p+1) for p, v in enumerate(bin_repr[::-1])])
    membership_colors.append(int(N*new_col))
  
  nx.draw(g, 
          pos=pos, 
          width=2,
          node_size=node_size,
          with_labels=False, 
          node_color=COL.get_rgb(membership_colors), 
          ax=ax)
  
  if labels is not None:
    nx.draw_networkx_labels(g, pos, dict(zip(range(N), labels)), 
      font_size=font_size, ax=ax)


def DrawModularityFigure(mod_ts, optmod_ts=None, data_ts=None, time=None, 
                         change_points=None, vis_change_points=None,
                         node_size=100,
                         y1label='Perturbation \n Modularity', 
                         y2label='Change in Phase', x1label='Time',
                         state_linewidth=2,
                         changepoint_linewidth=2,
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
  node_size : int (default 100)
      node size for partition graphs
  state_linewidth : int (default 2)
      Linewidth to use for the state plots.
  changepoint_linewidth : int (default 2)
      Linewidth to use for the changepoint lines.

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
  ax1.plot(time, mod_ts , color = blue_color, ls = '--', lw = 1)

  transFigure = fig.transFigure.inverted()
  
  cplist = sorted(change_points.keys())
  # get the mean position for a partition interval
  mean_cp = np.convolve(cplist + [time[-1]], 0.5*np.ones(2), 'same')
 
  if vis_change_points is None:
    vis_change_points = change_points.keys()

  vis_change_points = sorted(vis_change_points)
  num_vis_change_points = float(len(vis_change_points))

  for cp, mp in zip(cplist, mean_cp[1::]):
    coord3 = transFigure.transform(ax1.transData.transform([cp, 1.0]))
    coord4 = transFigure.transform(ax2.transData.transform([cp, ax2.get_ylim()[0]]))
    fig.lines.append(mpl.lines.Line2D((coord3[0],coord4[0]),(coord3[1],coord4[1]),color=black_color,
                               transform=fig.transFigure, lw=changepoint_linewidth))

  for cp, mp in zip(cplist, mean_cp[1::]):
    if cp not in vis_change_points:
      continue

    # this is an inset axes over the main axes
    a = plt.axes([0.1 + 0.85 * vis_change_points.index(cp) / num_vis_change_points, .8, 0.2, 0.2])
    plot_membership(change_points[cp], ax = a, colormap_name='Paired', node_size=node_size)
    a.text(0.0, 0.0, str(len(set(change_points[cp]))), fontsize = 20, ha='center', va='center', color = black_color)
    plt.setp(a, xticks=[], yticks=[])

    coord1 = transFigure.transform(a.transData.transform([0, -1.3]))
    coord2 = transFigure.transform(ax1.transData.transform([mp, 1.0]))
    
    fig.lines.append(mpl.lines.Line2D((coord1[0],coord2[0]),(coord1[1],coord2[1]),color=black_color,
                               transform=fig.transFigure, lw=2))

  ax1.plot(time, optmod_ts, color = red_color, lw = 3)


def plot_graph(G, pos, membership, nodelabels=None, nodesize=0.05, edgescale=1.0, nodeopts={}, labelopts={}, arrowopts={}, cmap='Paired'):

  def intersect_two_circles(xy1, r1, xy2, r2):
    d = np.linalg.norm(xy2-xy1)

    a = (r1**2 - r2**2 + d**2) / (2 * d)
    b = d - a
    h = np.sqrt(r1**2 - a**2)
    xy3 = xy1 + a * (xy2 - xy1) / d

    ix = xy3[0] + h * (xy2[1] - xy1[1]) / d
    iy = xy3[1] - h * (xy2[0] - xy1[0]) / d

    return np.array([ix, iy])

  cmap_helper = MplColorHelper(cmap, membership.min(), membership.max())
  
  bbox = plt.gca().get_window_extent() 
  asp_ratio = bbox.width/bbox.height

  xys = np.array([pos[ndx] for ndx in range(len(pos))])
  xys -= xys.min(axis=0)
  xys /= xys.max(axis=0)
  xys[:,0] *= asp_ratio
  
  plt.xlim([-(3*nodesize)*asp_ratio,(1+(3*nodesize))*asp_ratio])
  plt.ylim([-(3*nodesize),1+(3*nodesize)])
  plt.axis('off')
  
  try:
    edge_weights = nx.get_edge_attributes(G, 'weight')
  except KeyError:
    edge_weights = {edge:1.0 for edge in graph.edges()}
      
  for edge in G.edges():
      startxy = xys[edge[0],:].copy()
      endxy   = xys[edge[1],:].copy()
      arrowdict = dict(ec='k', fc='k', 
          lw=edge_weights[edge]*edgescale, 
          length_includes_head=True, 
          shape='full',
          head_length=nodesize*0.5,
          head_width=nodesize*0.5)
      for k, v in arrowopts.iteritems():
          fdict[k] = v

      if edge[0] == edge[1]:
          loopoffset = np.sign(startxy - xys.mean(axis=0)) * nodesize * 1.05
          cloop = plt.Circle(startxy + loopoffset, radius=nodesize*0.65, ec='k', fill=False, lw=edge_weights[edge]*edgescale)
          plt.gca().add_artist(cloop)
          arrowloc = intersect_two_circles(startxy, nodesize, startxy+loopoffset, nodesize*0.7)
          arrowlocstart = arrowloc + (arrowloc - startxy)*1e-5
          plt.arrow( arrowlocstart[0], arrowlocstart[1],  arrowloc[0]-arrowlocstart[0], arrowloc[1]-arrowlocstart[1], **arrowdict)

      else:
          angle   = np.arctan2(endxy[1]-startxy[1], endxy[0]-startxy[0])
          offset  = np.array([np.cos(angle),np.sin(angle)])*nodesize
          startxy += offset
          endxy   -= offset            
          
          if False and (not nx.is_directed(G) or G.has_edge(edge[1], edge[0])):
              if edge[0] > edge[1]: # will be drawn in the other direction
                  continue
              else:
                  midxy = (startxy + endxy) / 2.0
                  plt.arrow(midxy[0], midxy[1], endxy[0]-midxy[0], endxy[1]-midxy[1], **arrowdict)
                  plt.arrow(midxy[0], midxy[1], startxy[0]-midxy[0], startxy[1]-midxy[1], **arrowdict)
          else:
              plt.arrow(startxy[0], startxy[1], endxy[0]-startxy[0], endxy[1]-startxy[1], **arrowdict)
          
          #frac = 0.05 / np.linalg.norm(np.array([ex-sx,ey-sy]))
          #plt.annotate('', xytext=(sx,sy), xy=(ex+1e-4, ey), width=2, arrowprops=dict(arrowstyle=arrowstyle, frac=frac, facecolor='black'))
              
  
  for ndx, xy in enumerate(xys):  # Plot nodes
      cnode = plt.Circle((xy[0],xy[1]), ec='none',radius=nodesize, color=cmap_helper.get_rgb(membership[ndx]), **nodeopts)
      plt.gca().add_artist(cnode)
      if nodelabels is not None:
          plt.text(xy[0],xy[1], nodelabels[ndx], ha='center', va='center', **labelopts)
        

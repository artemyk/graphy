from __future__ import division, print_function, absolute_import
import six
range = six.moves.range

import numpy as np
import networkx as nx
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib import cm


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


def DrawModularityFigure(mod_ts, optmod_ts = None, data_ts = None, time = None, change_points = None, node_size = 100,
                        y1label = 'Perturbation \n Modularity', y2label = 'Change in Phase', x1label= 'Time', filename = None):
    """Make the Perturbation Modularity plots. 

  Parameters
  ----------
  mod_ts : array
      Modularity over time for partitions of interest (blue dashed).
      num_time x num_partitions
  optmod_ts : array
      Modularity over time for optimal partiton (red solid).
      num_time x 1
  data_ts : array
      Data timeseries for second plot (phase or change in phase, ect...)
      num_time x SystemSize
  time : array
      Time points for plots
  change_points : dict
      Dictionary of partitions keyed by the first time they become optimal
  node_size : int (default 100)
      node size for partition graphs
  filename : str (default None)
     If not None, saves the figure with this filename.

  """
    
    black_color = '#262626'
    blue_color = '#779ECB'
    red_color = '#C23B22'
    
    if time is None:
        time = range(len(mod_ts))
    
    fig, (ax1, ax2) = plt.subplots(2, sharex=True, sharey=False, figsize = (6,6))
    #plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.1, hspace=0.1)
    fig.subplots_adjust(top = 0.75, hspace=0.0, right = 0.95, left = 0.15)
    
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
    ax2.plot(time, data_ts, lw = 2)
    
    # plot the background modularity timeseries
    ax1.plot(time, mod_ts , color = blue_color, ls = '--', lw = 1)

    transFigure = fig.transFigure.inverted()
    
    num_change_points = float(len(change_points))
    cplist = sorted(change_points.keys())
    # get the mean position for a partition interval
    mean_cp = np.convolve(cplist + [time[-1]], 0.5*np.ones(2), 'same')
   
    for cp, mp in zip(cplist, mean_cp[1::]):
        # this is an inset axes over the main axes
        a = plt.axes([0.1 + 0.85 * cplist.index(cp) / num_change_points, .8, 1.0/num_change_points, 1.0/num_change_points])
        plot_membership(change_points[cp], ax = a, colormap_name='Paired', node_size=node_size)
        a.text(0.0, 0.0, str(len(set(change_points[cp]))), fontsize = 20, ha='center', va='center', color = black_color)
        plt.setp(a, xticks=[], yticks=[])

        coord1 = transFigure.transform(a.transData.transform([0, -1.3]))
        coord2 = transFigure.transform(ax1.transData.transform([mp, 1.0]))
        
        coord3 = transFigure.transform(ax1.transData.transform([cp, 1.0]))
        coord4 = transFigure.transform(ax2.transData.transform([cp, 0.0]))
    
        fig.lines.append(mpl.lines.Line2D((coord1[0],coord2[0]),(coord1[1],coord2[1]),color = black_color,
                                   transform=fig.transFigure, lw = 2))
        fig.lines.append(mpl.lines.Line2D((coord3[0],coord4[0]),(coord3[1],coord4[1]),color = black_color,
                                   transform=fig.transFigure, lw = 2))

    ax1.plot(time, optmod_ts, color = red_color, lw = 3)
    
    if not filename is None:
        plt.savefig(filename)


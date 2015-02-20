"""Module for running Louvain community detection algorithm.
"""


from __future__ import division, print_function, absolute_import
import six
range = six.moves.range

import numpy as np
import subprocess
import tempfile
import os

def optimize_modularity(conn_mx, debug=False):
  """Optimize directed, weighted Newman's modularity 
  using the Louvain algorithm.
  
  This uses C++ implementation from 
  https://github.com/raldecoa/SurpriseMe/tree/master/src/CPM

  For example:

  >>> import networkx as nx
  >>> import graphy
  >>> G = nx.karate_club_graph()
  >>> best_membership, q = graphy.louvain.optimize_modularity(nx.to_numpy_matrix(G))
  >>> print(best_membership, q)   # doctest: +SKIP


  Parameters
  ----------
  conn_mx : 2-dimensional np.array
    Connectivity matrix.
  debug : bool (default False)
    If True, prints various debugging information.

  Returns
  -------
  np.array
    Optimal membership vector indicating community membership of each node.
  float
    Modularity value corresponding to the optimal membership vector. Notice that
    because the modularity value is computed by adding up increments over 
    many moves, this may only be accurate to a few decimal places.

  """

  conn_mx = np.asarray(conn_mx)

  if conn_mx.shape[0] != conn_mx.shape[1]:
      raise ValueError('conn mx should be square')

  DIR = tempfile.gettempdir()
  NETWORK_FILE = os.path.join(DIR, 'net.txt')
  OUTPUT_FILE  = os.path.join(DIR, 'net.bin')
  NODEMAP_FILE = os.path.join(DIR, 'nmap.bin')
  CONF_FILE    = os.path.join(DIR, 'conf.bin')

  with open(NETWORK_FILE, 'w') as f:
      f.write('>\n')
      for ndx in range(conn_mx.shape[0]):
          f.write('%d %d\n' % (ndx, ndx))
      f.write('>\n0 0\n>\n')
      for ndx, r in enumerate(conn_mx):
          conns = np.flatnonzero(r)
          if not len(conns):
              f.write('%d %d %0.5f 0\n' % (ndx, 0, 0))
          for c in conns:
              f.write('%d %d %0.5f 0\n' % (ndx, c, r[c]))

  if debug:
    print("**** NETWORK FILE: ****")
    with open(NETWORK_FILE) as f:
      print(f.read())
    print()

  bin_dir = os.path.join(os.path.dirname(__file__),'..','external','SurpriseMeCPM','bin')
  subprocess.call([os.path.join(bin_dir,'slicer'), 
                   '-i', NETWORK_FILE, 
                   '-o', OUTPUT_FILE, 
                   '-n', NODEMAP_FILE, 
                   '-c', CONF_FILE], stderr=subprocess.PIPE)
  res = subprocess.Popen([os.path.join(bin_dir,'community'), 
                          OUTPUT_FILE, 
                          CONF_FILE],
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE)

  output, errors = map(lambda s: s.decode('ascii'), res.communicate())

  all_files = [NETWORK_FILE, OUTPUT_FILE, NODEMAP_FILE, CONF_FILE]
  if debug:
    print("**** OUTPUT: ****")
    print(output)
    print()
    print("**** STDERR: ****")
    print(errors)
    print("*****************")
    print("Leaving files in place:")
    print(all_files)

  else:
    for fname in all_files:
      os.remove(fname)

  q_value = float(errors.strip().split("\n")[-1].split(" ")[0])
  membership = np.zeros(conn_mx.shape[0], dtype='int')
  ndxs, vals = zip(*[map(int, l.strip().split("\t")) for l in output.strip().split("\n")])
  membership[list(ndxs)] = vals

  return membership, q_value


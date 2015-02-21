"""Module for running Louvain community detection algorithm.
"""


from __future__ import division, print_function, absolute_import
import six
range = six.moves.range

import numpy as np
import subprocess
import tempfile
import os
import scipy.sparse as sp

import scipy.sparse as sp

def optimize_modularity(conn_mx, rand_init=True, num_runs=1, debug=False):
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
  conn_mx : 2-dimensional np.array or scipy.sparse matrix
    Connectivity matrix.
  rand_init : bool (default True)
    Whether to randomly shuffle order of nodes (makes results non-deterministic)
  num_runs : int (default 1)
    How many runs to perform (highest quality run returned).  Only allow if 
    rand_init is True.
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
  # TODO: Implement multithreading?  Code seems to support it already

  if num_runs > 1 and not rand_init:
    raise ValueError('Multiple runs only makes sense when initial order of'
                     'nodes is randomized')

  is_sparse = sp.isspmatrix(conn_mx)
  if is_sparse:
      # transform to list of lists
      conn_mx = conn_mx.tolil()
  else:
      conn_mx = np.asarray(conn_mx)

  if conn_mx.shape[0] != conn_mx.shape[1]:
      raise ValueError('conn mx should be square')

  DIR = tempfile.gettempdir()
  pfx = str(os.getpid()) + '_'
  NETWORK_FILE = os.path.join(DIR, pfx+'net.txt')
  OUTPUT_FILE  = os.path.join(DIR, pfx+'net.bin')
  NODEMAP_FILE = os.path.join(DIR, pfx+'nmap.bin')
  CONF_FILE    = os.path.join(DIR, pfx+'conf.bin')

  with open(NETWORK_FILE, 'w') as f:
      f.write('>\n')
      for ndx in range(conn_mx.shape[0]):
          f.write('%d %d\n' % (ndx, ndx))
      f.write('>\n0 0\n>\n')
      for ndx in range(conn_mx.shape[0]):
          r = conn_mx[ndx,:]
          if is_sparse:
            conns = r.rows[0]
            weights = r.data[0]
          else:
            conns = np.flatnonzero(r)
            weights = r[conns]

          if not len(conns):
              f.write('%d %d %0.5f 0\n' % (ndx, 0, 0))
          for c,w in zip(conns, weights):
              f.write('%d %d %0.5f 0\n' % (ndx, c, w))

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

  call_opts = [os.path.join(bin_dir,'community'),]
  if rand_init:
    call_opts.append('-r')
  call_opts.append(OUTPUT_FILE)
  call_opts.append(CONF_FILE)
    
  best_membership, best_q = None, None
  for run_ndx in range(num_runs):
    res = subprocess.Popen(call_opts, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    output, errors = map(lambda s: s.decode('ascii'), res.communicate())

    all_files = [NETWORK_FILE, OUTPUT_FILE, NODEMAP_FILE, CONF_FILE]
    if debug:
      print("******* RUN %d *******" % run_ndx)
      print("**** OUTPUT: ****")
      print(output)
      print()
      print("**** STDERR: ****")
      print(errors)
      print("*****************")
      print("Leaving files in place:")
      print(all_files)

    q_value = float(errors.strip().split("\n")[-1].split(" ")[0])
    membership = np.zeros(conn_mx.shape[0], dtype='int')
    ndxs, vals = zip(*[map(int, l.strip().split("\t")) for l in output.strip().split("\n")])
    membership[list(ndxs)] = vals

    if best_q is None or q_value > best_q:
      best_membership, best_q = membership, q_value

  if not debug:
    for fname in all_files:
      os.remove(fname)

  return best_membership, best_q

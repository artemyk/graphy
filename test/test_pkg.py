from __future__ import division, print_function, absolute_import
import six
range = six.moves.range

import networkx as nx
import functools
import numpy as np
import scipy.sparse as sp

import graphy

import nose

def test_partition_search():
	#mx = np.array(nx.to_numpy_matrix(nx.karate_club_graph()),dtype='float')
	
	ground_truth = np.asarray([0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1])
	mxp = graphy.graphgen.get_weighted_block_matrix(ground_truth,0.7,0.1)
	mx = graphy.graphgen.sample_connection_matrix(mxp)
	
	qualityObj = graphy.qualityfuncs.Modularity(mx)

	found_membership = graphy.partitions.find_optimal(qualityObj)

	remapped = graphy.partitions.remap2match(found_membership, ground_truth)
	assert((found_membership == ground_truth).all)


def test_gen():
	mx = graphy.graphgen.gen_hierarchical_weighted_block_matrix(5, 2, 1, [0.5,0.1])
	mx2 = graphy.graphgen.get_weighted_block_matrix([0,0,0,0,0,1,1,1,1,1],0.5,0.1)
	np.testing.assert_allclose(mx, mx2)

	mx3 = graphy.graphgen.gen_hierarchical_weighted_block_matrix(2, 2, 2, [0.3,0.2,0.1])
	trg = np.array(
		[[0.  ,0.3 ,0.2 ,0.2 ,0.1 ,0.1 ,0.1 ,0.1],
		 [0.3 ,0.  ,0.2 ,0.2 ,0.1 ,0.1 ,0.1 ,0.1],
		 [0.2 ,0.2 ,0.  ,0.3 ,0.1 ,0.1 ,0.1 ,0.1],
		 [0.2 ,0.2 ,0.3 ,0.  ,0.1 ,0.1 ,0.1 ,0.1],
		 [0.1 ,0.1 ,0.1 ,0.1 ,0.  ,0.3 ,0.2 ,0.2],
		 [0.1 ,0.1 ,0.1 ,0.1 ,0.3 ,0.  ,0.2 ,0.2],
		 [0.1 ,0.1 ,0.1 ,0.1 ,0.2 ,0.2 ,0.  ,0.3],
		 [0.1 ,0.1 ,0.1 ,0.1 ,0.2 ,0.2 ,0.3 ,0. ]]
	)
	np.testing.assert_allclose(mx3, trg)

def test_louvain():
	def _run_test(conn_mx, groundtruth=None):

		dirMod = graphy.qualityfuncs.DirectedModularity(conn_mx)
		best_membership, q = graphy.louvain.optimize_modularity(conn_mx)

		assert(np.abs( q - dirMod.quality(best_membership) ) < 1e-2)

		if groundtruth is not None:
			best_membership = graphy.partitions.remap2match(best_membership,groundtruth)
			assert(np.array_equal(groundtruth, best_membership))

		best_membership_sp, q_sp = graphy.louvain.optimize_modularity(sp.csr_matrix(conn_mx))
		assert(np.abs( q - q_sp ) < 1e-2)

		assert(np.array_equal(best_membership, best_membership_sp))

	c1 = graphy.graphgen.gen_hierarchical_weighted_block_matrix(5, 2, 2, [1, 0.1, 0.01])

	graphy.louvain.optimize_modularity(c1, num_runs=2)
	graphy.louvain.optimize_modularity(c1, rand_init=False)

	_run_test(c1, np.array([0,0,0,0,0,1,1,1,1,1,2,2,2,2,2,3,3,3,3,3]))

	c1[0:5,:] = 0
	_run_test(c1)

@nose.tools.raises(ValueError)
def test_louvain_error():
	graphy.louvain.optimize_modularity(np.eye(10), num_runs=2, rand_init=False)

def test_plotting_of_gen_graph():
	G = graphy.graphgen.gen_hierarchical_net(5, 2)
	pos = graphy.graphgen.get_hierarchical_net_pos(G)
	graphy.plotting.plot_graph(G, pos=pos)

def test_louvain_randomization():
	randgraph = (np.random.rand(50,50) > 0.7).astype('int')
	_, q1 = graphy.louvain.optimize_modularity(randgraph, rand_init=True)
	_, q2 = graphy.louvain.optimize_modularity(randgraph, rand_init=True)
	assert(q1!=q2)

	_, q3 = graphy.louvain.optimize_modularity(randgraph, rand_init=False)
	_, q4 = graphy.louvain.optimize_modularity(randgraph, rand_init=False)
	assert(q3==q4)

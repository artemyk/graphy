from __future__ import division, print_function, absolute_import
import six
range = six.moves.range

import networkx as nx
import functools
import numpy as np
import scipy.sparse as sp

import graphy


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
	_run_test(c1, np.array([0,0,0,0,0,1,1,1,1,1,2,2,2,2,2,3,3,3,3,3]))

	c1[0:5,:] = 0
	_run_test(c1)

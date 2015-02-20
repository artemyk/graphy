from __future__ import division, print_function, absolute_import
import six
range = six.moves.range

import networkx as nx
import functools
import numpy as np

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
	conn_mx = graphy.graphgen.gen_hierarchical_weighted_block_matrix(5, 2, 2, [1, 0.1, 0.01])

	dirMod = graphy.qualityfuncs.DirectedModularity(conn_mx)
	best_membership, q = graphy.louvain.optimize_modularity(conn_mx)

	assert(np.abs( q - dirMod.quality(best_membership) ) < 1e-2)

	groundtruth = np.array([0,0,0,0,0,1,1,1,1,1,2,2,2,2,2,3,3,3,3,3])
	best_membership = graphy.partitions.remap2match(best_membership,groundtruth)

	assert(np.array_equal(groundtruth, best_membership))

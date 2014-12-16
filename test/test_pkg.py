from __future__ import division, print_function, absolute_import
import six
range = six.moves.range

import networkx as nx
import functools
import numpy as np

import graphy


def test_partition_search_karate():
	#mx = np.array(nx.to_numpy_matrix(nx.karate_club_graph()),dtype='float')
	
	ground_truth = np.asarray([0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1])
	mx = graphy.graphgen.get_block_matrix(ground_truth, 0.7, 0.1)

	qualityObj = graphy.qualityfuncs.Modularity(mx)

	found_membership = graphy.partitions.find_optimal(mx.shape[0], qualityObj)

	remapped = graphy.partitions.remap2match(found_membership, ground_truth)
	assert((found_membership == ground_truth).all)



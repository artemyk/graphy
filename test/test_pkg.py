from __future__ import division, print_function, absolute_import
import six
range = six.moves.range

import networkx as nx
import functools
import numpy as np

import graphy


def modularity(membership, mx, total_stubs, ks):
	q = 0
	for comm in set(membership):
		ixs = membership == comm
		q += mx[ixs,:][:,ixs].sum() - np.outer(ks[ixs], ks[ixs]/total_stubs).sum()
	return q

def test_partition_search_karate():
	mx = np.array(nx.to_numpy_matrix(nx.karate_club_graph()),dtype='float')
	ground_truth = np.asarray([0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1])
	mx = graphy.graphgen.get_block_matrix(ground_truth, 0.7, 0.1)
	total_stubs = mx.sum()
	ks = mx.sum(axis=1)

	modularity_as_func_of_membership = functools.partial(modularity, mx=mx, total_stubs=total_stubs, ks=ks)

	opt_obj = graphy.partitions.FindOptimal(mx.shape[0], modularity_as_func_of_membership)
	found_membership = opt_obj.run()

	remapped = graphy.partitions.remap2match(found_membership, ground_truth)
	assert((found_membership == ground_truth).all)



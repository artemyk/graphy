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

    ground_truth = np.asarray([0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1])
    mxp = graphy.graphgen.get_weighted_block_matrix(ground_truth,0.7,0.1)
    mx = graphy.graphgen.sample_connection_matrix(mxp)

    for qualityObj in [graphy.qualityfuncs.Modularity(mx), 
                       graphy.qualityfuncs.InfoMapCodeLength(mx)]:
        found_membership, q = qualityObj.find_optimal()
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

def test_louvain_sparse():
    """
    Test whether optimize modularity returns same results if passing sparse/dense matrix
    """
    from numpy.random import seed
    A = sp.rand(100, 100, .1, 'lil')
    Ad = A.todense()
    m1, q1 = graphy.louvain.optimize_modularity(A, rand_init=False)
    m2, q2 = graphy.louvain.optimize_modularity(Ad, rand_init=False)
    assert np.allclose(m1, m2), "memberships differ"
    assert np.isclose(q1, q2), "modularity differs"


def test_modularity_vals():
    """
    Test returned modularity values
    """
    G=nx.to_numpy_matrix(nx.karate_club_graph())
    membership = np.zeros(G.shape[0])

    # Communities reported by Zachary
    membership[[0,1,2,3,4,5,6,7,10,11,12,13,16,17,19,21]] = 1

    qObj = graphy.qualityfuncs.Modularity(G)
    q = qObj.quality(membership)
    
    # 0.371 value reported in:    
    # Donetti, Munoz, Detecting Network Communities: a new systematic 
    #   and efficient algorithm, 2004, http://arxiv.org/pdf/cond-mat/0404652.pdf
    assert(np.isclose(q, 0.371, atol=5e-4))

    # Directed modularity should still work for undirected graphs
    qObj2 = graphy.qualityfuncs.DirectedModularity(G)
    q2 = qObj2.quality(membership)
    assert(np.isclose(q2, 0.371, atol=5e-4))

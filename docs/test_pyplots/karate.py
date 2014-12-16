import numpy as np
import networkx as nx
import graphy

mx = np.array(nx.to_numpy_matrix(nx.karate_club_graph()))

qualityObj = graphy.qualityfuncs.Modularity(mx)

best_membership = graphy.partitions.find_optimal(mx.shape[0], qualityObj)




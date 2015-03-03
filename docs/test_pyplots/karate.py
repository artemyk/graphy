import numpy as np
import networkx as nx
import graphy

G = nx.karate_club_graph()
qualityObj = graphy.qualityfuncs.Modularity(nx.to_numpy_matrix(G))
best_membership, q = qualityObj.find_optimal()
graphy.plotting.plot_graph(G, pos=nx.spring_layout(G), colors=best_membership)


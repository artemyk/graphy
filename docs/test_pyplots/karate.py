mx = np.array(nx.to_numpy_matrix(nx.karate_club_graph()))

qualityObj = graphy.costfunctions.Modularity(mx)

opt_obj = graphy.partitions.FindOptimal(mx.shape[0], qualityObj)

found_membership = opt_obj.run()




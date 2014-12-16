import numpy as np


class CostFunction(object):
	def __init__(self):
		# precompute quantities here
		pass

	def cost(self, membership):
		raise NotImplementedError
		
class Modularity(CostFunction):
	def __init__(self, mx):
		self.mx = mx.astype('float')
		self.total_stubs = self.mx.sum()
		self.ks = self.mx.sum(axis=1)

	def quality(self, membership):
		q = 0
		for comm in set(membership):
			ixs = membership == comm
			q += self.mx[ixs,:][:,ixs].sum() - np.outer(self.ks[ixs], self.ks[ixs]/self.total_stubs).sum()
		return q


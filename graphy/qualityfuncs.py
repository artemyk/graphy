"""Module implements quality functions for graph decompositions.
"""

from __future__ import division, print_function, absolute_import
import six
range = six.moves.range

import numpy as np


class QualityFunction(object):
	def __init__(self):
		"""Implements a quality function.
		"""
		# precompute quantities here
		pass

	def quality(self, membership):
		"""Returns a quality score corresponding to membership vector. 

		Parameters
		----------
		membership : np.array
			Membership vector

		Returns
		-------
		float
			Quality

		"""
		raise NotImplementedError
		
class Modularity(QualityFunction):
	def __init__(self, mx):
		"""Class that implements Newman's modularity quality function.
		"""
		self.mx = mx.astype('float')
		self.total_stubs = self.mx.sum()
		self.ks = self.mx.sum(axis=1)

	def quality(self, membership):
		q = 0
		for comm in set(membership):
			ixs = membership == comm
			q += self.mx[ixs,:][:,ixs].sum() - np.outer(self.ks[ixs], self.ks[ixs]/self.total_stubs).sum()
		return q


# -*- coding: utf-8 -*-

""" Functions for the management of vectors """

import numpy as np

def uniquec(l):
	""" Counts the instances of the uniqued integers in l. """

	# Possible integer values
	possible = range(max(l) + 1)

	# Count elements occurrences
	counts = [0 for i in possible]
	for n in l:
		counts[n] += 1

	# Return tupled counts
	return [(i, counts[i]) for i in possible if counts[i]]

def flatten_and_select(v, s):
	"""
	Flatten the array v and select rows based on s.

	@param:
	 - v <np.array>
	 - s <list> list of v rows indexes
	"""

	v = v.reshape((np.prod(v.shape))).tolist()
	v = [v[i] for i in s]
	return(v)

def rm_from_mask(L, torm):
	"""
	Remove elements from a mask.

	@param:
	 - L <np.array[int]> labelled objects
	 - torm <list> list of objects indexes (to remove)
	"""

	if len(torm) <= L.max() - len(torm):
		# Update list of objects to be discarded
		torm = [e + 1  for e in torm]

		# Identify which objects to discard
		rm_mask = np.vectorize(lambda x: x in torm)(L)

		# Discard and re-label
		L[rm_mask] = 0
	else:
		# Select objects to be kept
		tokeep = [e + 1 for e in range(L.max()) if not e in torm]

		# Identify which objects to discard
		rm_mask = np.vectorize(lambda x: not x in tokeep)(L)

		# Discard and re-label
		L[rm_mask] = 0

	# Output
	return(L > 0)

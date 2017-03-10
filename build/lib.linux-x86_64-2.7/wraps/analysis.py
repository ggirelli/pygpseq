# -*- coding: utf-8 -*-

from .. import const
from ..tools import path as pt
from ..tools import io as iot

class Analyzer(iot.ioInterface):
	""" GPSeq image analyzer """

	# Class-bound attributes
	__version__ = const.VERSION

	# Channel names
	dna_names = ('dapi')
	sig_names = ('tmr', 'cy5')

	# Adaptive threshold combination
	# If even, increased by one
	# If <= 1, turned off
	adp_thr = 101

	# Voxel ratio
	aspect = (1., 1., 1.)

	# Nuclear radius interval
	radius_interval = (10., float('inf'))

	# Minimum Z size of the nucleus
	# If greater than 1, the ceil is taken as the minimum number of slices
	# If lower than 1, the float is taken as the fraction of the stack
	min_z_size = .25

	# Segmentation/Analysis type
	seg_type = const.SEG_SUM_PROJ
	an_type = const.AN_SUM_PROJ

	# Nuclear selection features
	# Uses at least two, up to three, features
	# The first two features are used in the summary scatterplot, order matters
	nsf = (const.NSEL_SIZE, const.NSEL_SUMI, const.NSEL_SHAPE)

	# Bounding box [Z Y X]
	offset = (0, 5, 5)

	# Partial nucleus erosion
	part_n_erosion = .5

	# Perform nuclear surface calculation
	# Otherwise will store the size in the surface field
	calc_n_surface = False

	# Sigma for smoothing
	sigma = .1

	# Number of bins for profile calculation
	nbins = 200

	# Remove nuclei touching the top/bottom of the stack
	rm_z_tips = False

	# Rescale deconvolved images?
	rescale_deconvolved = False

	# Perform chromatic aberration correction
	correctCA = False

	# Use normalized distance for profile study
	normalize_distance = True

	# Dictionary condition-description
	cdescr = {}

	def __init__(self):
		super(Analyzer, self).__init__()

	def check_anseg_types(self):
		""" Checks seg_type and an_type """

		# 2D segmentation does not allow 3D analysis
		no3d_cond = self.an_type == const.AN_3D
		no3d_cond = no3d_cond and self.seg_type != const.SEG_3D
		if no3d_cond:
			# Revert analysis to default
			msg = '3D analysis is not available for 2D segmentation.\n'
			msg += 'Using sum z-projection analysis instead...'
			self.printout(msg, -1)
			self.an_type = const.SEG_SUM_PROJ

		# 2D segmentation does not allow mid-section analysis
		nomid_cond = self.an_type == const.AN_MID
		nomid_cond = nomid_cond and self.seg_type != const.SEG_3D
		if nomid_cond:
			# Revert analysis to default
			msg = 'Mid-section analysis is not available for 2D segmentation.\n'
			msg += 'Using sum z-projection analysis instead...'
			self.printout(msg, -1)
			self.an_type = const.SEG_SUM_PROJ

	def __setattr__(self, name, value):
		""" Check the attribute and set it. """

		# Check the attribute
		check = self.check_attr(name, value)

		if True == check:
			# Set the attribute
			return(super(Analyzer, self).__setattr__(name, value))
		else:
			# Don't set the attribute
			return(None)

	def check_attr(self, name, value):
		""" Check attribute format and value. """

		# Default answer
		checked = True

		if name in ['dna_names', 'sig_names']:
			# Require tuple
			if not type(()) == type(value):
				checked = False
			# Require non-empty tuple
			elif 0 == len(value):
				checked = False
			else:
				# Require tuple of strings
				types = [type('') == type(s) for s in value]
				types = all(types)
				if not types:
					checked = False

			if not checked:
				msg = '"' + name + '" must be a non-empty tuple of strings.\n'
				msg += 'Keeping previous value.'
				self.printout(msg, -1)
		elif name in ['aspect', 'radius_interval']:
			# Require tuple
			if not type(()) == type(value):
				checked = False
			# Require non-empty tuple
			elif 0 == len(value):
				checked = False
			else:
				# Require tuple of floats
				types = [type(1.) == type(s) for s in value]
				types = all(types)
				if not types:
					checked = False

			if not checked:
				msg = '"' + name + '" must be a non-empty tuple of floats.\n'
				msg += 'Keeping previous value.'
				self.printout(msg, -1)
		elif 'min_z_size' == name:
			# Allow float or integer
			if not type(value) in [type(0), type(0.)]:
				checked = False

			if not checked:
				msg = '"' + name + '" must be a float lower than 1 or'
				msg += ' an integer greater than 1.\n'
				msg += 'Keeping previous value.'
				self.printout(msg, -1)
		elif 'seg_type' == name:
			# Check that it is one of the allowed constants
			seg_types = [const.SEG_SUM_PROJ, const.SEG_MAX_PROJ, const.SEG_3D]
			if not value in seg_types:
				checked = False

			if not checked:
				msg = '"' + name + '" must be one of the following values:\n'
				msg += str(seg_types) + '\n'
				msg += 'Keeping previous value.'
				self.printout(msg, -1)
		elif 'an_type' == name:
			# Check that it is one of the allowed constants
			an_types = [const.AN_SUM_PROJ, const.AN_MAX_PROJ,
				const.AN_3D, const.AN_MID]
			if not value in an_types:
				checked = False

			if not checked:
				msg = '"' + name + '" must be one of the following values:\n'
				msg += str(an_types) + '\n'
				msg += 'Keeping previous value.'
				self.printout(msg, -1)
		elif 'nsf' == name:
			# Require tuple
			if not type(()) == type(value):
				checked = False
			# Require non-empty tuple
			elif 0 == len(value):
				checked = False
			# Check that it is one of the allowed constants
			elif not value in range(len(const.NSEL_FIELDS)):
					checked = False

			if not checked:
				msg = '"' + name + '" must be a non-empty tuple with at least 2'
				msg += ' of the following values:\n'
				msg += str(range(len(const.NSEL_FIELDS))) + '\n'
				msg += 'Keeping previous value.'
				self.printout(msg, -1)
		elif 'offset' == name:
			# Require tuple
			if not type(()) == type(value):
				checked = False
			# Require non-empty tuple
			elif 0 == len(value):
				checked = False
			else:
				# Require tuple of integers
				types = [type(0) == type(s) for s in value]
				types = all(types)
				if not types:
					checked = False

			if not checked:
				msg = '"' + name + '" must be a non-empty tuple of integers.\n'
				msg += 'Keeping previous value.'
				self.printout(msg, -1)
		elif 'sigma' == name:
			# Require float
			if not type(1.) == type(value):
				checked = False
			# Require positive float
			elif 0 > value:
					checked = False

			if not checked:
				msg = '"' + name + '" must be a positive float.\n'
				msg += 'Keeping previous value.'
				self.printout(msg, -1)
		elif name in ['rm_z_tips', 'rescale_deconvolved', 'correctCA']:
			# Require boolean
			if not type(True) == type(value):
				checked = False
			
				msg = '"' + name + '" must be a Boolean.\n'
				msg += 'Keeping previous value.'
				self.printout(msg, -1)
		elif 'cdescr' == name:
			# Require a dictionary
			if not type({}) == type(value):
				checked = False
			else:
				# Require a dictionary with string values
				types = [type('') == type(v) for v in value.values()]
				types = all(types)
				if not types:
					checked = False

				msg = '"' + name + '" must be a dictionary with string values.\n'
				msg += 'Keeping previous value.'
				self.printout(msg, -1)

		# Output
		return(checked)

# -*- coding: utf-8 -*-

"""Image binarization module.

This module contains code to binarize images.

"""

import math

import numpy as np
from scipy import ndimage as ndi
from skimage.filters import threshold_otsu
from skimage.measure import label

from .. import const

from . import image as imt
from . import io as iot
from . import stat as stt
from . import vector as vt

class Binarize(iot.IOinterface):
	"""Image binarization class.

	If the class has public attributes, they may be documented here
	in an ``Attributes`` section and follow the same formatting as a
	function's ``Args`` section. Alternatively, attributes may be documented
	inline with the attribute's declaration (see __init__ method below).

	Properties created with the ``@property`` decorator should be documented
	in the property's getter method.

	Attributes:
		an_type (int): analysis type accordin to `pygpseq.const`.
		seg_type (int): segmentation type accordin to `pygpseq.const`.
		do_global_thr (bool): True to apply global threshol (Otsu).
		do_adaptive_thr (bool): True to apply local (adaptive) threshold.
		adaptive_neighbourhood (int): neighbourhood square side for local thr.
		do_clear_borders (bool): True to remove objects touching the borders.
		do_clear_Z_borders (bool): True to remove objects touching Z borders.
		do_fill_holes (bool): True to fill holes (both 2D and 3D).
		radius_interval (tuple[float]): object radius interval.
		min_z_size (float): minimum Z (relative) size of the objects.

	"""

	an_type = 0
	seg_type = 0
	do_global_thr = True
	do_adaptive_thr = True
	adaptive_neighbourhood = 101
	do_clear_borders = True
	do_clear_Z_borders = False
	do_fill_holes = True
	radius_interval = (10., float('inf'))
	min_z_size = .25

	def __init__(self, **kwargs):
		"""Initialize binarization settings all at once with kwargs.

		Args:
			**kwargs: arbitrary keyword arguments stored in the class.
		"""

		# Run IOinterface __init__ method
		super(Binarize, self).__init__()

		# Store provided kwargs in the current instance.
		for k in kwargs.keys():
			if not k == 'logpath':
				self[k] = kwargs[k]

	def run(self, im):
		"""Binarize image with current instance settings.

		Args:
			im (np.array): image to be thresholded

		Returns:
			tuple: binarized image, Otsu's threshold value and log string
		"""

		log = ''

		# If no threshold is requested, return the image
		if not self.do_global_thr and not self.do_adaptive_thr:
			msg = 'No threshold applied.'
			log += self.printout(msg, -1)

			return((im, log))

		# Make Z-projection ----------------------------------------------------
		if const.SEG_3D != self.seg_type and 2 != len(im.shape):
			msg = 'Generating Z-projection ['
			msg += const.SEG_LABELS[self.seg_type] + ']...'
			log += self.printout(msg, 2)
			im = imt.mk_z_projection(im, self.seg_type)

		# Binarize images ------------------------------------------------------
		mask = []

		# Perform global threshold
		if self.do_global_thr:
			thr = threshold_otsu(im)
			msg = 'Thresholding image, global thr: ' + str(thr)
			log += self.printout(msg, 2)
			mask.append(imt.binarize(im, thr))

		# Perform adaptive threshold
		if self.do_adaptive_thr and 1 < self.adaptive_neighbourhood:
			msg = 'Applying adaptive threshold to neighbourhood: '
			msg += str(self.adaptive_neighbourhood)
			log += self.printout(msg, 2)
			mask.append(imt.threshold_adaptive(im, self.adaptive_neighbourhood))

		# Combine masks
		if len(mask) == 2:
			mask = np.logical_and(mask[0], mask[1])
		else:
			mask = mask[0]

		# Remove objects touching borders --------------------------------------
		if self.do_clear_borders:
			msg = 'Removing objects touching the image border...'
			log += self.printout(msg, 2)
			mask = imt.clear_borders(mask, self.do_clear_Z_borders)
		
		# Fill holes -----------------------------------------------------------
		if self.do_fill_holes:
			log += self.printout('Filling holes...', 2)
			mask = ndi.binary_fill_holes(mask)

			if 3 == len(mask.shape):
				# Single slice filling
				for sliceid in range(mask.shape[0]):
					slide = mask[sliceid, :, :]
					mask[sliceid, :, :] = ndi.binary_fill_holes(slide)

		# Output ---------------------------------------------------------------
		return((mask, thr, log))
	
	def filter_obj_XY_size(self, mask):
		"""Filter objects XY size.

		Note:
			Uses self.radius_interval to filter the objects in the provided
			mask based on the selected segmentation type.

		Args:
			mask (np.array): binary image

		Returns:
			tuple: filtered binary image and log string
		"""

		# Start logging
		log = ''
		log += self.printout('Filtering objects XY size...', 2)

		# From radius to size
		sinter = stt.r_to_size(self.radius_interval, self.seg_type)
		msg = 'Allowed size interval: '
		msg += '[' + str(round(sinter[0])) + ', '
		msg += str(round(sinter[1])) + '] [' + imt.get_unit(mask.shape) + ']'
		log += self.printout(msg, 3)

		# Identify objects XY size
		log += self.printout('Retrieving objects XY size...', 3)
		L = label(mask)
		xysizes = imt.get_objects_xysize(L)
		log += self.printout('Found ' + str(L.max()) + ' objects.', 3)
		
		# Select objects to be discarded
		torm = np.logical_or(xysizes < sinter[0], xysizes > sinter[1])
		torm = [ii for ii, x in enumerate(torm) if x]
		log += self.printout('Discarding ' + str(len(torm)) + ' objects.', 3)

		# Remove objects outside of size interval
		mask = vt.rm_from_mask(L, torm)
		L = label(mask)

		# Output
		return((mask, log))

	def filter_obj_Z_size(self, mask):
		"""Filter objects Z size.

		Note:
			Uses self.min_z_size to filter the objects in the provided mask.

		Args:
			mask (np.array): binary image

		Returns:
			tuple: filtered binary image and log string
		"""

		# Start logging
		log = ''
		log += self.printout('Filtering objects Z size...', 2)

		# If not a stack, return the mask
		if 3 > len(mask.shape):
			return((mask, log))

		# Check provided conditions
		doFilterZsize = 0 != int(math.ceil(self.min_z_size))
		doFilterZsize = doFilterZsize and self.an_type == const.AN_3D
		if not doFilterZsize:
			return((mask, log))

		# From size to number of slices
		if self.min_z_size > 1:
			self.min_z_size = int(math.ceil(self.min_z_size))
		else:
			self.min_z_size = self.min_z_size * mask.shape[0]
			self.min_z_size = int(math.ceil(self.min_z_size))
		msg = 'Minimum ' + str(self.min_z_size) + ' slices.'
		log += self.printout(msg, 3)

		# Identify objects Z size
		log += self.printout('Retrieving objects Z size...', 3)
		L = label(mask)
		zsizes = imt.get_objects_zsize(L)
		log += self.printout('Found ' + str(L.max()) + ' objects.', 3)
		
		# Select objects to be discarded
		torm = np.array(zsizes) < self.min_z_size
		torm = [ii for ii, x in enumerate(torm) if x]
		msg = 'Discarding ' + str(len(torm)) + ' objects.'
		log += self.printout(msg, 3)

		# Remove objects lower than minimum size
		mask = vt.rm_from_mask(L, torm)
		L = label(mask)

		# Output
		return((mask, log))

	def __getitem__(self, key):
		""" Allow get item. """
		if key in dir(self):
			return(getattr(self, key))
		else:
			return(None)

	def __setitem__(self, key, value):
		""" Allow set item. """
		if key in dir(self):
			self.__setattr__(key, value)

	def __setattr__(self, name, value):
		""" Check the attribute and set it. """

		# Check the attribute
		check = self.check_attr(name, value)

		if True == check:
			# Set the attribute
			return(super(Binarize, self).__setattr__(name, value))
		else:
			# Don't set the attribute
			return(None)

	def check_attr(self, name, value):
		"""Check attribute format and value.

		Args:
			name (string): attribute name
			value: attribute value

		Returns:
			bool: whether the provided (name,value) couple passed its
				name-specific test
		"""

		# Default answer
		checked = True

		if name in ['do_global_thr', 'do_adaptive_thr', 'do_clear_borders',
			'do_clear_Z_borders', 'do_fill_holes']:
			# Require boolean
			if not type(True) == type(value):
				checked = False
			
				msg = '"' + name + '" must be a Boolean.\n'
				msg += 'Keeping previous value [' + str(self[name]) + '].'
				self.printout(msg, -1)

		elif name == 'adaptive_neighbourhood':
			# Require boolean
			if not type(0) == type(value):
				checked = False
			
				msg = '"' + name + '" must be an integer.\n'
				msg += 'Keeping previous value [' + str(self[name]) + '].'
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
				msg += 'Keeping previous value [' + str(self[name]) + '].'
				self.printout(msg, -1)

		elif 'seg_type' == name:
			# Check that it is one of the allowed constants
			seg_types = [const.SEG_SUM_PROJ, const.SEG_MAX_PROJ, const.SEG_3D]
			if not value in seg_types:
				checked = False

			if not checked:
				msg = '"' + name + '" must be one of the following values:\n'
				msg += str(seg_types) + '\n'
				msg += 'Keeping previous value [' + str(self[name]) + '].'
				self.printout(msg, -1)

		return(checked)

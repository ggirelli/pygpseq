# -*- coding: utf-8 -*-

import numpy as np
from scipy import ndimage as ndi
from scipy.ndimage.morphology import distance_transform_edt
import skimage.io as io
from skimage.measure import label, marching_cubes, mesh_surface_area

from .. import const
from ..tools import io as iot, image as imt, stat as stt, string as st
from ..tools import vector as vt

class Nucleus(iot.ioInterface):
	"""docstring for nucleus"""

	# Class-bound attributes
	__version__ = const.VERSION

	# Series id (1-indexed)
	s = 0

	# Nucleus id (1-indexed)
	n = 0

	# Box corner coordinates
	box = []

	# Voxel ratio
	aspect = [1, 1, 1]

	# Background
	dna_bg = 0
	sig_bg = 0

	# Size in px or vx
	flat_size = 0
	size = 0
	unit = ''

	# Nuclear surface
	surf = 0

	# Sum/Mean intensity
	sumI = 0
	meanI = 0

	# Circularity/sphearicity
	shape = 0

	# Threshold value
	thr = 0

	def __init__(self, logpath, n, series_id, mask, i, thr, offset, aspect,
		dna_bg, sig_bg, calc_n_surface = None, **kwargs):
		"""
		@param:
		 - logpath <string> path to the log file
		 - n <int> nucleus id (1-indexed)
		 - series_id <int> series id (1-indexed)
		 - mask <numpy.array[boolean]> binary image
		 - i <numpy.array[uint16]> image
		 - thr <uint16> threshold obtained with Otsu's method
		 - offset <tuple[int]> dimensions box/square offset
		 - aspect <tuple[float]> pixel/voxel dimension proportion
		 - dna_bg <uint16> median background for DNA channel
		 - sig_bg <uint16> median background for Signal channel
		"""
		
		# parent class __init__
		super(Nucleus, self).__init__(path = logpath, append = True)

		# Default values
		if None == calc_n_surface:
			calc_n_surface = True

		# Store parameters locally
		self.s = series_id
		self.n = n
		self.box = self.get_bounding_box(mask, offset)
		self.thr = thr
		self.dna_bg = dna_bg
		self.sig_bg = sig_bg
		self.aspect = aspect

		# Apply box selection to the image
		i = imt.apply_box(i, self.box)
		mask = imt.apply_box(mask, self.box)

		# Nuclear measurements
		self.size = mask.sum()
		if 3 == len(i.shape):
			self.flat_size = mask.max(0).sum()
		else:
			self.flat_size = self.size
		self.unit = imt.get_unit(i.shape)
		self.sumI = i[mask == 1].sum()
		self.meanI = self.sumI / self.size
		self.shape = imt.describe_shape(mask, self.aspect)
		if 3 == len(mask.shape) and calc_n_surface:
			self.surf = imt.calc_surface(mask, self.aspect)
		else:
			self.surf = self.size

	def get_bounding_box(self, mask, offset = None):
		"""
		Return the bounding box (2d or 3d) of the object in mask.
		An offset can be specified for each dimension.
		If no offset is specified, it defaults to 0.
		If only one offset is specified, it will be used for every dimension.
		If the number of offsets specified does not match the number of
		dimensions, onlyt the first will be used for every dimension.
		"""

		if 2 == len(mask.shape):
			return(self.get_2d_bounding_box(mask, offset))
		elif 3 == len(mask.shape):
			return(self.get_3d_bounding_box(mask, offset))

	def get_2d_bounding_box(self, mask, offset = None):
		"""
		Return the bounding box (2d) of the object in mask.
		An offset can be specified for each dimension.
		If no offset is specified, it defaults to 0.
		If only one offset is specified, it will be used for every dimension.
		If the number of offsets specified does not match the number of
		dimensions, onlyt the first will be used for every dimension.
		"""

		# Check provided offset
		offset = self.check_box_offset(mask.shape, offset)

		# Binarize mask if it is not
		mask = mask.astype('bool').astype('uint8')

		box = []

		# Y-side boundaries
		vy = mask.max(1).tolist()
		vy_min = max(0, vy.index(1) - offset[1])
		vy.reverse()
		vy_max = min(mask.shape[1] - 1, len(vy) - vy.index(1) - 1 + offset[1])
		box.append((vy_min, vy_max))

		# X-side boundaries
		vx = mask.max(0).tolist()
		vx_min = max(0, vx.index(1) - offset[0])
		vx.reverse()
		vx_max = min(mask.shape[0] - 1, len(vx) - vx.index(1) - 1 + offset[0])
		box.append((vx_min, vx_max))

		return(box)

	def get_3d_bounding_box(self, mask, offset = None):
		"""
		Return the bounding box (3d) of the object in mask.
		An offset can be specified for each dimension.
		If no offset is specified, it defaults to 0.
		If only one offset is specified, it will be used for every dimension.
		If the number of offsets specified does not match the number of
		dimensions, onlyt the first will be used for every dimension.
		"""

		# Check provided offset
		offset = self.check_box_offset(mask.shape, offset)

		# Binarize mask if it is not
		mask = mask.astype('bool').astype('uint8')

		# Retrieve 2D bounding box
		box = [()]
		box.extend(self.get_2d_bounding_box(mask.max(0), offset[1:2]))

		# Z-side boundaries
		vz = mask.max(1).max(1).tolist()
		vz_min = max(0, vz.index(1) - offset[0])
		vz.reverse()
		vz_max = min(mask.shape[0] - 1, len(vz) - vz.index(1) - 1 + offset[0])
		box[0] = (vz_min, vz_max)

		return(box)

	def check_box_offset(self, shape, offset = None):
		"""
		Check bounding box offset.
		If no offset is specified, it defaults to 0.
		If only one offset is specified, it will be used for every dimension.
		If the number of offsets specified does not match the number of
		dimensions, onlyt the first will be used for every dimension.
		"""

		if None == offset:
			offset = 0

		# Make offset into a list
		if type([]) != type(offset):
			offset = list(offset)

		# Identify the offset for every dimension
		if len(offset) != len(shape):
			offset = [offset[0] for d in shape]

		return(offset)

	def get_data(self, dna_ch, sig_ch, an_type, aspect, debugging,
		part_n_erosion, **kwargs):
		"""
		Get nuclear data.

		@param:
		 - dna_ch <numpy.array> image (dimensionality based on an_type)
		 - sig_ch <numpy.array> image (dimensionality based on an_type)
		 - an_type <pyGPSeq.const> analysis type
		 - aspect <tuple[float]> pixel/voxel dimension proportion
		"""

		# Set output suffix
		if not 'suffix' in kwargs.keys():
			suffix = ''
		else:
			suffix = st.add_leading_dot(kwargs['suffix'])

		# Set plotting
		if not 'plotting' in kwargs.keys():
			kwargs['plotting'] = True

		# Start log
		log = ""
		msg = 'Retrieving single-pixel data for nucleus #' + str(self.n) + '...'
		log += self.printout(msg, 3)

		# Apply box selection to channels
		dna = imt.apply_box(dna_ch, self.box)
		sig = imt.apply_box(sig_ch, self.box)

		if not 'mask' in kwargs.keys():
			# Re-build mask
			mask = dna >= self.thr
		else:
			mask = imt.apply_box(kwargs['mask'], self.box)

		# Fill mask holes
		mask = ndi.binary_fill_holes(mask)
		if 3 == len(mask.shape):
			for sliceid in range(mask.shape[0]):
				mask[sliceid, :, :] = ndi.binary_fill_holes(mask[sliceid, :, :])

		# Select largest object only
		L = label(mask)
		if 1 < L.max():
			sizes = imt.get_objects_xysize(L)
			mask = L == sizes.index(max(sizes)) + 1
		elif 0 == L.max():
			msg = 'Found empty nucleus'
			msg += ' [' + str(self.s) + '.' + str(self.n) + '].'
			self.printout(msg, -1)

		# Apply mask to boxes
		dna[mask == 0] = 0
		sig[mask == 0] = 0

		if const.AN_MID == an_type:
			# Identify middle section
			mid = [dna[sliceid, :, :].sum() for sliceid in range(mask.shape[0])]
			mid = mid.index(max(mid))

			# Select only mid-section
			mask = mask[mid, :, :]
			dna = dna[mid, :, :]
			sig = sig[mid, :, :]

		# Perform distance transform
		D = distance_transform_edt(mask, aspect[3-len(mask.shape):])

		# Export single-nucleus images in debugging mode
		if debugging:
			fname = kwargs['out_dir'] + const.OUTDIR_DEBUG
			fname += 's' + str(self.s) + 'n' + str(self.n)
			if kwargs['plotting']:
				io.imsave(fname + suffix + '.tif', mask.astype('u4'))
			if kwargs['plotting']:
				io.imsave(fname + '.dist' + suffix + '.tif', D.astype('u4'))
			if kwargs['plotting']:
				io.imsave(fname + '.dna' + suffix + '.tif', dna.astype('u4'))
			if kwargs['plotting']:
				io.imsave(fname + '.sig' + suffix + '.tif', sig.astype('u4'))

		# Select pixels for partial 3D nuclear analysis
		sm = np.zeros(mask.shape, dtype = 'u4')
		if const.AN_3D == an_type:
			# Get selection mask
			sm, e = imt.get_partial_nuclear_volume(mask, dna, part_n_erosion)

			# Log error message
			if not '' == e:
				log += self.printout(e, -1)
		
		# Convert image into a list
		mask_flat = mask.reshape([np.prod(mask.shape)])
		mask_flat = mask_flat.tolist()
		mask_flat = [i for i in range(len(mask_flat)) if 1 == mask_flat[i]]

		# Prepare output
		data = np.zeros(len(mask_flat), dtype = const.DTYPE_NUCLEAR_DATA)

		# Flatten data for export
		data['dna'] = vt.flatten_and_select(dna, mask_flat)
		data['sig'] = vt.flatten_and_select(sig, mask_flat)
		data['d'] = vt.flatten_and_select(D, mask_flat)
		data['part'] = vt.flatten_and_select(sm, mask_flat)
		data['n'] = [self.n for i in data['dna']]

		# Remove background
		data['dna'] = np.array(data['dna'])
		data['dna'][data['dna'] < self.dna_bg] = self.dna_bg
		data['dna'] = np.array(data['dna']) - self.dna_bg
		data['sig'] = np.array(data['sig'])
		data['sig'][data['sig'] < self.sig_bg] = self.sig_bg
		data['sig'] = np.array(data['sig']) - self.sig_bg

		# Add normalized distance
		data['dnorm'] = data['d'] - min(data['d'])
		if not 0. == max(data['dnorm']):
			data['dnorm'] = data['dnorm'] / max(data['dnorm'])
		else:
			msg = 'Found 0 maximum distance at '
			msg += str(self.s) + '.' + str(self.n) + '\n\n'
			log += self.printout(msg, -1)

		# Output
		return((data, log))

	def get_summary(self):
		""" Get nuclear summary. """

		# Output
		data = np.array(
			(self.s, self.n, self.flat_size, self.size, self.surf, self.sumI,
			self.meanI, self.shape), dtype = const.DTYPE_NUCLEAR_SUMMARY)

		return(data)

	def export(self, **kwargs):
		""" Export nuclear data. """

		# Set output suffix
		if not 'suffix' in kwargs.keys():
			suffix = ''
		else:
			suffix = st.add_leading_dot(kwargs['suffix'])

		# Get nuclear data
		data, log = self.get_data(**kwargs)

		# Export as csv file
		out_fname = kwargs['series_name'] + '.nucleus' + str(self.n)
		out_fname += suffix + '.csv'
		np.savetxt(kwargs['out_dir'] + fname, data,
			header = ",".join([h for h in data.dtype.names]),
			delimiter = ',', comments = '')

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

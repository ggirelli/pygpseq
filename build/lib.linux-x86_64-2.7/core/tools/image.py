## -*- coding: utf-8 -*-

""" Functions for the management of images """

import os

import numpy as np
from scipy.ndimage.morphology import distance_transform_edt
from skimage import filters
from skimage.measure import label, marching_cubes, mesh_surface_area
from skimage.morphology import closing, convex_hull_image, cube
from skimage.morphology import square
from skimage.segmentation import clear_border

from . import vector as vt
from .. import const

def autoselect_time_frame(im):
	""" Selects the first non-empty time frame found. """

	if 4 == len(im.shape):
		if im.shape[0] == 1:
			return(im)
		
		selected = None

		if 4 == len(im.shape):
			for i in range(im.shape[3]):
				if 0 != im[:, :, :, i].max():
					selected = i
					break

		return(im[:, :, :, selected])
	else:
		return(im)

def mk_z_projection(i, p_type):
	"""
	Make Z-projection. Allowed p_types: SUM and MAX.
	Otherwise return the original image.
	"""
	if const.SUM_PROJ == p_type:
		i = i.sum(0).astype(i.dtype)
	elif const.MAX_PROJ == p_type:
		i = i.max(0).astype(i.dtype)
	return(i)

def binarize(i, thr):
	""" Binarizes and image using the provided threshold """
	if 2 == len(i.shape):
		i = closing(i > thr, square(3))
	elif 3 == len(i.shape):
		i = closing(i > thr, cube(3))
	return(i)

def clear_borders(i, clean_z = None):
	""" Remove objects touching the borders of the image """
	if 2 == len(i.shape):
		i = clear_border(i)
	elif 3 == len(i.shape):
		for slide_id in range(i.shape[0]):
			i[slide_id, :, :] = clear_border(i[slide_id, :, :])
		if True == clean_z:
			for slide_id in range(i.shape[1]):
				i[:, slide_id, :] = clear_border(i[:, slide_id, :])
	return(i)

def threshold_adaptive(i, block_size):
	"""
	Adaptive threshold.

	@param:
	 - i <np.array> image
	 - block_size <int> neighbourhood, if even is incremented
	"""

	# Increment neighbourhood size
	if 0 == block_size % 2:
		block_size += 1

	# Local threshold mask
	lmask = np.zeros(i.shape)

	# Apply threshold per slice
	if 3 == len(i.shape):
		for slice_id in range(i.shape[0]):
			lmask[slice_id, :, :] = filters.threshold_adaptive(
				i[slice_id, :, :], block_size)
		lmask = closing(lmask, cube(3))
	else:
		lmask = filters.threshold_adaptive(i, block_size)
		lmask = closing(lmask, square(3))

	# Output
	return(lmask)

def get_unit(shape):
	""" Returns "vx" for 3d images, "px" for 2d images """
	if 2 == len(shape):
		return('px')
	elif 3 == len(shape):
		return('vx')
	else:
		return('')

def apply_box(i, box):
	""" Applies square/box selection to an image. """

	# Check box
	box = check_box(i.shape, box)

	# Apply box
	return(i[np.ix_(*[tuple(range(t[0], t[1]+1)) for t in box])])

def check_box(shape, box):
	"""
	Check if a square/box selection can be applied to an image.
	If not enough corners are specified, the whole dimension is selected.
	"""

	# Every dimension in the box should have two coordinates
	if len(shape) >= len(box):
		for i in range(len(box)):
			if len(box[i]) != 2:
				box[i] = (0, shape[i] - 1)

	# Fill missing dimensions
	if len(shape) > len(box):
		[box.extend([0, shape[d]]) for d in range(int(len(box)/2), len(shape))]

	# Remove extra dimensions
	box = box[:len(shape)]

	return(box)

def describe_shape(mask, spacing = None):
	"""
	Calculates sphericity (3d) or solidity (2d) of the provided mask.
	Mask is a binary image with only one object.
	"""

	# Aspect ratio for 3d surface calculation
	if None == spacing:
		spacing = [1.0 for d in mask.shape]

	# Force binary type
	mask = mask.astype('bool')

	# Check number of objects
	if 1 != label(mask).max():
		return(0)

	# Calculate shape descriptor
	if 2 == len(mask.shape):
		# Calculate solidity
		return(float(mask.sum()) / convex_hull_image(mask).sum())
	elif 3 == len(mask.shape):
		# Add top/bottom slices
		shape = list(mask.shape)
		shape[0] = 1
		mask = np.vstack((np.zeros(shape), mask, np.zeros(shape)))

		# Calculate sphericity
		s = marching_cubes(mask, 0.0, spacing)
		s = mesh_surface_area(s[0], s[1])
		return((np.pi * (6.0 * mask.sum())**2)**(1/3.0) / s)
	else:
		return(0)

def calc_surface(mask, spacing = None):
	""" Calculates the surface of a binary mask containing only one object. """


	# Aspect ratio for 3d surface calculation
	if None == spacing:
		spacing = [1.0 for d in mask.shape]

	# Force binary type
	mask = mask.astype('bool')

	# Check number of objects
	if 1 != label(mask).max():
		return(0)

	# Add top/bottom slices
	shape = list(mask.shape)
	shape[0] = 1
	mask = np.vstack((np.zeros(shape), mask, np.zeros(shape)))

	# Calculate sphericity
	s = marching_cubes(mask, 0.0, spacing)
	return(mesh_surface_area(s[0], s[1]))

def slice_k_d_img(img, k, slices = None):
	""" Select one 3d image from a k-d image. """

	if k > len(img.shape):
		return(img)

	idxs = [(0,) for i in range(len(img.shape) - k)]
	for i in range(k):
		idxs.append(tuple(range(img.shape[len(img.shape) - k + i])))
	img = img[np.ix_(*idxs)].reshape(img.shape[len(img.shape) - k:])

	return(img)

def get_rescaling_factor(path, **kwargs):
	""" Get rescaling factor for deconvolved. """

	# Build proper path to the deconvolution log file
	path = kwargs['basedir'] + path[0]
	path = path.split('.')
	path[len(path) - 2] = path[len(path) - 2] + '_history'
	path[len(path) - 1] = 'txt'
	path = '.'.join(path)
	
	if not os.path.exists(path):
		# Means that the image was not deconvolved
		factor = 1
	else:
		# Identify line with scaling factor
		fhistory = open(path, 'r')
		frows = fhistory.readlines()
		fhistory.close()

		# Retrieve factor
		needle = 'Stretched to Integer type'
		factor = filter(lambda x: needle in x, frows)[0]
		factor = factor.strip()
		factor = factor.split(' ')
		factor = float(factor[len(factor) - 1])

	return(factor)

def get_partial_nuclear_volume(mask, i, erosion):
	"""
	Retrieve the partial volume of a nucleus:
	 - top half (from highest intensity sum slice)
	 - central portion (through erosion)

	@param:
	 - mask <np.array[boolean]> the original mask
	 - i <np.array[uint]> the actual image
	 - erosion <float> maximum distance allowed

	@return:
	 sel_mask <np.array[boolean]> partial volume mask
	 err_log <string> error log if anything bad occurred
	"""

	# Empty error log
	err_log = ''

	# Copy original mask
	sel_mask = mask.copy()

	# Identify middle section
	mid = [i[sliceid, :, :].sum() for sliceid in range(mask.shape[0])]
	mid = mid.index(max(mid))

	# Select top half
	sel_mask[0:mid, :, :] = False

	# Select erosion on max projection
	d2d = distance_transform_edt(sel_mask.max(0))

	# Normalize distance
	d2d = d2d - d2d.min()
	if not 0 == d2d.max():
		d2d = d2d / float(d2d.max())
	else:
		# Log error
		err_log += 'Found 0 maximum distance at '
		err_log += str(self.s) + '.' + str(self.n) + '\n\n'

	# Make 2d mask
	mask2d = d2d >= erosion

	# Apply 2d erosion mask to every slice
	for i in range(sel_mask.shape[0]):
		sel_mask[i, :, :] = np.logical_and(sel_mask[i, :, :], mask2d)

	# Convert in unsigned integers
	sel_mask = sel_mask.astype('u4')

	# Output
	return((sel_mask, err_log))

def estimate_background(i, mask, seg_type):
	"""
	Estimates background median.

	@param:
	 - i <np.array> the image
	 - seg_type <string> one of the seg_type pyGPSeq constants
	"""

	if const.SEG_3D != seg_type and 2 != len(i.shape):
		bg = np.median(mk_z_projection(i, seg_type)[mask == 0])
	else:
		bg = np.median(i[mask == 0])

	return(bg)

def get_objects_xysize(L):
	""" Retrieve objects size (2/3D). """
	Larray = L.reshape([np.prod(L.shape)]).tolist()
	sizes = [t[1] for t in vt.uniquec(Larray) if t[0] != 0]
	return(sizes)

def get_objects_zsize(L):
	""" Retrieve objects size (2/3D). """
	sizes = [(L == i).astype('int').sum(0).max() for i in range(1, L.max())]
	return(sizes)

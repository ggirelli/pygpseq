# -*- coding: utf-8 -*-

'''
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
@description: methods for FISH image manipulation.
'''

# DEPENDENCIES =================================================================

import numpy as np
import os

from skimage.measure import label
from skimage.morphology import dilation

from pygpseq.fish import dot
from pygpseq.fish import nucleus

from pygpseq.tools import Binarize
from pygpseq.tools import image as imt
from pygpseq.tools.io import IOinterface, printout
from pygpseq.tools import plot

# FUNCTIONS ====================================================================

def analyze_field_of_view(sid, data, im2fov,
	dilate_factor, aspect,
	mask_dir, mask_prefix,
	plotCompartments, pole_fraction, outdir, noplot,
	labeled, compressed, istruct,
	an_type, seg_type, # Required by the Binarize class
	verbose = False):
	'''

	Args:
		sid (int): series ID.
		data (pd.DataFrame): FISH data table.
		im2fov (dict): sid-to-absPath dictionary.
	'''

	v = verbose
	msg = printout("Job '%s'..." % (im2fov[sid],), 1, v)
	subt = data.loc[np.where(data['File'] == sid)[0], :]

	# INPUT ====================================================================

	# Read image
	msg += printout("Reading image ...", 2, v)
	im = imt.read_tiff(im2fov[sid])

	# Re-slice
	msg += printout("Re-slicing ...", 2, v)
	im = imt.slice_k_d_img(im, 3)

	# Get DNA scaling factor and rescale
	sf = imt.get_rescaling_factor([im2fov[sid]],
		basedir = os.path.dirname(im2fov[sid]))
	msg += printout("Re-scaling with factor %f..." % (sf,), 2, v)
	im = (im / sf).astype('float')

	# Pick first XY plane
	if 3 == len(im.shape) and 1 == im.shape[0]: im = im[0]

	# SEGMENTATION =============================================================
	
	Segmenter = Binarize(an_type = an_type, seg_type = seg_type,
		verbose = verbose)
	
	# Check if already segmented
	already_segmented = False
	if not type(None) == type(mask_dir):
		mpath = os.path.join(mask_dir, mask_prefix + im2fov[sid])
		already_segmented = os.path.isfile(mpath)

	# Skip or binarize
	if already_segmented:
		msg += printout("Skipped binarization, using provided mask.", 3, v)
		imbin = imt.read_tiff(mpath) != 0 # Read and binarize
		thr = 0
	else:
		msg += printout("Binarizing...", 2, v)
		(imbin, thr, log) = Segmenter.run(im)
		msg += log

	# Estimate background
	dna_bg = imt.estimate_background(im, imbin, seg_type)
	msg += printout("Estimated background: %.2f a.u." % (dna_bg,), 3, v)

	# Filter based on object size
	imbin, tmp = Segmenter.filter_obj_XY_size(imbin)
	imbin, tmp = Segmenter.filter_obj_Z_size(imbin)

	# NUCLEI ===================================================================
	
	msg += printout("Retrieving nuclei...", 2, v)
	L = label(imbin)

	# Save mask ----------------------------------------------------------------
	
	# Export binary mask as TIF
	if not type(None) == type(mask_dir) and not already_segmented:
		msg += printout("Exporting mask as tif...", 4, v)
		if not os.path.isdir(mask_dir): os.mkdir(mask_dir)
		if labeled: # Export labeled mask
			plot.save_tif(mpath, L, 'uint8', compressed)
		else: # Export binary mask (min/max)
			L[np.nonzero(L)] = 255
			plot.save_tif(mpath, L, 'uint8', compressed)
			L = label(imbin)

	# Export mask as PNG
	if not noplot:
		# Create png masks output directory
		maskdir = os.path.join(outdir, 'masks/')
		if not os.path.isdir(maskdir): os.mkdir(maskdir)
		imbname = os.path.splitext(os.path.basename(im2fov[sid]))[0]

		# Save default mask
		msg += printout("Saving default binary mask...", 3, v)
		plot.export_mask_png( "%smask.%s.default.png" % (maskdir,
		imbname), imbin, im2fov[sid], "Default mask.")

		# Export dilated mask
		if 0 != dilate_factor:
			msg += printout("Saving dilated mask...", 3, v)
			plot.export_mask_png("%smask.%s.dilated%d.png" % (maskdir,
				imbname, dilate_factor), dilation(imbin, istruct), im2fov[sid],
				"Dilated mask, %d factor." % (dilate_factor,))

		# Export labeled mask
		msg += printout("Saving nuclear ID mask...", 3, v)
		plot.export_mask_png("%smask.%s.nuclei.png" % ( maskdir, imbname),
			L, im2fov[sid], 'Nuclei in "%s" [%d objects]' % (
			os.path.basename(im2fov[sid]), L.max()))

	# Store nuclei -------------------------------------------------------------
	msg, curnuclei = nucleus.build_nuclei(msg, L, dilate_factor,
		series_id = sid, thr = thr,
		dna_bg = dna_bg, sig_bg = 0,
		aspect = aspect, offset = (1, 1, 1),
		logpath = IOinterface().logpath, i = im, istruct = istruct)

	# ANALYSIS =================================================================
	
	msg += printout("Analyzing...", 2, v)

	# Assign dots to cells -----------------------------------------------------
	
	msg += printout("Assigning dots to cells...", 3, v)
	subt = dot.dots2cells(subt, curnuclei, dilate_factor)

	# Distances ----------------------------------------------------------------
	
	msg += printout("Calculating lamina distance...", 3, v)
	subt, msg = dot.calc_dot_distances(msg, subt, curnuclei, aspect)

	# Compartments -------------------------------------------------------------

	msg += printout("Annotating compartments...", 3, v)
	compdir = None

	if plotCompartments: # Create compartments output directory
		compdir = os.path.join(outdir, 'compartments/')
		if not os.path.isdir(compdir): os.mkdir(compdir)

	# Perform annotation
	subt, tvcomp, msg = nucleus.annotate_compartments(
		msg, subt, curnuclei, compdir, pole_fraction, aspect)

	# CONCLUDE =================================================================

	# Remove masks from curnuclei to free some memory
	for k in curnuclei.keys(): del curnuclei[k].mask

	# Output
	msg += printout("< Finished job.", 0, v)
	return((curnuclei, subt, tvcomp, msg))

# END ==========================================================================

################################################################################

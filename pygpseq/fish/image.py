# -*- coding: utf-8 -*-

'''
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
@description: methods for FISH image manipulation.
'''

# DEPENDENCIES =================================================================

import numpy as np
import os
import pandas as pd
import tifffile

from skimage.measure import label
from skimage.morphology import dilation

from ..tools import Binarize
from ..tools import image as imt
from ..tools import io as iot
from ..tools import plot

from .dot import dots2cells, calc_dot_distances
from .nucleus import build_nuclei, annotate_compartments

# FUNCTIONS ====================================================================

def analyze_field_of_view(ii, imfov, imdir, an_type, seg_type,
	maskdir, dilate_factor, aspect, t, main_mask_dir, main_mask_prefix,
	doCompartments, plotCompartments, pole_fraction, outdir, noplot):
	
	# Logger for logpath
	logger = iot.IOinterface()

	istruct = imt.mkIsoStruct(dilate_factor, aspect)

	idx = ii
	impath = imfov[ii]
	msg = "> Job '%s'...\n" % (impath,)
	subt_idx = np.where(t['File'] == idx)[0]

	# Read image
	msg += "   - Reading ...\n"
	im = tifffile.imread(os.path.join(imdir, impath))
	if 1 == im.shape[0]:
		im = im[0]

	# Re-slice
	msg += "    > Re-slicing ...\n"
	im = imt.autoselect_time_frame(im)
	im = imt.slice_k_d_img(im, 3)

	# Get DNA scaling factor and rescale
	sf = imt.get_rescaling_factor([impath], basedir = imdir)
	im = (im / sf).astype('float')
	msg += "    > Re-scaling with factor %f...\n" % (sf,)

	# Pick first timeframe
	if 3 == len(im.shape) and 1 == im.shape[0]:
		im = im[0]

	# Binarize image -----------------------------------------------------------
	binarization = Binarize(
		an_type = an_type, seg_type = seg_type, verbose = False)
	
	# Check if already segmented
	already_segmented = False
	if not type(None) == type(main_mask_dir):
		mpath = os.path.join(main_mask_dir, main_mask_prefix + impath)
		if os.path.isfile(mpath):
			already_segmented = True

	# Skip or binarize
	if already_segmented:
		msg += "   - Skipped binarization, using provided mask.\n"
		imbin = tifffile.imread(mpath) != 0
		thr = 0
	else:
		msg += "   - Binarizing...\n"
		(imbin, thr, log) = binarization.run(im)
		if not type(None) == type(main_mask_dir):
			if os.path.isdir(main_mask_dir):
				msg += "   >>> Exporting mask as tif...\n"
				if labeled:
					plot.save_tif(mpath, label(imbin), 'uint8', compressed)
				else:
					plot.save_tif(mpath, imbin, 'uint8', compressed)
		msg += log

	# Find nuclei --------------------------------------------------------------
	msg += "   - Retrieving nuclei...\n"

	# Estimate background
	dna_bg = imt.estimate_background(im, imbin, seg_type)
	msg += "    > Estimated background: %.2f a.u.\n" % (dna_bg,)

	# Filter object size
	imbin, tmp = binarization.filter_obj_XY_size(imbin)
	imbin, tmp = binarization.filter_obj_Z_size(imbin)

	# Save default mask
	msg += "   - Saving default binary mask...\n"
	outname = "%smask.%s.default.png" % (maskdir, os.path.splitext(impath)[0])
	if not noplot:
		plot.export_mask_png(outname, imbin, impath, "Default mask.")

	# Export dilated mask
	if not noplot and 0 != dilate_factor:
		msg += "   - Saving dilated mask...\n"
		imbin_dil = dilation(imbin, istruct)
		title = "Dilated mask, %d factor." % (dilate_factor,)
		outname = "%smask.%s.dilated%d.png" % (maskdir,
			os.path.splitext(impath)[0], dilate_factor)
		plot.export_mask_png(outname, imbin_dil, impath, title)

	# Identify nuclei
	L = label(imbin)
	seq = range(1, L.max() + 1)

	# Save mask ----------------------------------------------------------------
	msg += "   - Saving nuclear ID mask...\n"
	title = 'Nuclei in "%s" [%d objects]' % (impath, L.max())
	outpath = "%smask.%s.nuclei.png" % (maskdir, os.path.splitext(impath)[0])
	if not noplot:
		plot.export_mask_png(outpath, L, impath, title)

	# Store nuclei -------------------------------------------------------------
	msg, curnuclei = build_nuclei(msg, L, dilate_factor,
		series_id = ii, thr = thr,
		dna_bg = dna_bg, sig_bg = 0,
		aspect = aspect, offset = (1, 1, 1),
		logpath = logger.logpath, i = im)

	# Assign dots to cells -----------------------------------------------------
	msg += "   - Analysis...\n"
	msg += "    > Assigning dots to cells...\n"
	subt = dots2cells(t.loc[subt_idx, :], curnuclei, dilate_factor)

	# Distances ----------------------------------------------------------------
	msg += "    > Calculating lamina distance...\n"
	subt, msg = calc_dot_distances(msg, subt, curnuclei, aspect)

	# Compartments -------------------------------------------------------------

	if doCompartments:
		msg += "    > Annotating compartments...\n"
		compdir = None

		if plotCompartments:
			# Create compartments output directory
			compdir = os.path.join(outdir, 'compartments/')
			if not os.path.isdir(compdir):
				os.mkdir(compdir)

		# Perform annotation
		subt, tvcomp, msg = annotate_compartments(
			msg, subt, curnuclei, compdir, pole_fraction)
	else:
		msg += "    > Skipped compartments annotation.\n"

	# Clean and output ---------------------------------------------------------

	# Remove masks from curnuclei
	for k in curnuclei.keys():
		del curnuclei[k].mask

	# Output
	msg += "< Finished job.\n"
	#print(msg)
	return((curnuclei, subt, subt_idx, tvcomp, msg))

# END ==========================================================================

################################################################################

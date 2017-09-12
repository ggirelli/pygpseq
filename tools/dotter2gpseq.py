#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 1.1.0
# Date: 20170718
# Project: GPSeq
# Description: Calculate radial position of dots in cells
# 
# Changelog:
#  v1.1.0 - 20170830: added G1 cells selection.
#  v1.0.0 - 20170718: first implementation.
#  
# Todo:
#  - Parallelize when possible.
#  - Allow nucleus dilation or distance from lamina for out-of-nucleus dots.
# 
# ------------------------------------------------------------------------------



# DEPENDENCIES =================================================================

import argparse
import math
import numpy as np
import os
import pandas as pd
from scipy.ndimage.morphology import distance_transform_edt
import skimage.io as io
from skimage.measure import label

import pygpseq as gp
from pygpseq.wraps.nucleus import Nucleus
from pygpseq.tools import image as imt
from pygpseq.tools import io as iot
from pygpseq.tools import stat as stt

# PARAMETERS ===================================================================

# Add script description
parser = argparse.ArgumentParser(description = '''
Calculate radial position of dots in cells. Use the -s option to trigger nuclei
recognition and G1 selection: a column will be added flagging dots belonging to
cells expected to be in G1. The G1 selection is actually a selection of the most
represented cell sub-population based on flatten area and integral of DNA stain
intensity. In other words, it will selected the most represented cell cycle
phase in your cell population (generally, G1). Please, note that nuclei
recognition and G1 selection are time consuming steps and require a large number
of nuclei to work properly (i.e., more than 300). Images are expected to follow
DOTTER filename notation: 'channel_series.tif'.
''')

# Add mandatory arguments
parser.add_argument('dotCoords', type = str, nargs = 1,
	help = 'Dot coordinates table generated by DOTTER.')
parser.add_argument('imgFolder', type = str, nargs = 1,
	help = 'Folder with tiff images.')
parser.add_argument('output', type = str, nargs = 1,
	help = 'Output tsv table.')

# Optional parameters
parser.add_argument('-a', '--aspect', type = float, nargs = 3,
	help = """Physical size of Z, Y and X voxel sides.
	Default: 300.0 130.0 130.0""",
	metavar = ('Z', 'Y', 'X'), default = [300., 130., 130.])
parser.add_argument('-d', '--delim', type = float, nargs = 1,
	help = """Input table delimiter. Default: ','""", default = [','])

# Add flags
parser.add_argument('-s',
	action = 'store_const', dest = 'sel',
	const = True, default = False,
	help = 'Perform nuclei recognition and G1 selection (time consuming).')

# Parse arguments
args = parser.parse_args()

# Assign to in-script variables
dot_table_name = args.dotCoords[0]
imdir = args.imgFolder[0]
aspect = args.aspect
(az, ay, ax) = aspect
outpath = args.output[0]
doSel = args.sel
delim = args.delim[0]

# Params
seg_type = gp.const.SEG_3D
an_type = gp.const.AN_3D

# FUNCTIONS ====================================================================

def in_box(coords, box):
	''''''
	c = True
	for dim in range(len(coords)):
		c = c and coords[dim] >= box[dim][0] and coords[dim] <= box[dim][1]

	return c

def in_nucleus(n, s, coords):
	''''''
	if not n.s == s:
		return False
	
	return in_box(coords, n.box)


# RUN ==========================================================================

# Read table
t = pd.read_csv(dot_table_name, delim)
t['cell_ID'] = np.zeros(len(t.index))
t['lamin_dist'] = np.zeros(len(t.index))
t['lamin_dist_norm'] = np.zeros(len(t.index))
t['centr_dist'] = np.zeros(len(t.index))
t['centr_dist_norm'] = np.zeros(len(t.index))

# Extract FoV number
t['File'] = [int(f.split('/')[-1].split('.')[0]) for f in t['File']]

# Identify tiff images
flist = []
for (dirpath, dirnames, filenames) in os.walk(imdir):
    flist.extend(filenames)
    break
imlist = [f for f in flist if 'tif' in f]

# Assign field of views to images
imfov = {}
for i in set(t['File']):
	imfov[i] = [im for im in imlist if "%03d" % (i,) in im][0]

# Nuclei container
nuclei = []

# Logger for logpath
logger = iot.IOinterface()

# Cycle through
for ii in range(len(imfov.keys())):
	(idx, impath) = list(imfov.items())[ii]
	print("  · '%s'..." % (impath,))
	subt_idx = np.where(t['File'] == idx)[0]

	# Read image
	print("   - Reading ...")
	im = io.imread(os.path.join(imdir, impath))[0]

	# Re-slice
	print("    > Re-slicing ...")
	im = imt.autoselect_time_frame(im)
	im = imt.slice_k_d_img(im, 3)

	# Get DNA scaling factor and rescale
	sf = imt.get_rescaling_factor([impath], basedir = imdir)
	im = (im / sf).astype('float')
	print("    > Re-scaling with factor %f..." % (sf,))

	# Pick first timeframe
	if 3 == len(im.shape) and 1 == im.shape[0]:
		im = im[0]

	# Binarize images
	print("   - Binarizing...")
	binarization = gp.tools.binarize.Binarize(
		an_type=an_type,
		seg_type=seg_type
	)
	(imbin, thr, log) = binarization.run(im)

	# Find nuclei
	if doSel:
		print("   - Retrieving nuclei...")

		# Estimate background
		dna_bg = imt.estimate_background(im, imbin, seg_type)
		print("    > Estimated background: %.2f a.u." % (dna_bg,))

		# Filter object size
		imbin, tmp = binarization.filter_obj_XY_size(imbin)
		imbin, tmp = binarization.filter_obj_Z_size(imbin)

		# Identify and store nuclei
		L = label(imbin)
		seq = range(1, L.max() + 1)
		kwargs = {
			'series_id' : ii, 'thr' : thr,
			'dna_bg' : dna_bg, 'sig_bg' : 0,
			'aspect' : aspect, 'offset' : (1, 1, 1),
			'logpath' : logger.logpath, 'i' : im
		}

		print("    > Saving %d nuclei..." % (L.max(),))
		nuclei.extend([Nucleus(n = n, mask = L == n, **kwargs) for n in seq])

	# Calculate distance from lamina
	print("   - Analysis...")
	print("    > Calculating lamina distance...")

	# Add empty top/bottom slides
	imbin_tops = np.zeros((imbin.shape[0]+2, imbin.shape[1], imbin.shape[2]))
	imbin_tops[1:(imbin.shape[0]+1),:,:] = imbin

	# Extract cell ID
	print("    > Assigning dots to cells...")
	L = label(imbin_tops)[1:(imbin.shape[0]+1),:,:]
	t.loc[subt_idx, 'cell_ID'] = L[t['z'], t['x'], t['y']][subt_idx]

	# Calculate distance and store it
	print("    > Calculating distances...")
	D = distance_transform_edt(imbin_tops, aspect)[1:(imbin.shape[0]+1),:,:]
	t.loc[subt_idx, 'lamin_dist'] = D[t['z'], t['x'], t['y']][subt_idx]

	# Retrieve max lamin dist per cell
	cell_max_lamin_dist = {}
	for cellID in set(t.loc[subt_idx, 'cell_ID'].tolist()):
		cell_max_lamin_dist[cellID] = np.max(
			D[np.where(L == cellID)])

	# Normalize lamin_dist
	fnorm = [cell_max_lamin_dist[cellID]
		for cellID in t.loc[subt_idx, 'cell_ID'].tolist()]
	t.loc[subt_idx, 'lamin_dist_norm'] = t.loc[subt_idx, 'lamin_dist'] / fnorm

	# Normalized centr_dist
	t.loc[subt_idx, 'centr_dist_norm'] = 1 - t.loc[subt_idx, 'lamin_dist_norm']

	# Calculate centr_dist
	t.loc[subt_idx, 'centr_dist'] = t.loc[subt_idx, 'centr_dist_norm'] * fnorm

# Identify G1 cells
if doSel:
	print("  - Flagging G1 cells...")

	# Retrieve nuclei summaries
	print('   > Retrieving nuclear summary...')
	summary = np.zeros(len(nuclei),
		dtype = gp.const.DTYPE_NUCLEAR_SUMMARY)
	for i in range(len(nuclei)):
		summary[i] = nuclei[i].get_summary()

	# Filter nuclei
	print('   > Filtering nuclei based on flatten size and intensity...')
	cond_name = 'none'
	sigma = .1
	nsf = (gp.const.NSEL_FLAT_SIZE, gp.const.NSEL_SUMI)
	out_dir = '.'

	# Filter features
	sel_data = {}
	plot_counter = 1
	for nsfi in nsf:
		# Identify Nuclear Selection Feature
		nsf_field = gp.const.NSEL_FIELDS[nsfi]
		nsf_name = gp.const.NSEL_NAMES[nsfi]
		print('   >> Filtering %s...' % (nsf_name,))

		# Start building output
		d = {'data' : summary[nsf_field]}

		# Calculate density
		d['density'] = stt.calc_density(d['data'], sigma = sigma)

		# Identify range
		args = [d['density']['x'], d['density']['y']]
		d['fwhm_range'] = stt.get_fwhm(*args)

		# Plot
		sel_data[nsf_field] = d

	# Select based on range
	f = lambda x, r: x >= r[0] and x <= r[1]
	for nsfi in nsf:
		nsf_field = gp.const.NSEL_FIELDS[nsfi]
		nsf_name = gp.const.NSEL_NAMES[nsfi]
		print("   > Selecting range for %s ..." % (nsf_name,))

		# Identify nuclei in the FWHM range
		nsf_data = sel_data[nsf_field]
		nsf_data['sel'] = [f(i, nsf_data['fwhm_range'])
			for i in nsf_data['data']]
		sel_data[nsf_field] = nsf_data
	
	# Select those in every FWHM range
	print("   > Applying selection criteria")
	nsfields = [gp.const.NSEL_FIELDS[nsfi] for nsfi in nsf]
	selected = [sel_data[f]['sel'] for f in nsfields]
	g = lambda i: all([sel[i] for sel in selected])
	selected = [i for i in range(len(selected[0])) if g(i)]
	sub_data = np.array(summary[selected])

	# Identify selected nuclei objects
	sel_nucl = []
	for n in nuclei:
		if n.n in sub_data['n'][np.where(n.s == sub_data['s'])[0]]:
			sel_nucl.append(n)

	# Check which dots are in which nucleus and update flag
	print("   > Matching DOTTER cells with GPSeq cells...")
	t['G1'] = np.zeros((t.shape[0],))
	for ti in t.index:
		for ni in range(len(sel_nucl)):
			n = sel_nucl[ni]
			if in_nucleus(n, int(t.ix[ti, 0]-1), tuple(t.ix[ti, [5, 3, 4]])):
				t.ix[ti, 'G1'] = 1
				break

# Write output
t.to_csv(outpath, sep = '\t', index = False)

# END ==========================================================================

################################################################################

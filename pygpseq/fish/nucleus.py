# -*- coding: utf-8 -*-

'''
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
@description: methods for nuclear channels manipulation.
'''

# DEPENDENCIES =================================================================

import matplotlib
import matplotlib.pyplot as plt

import math
import numpy as np
import os
import pandas as pd
from scipy.ndimage.measurements import center_of_mass
from skimage import draw
from skimage.morphology import dilation

from .. import const
from ..anim import Nucleus
from ..tools import image as imt
from ..tools import plot
from ..tools import stat as stt

# FUNCTIONS ====================================================================

def annotate_compartments(msg, t, nuclei, outdir, pole_fraction):
	'''
	Add compartment status to dots table (by DOTTER).
	For each nucleus: the major three axes are identified, the nucleus is
	centered and rotated. Then, the dots are also centered and rotated,
	and assigned to different compartments based on the fitted ellipsoid.
	Information on the goodness of ellipsoid fit is added to the main log and
	can be extracted by grepping lines starting with "   >>>> GoF_ellipse:".

	Args:
	   t (pd.DataFrame): DOTTER output table.
	   msg (string): log message, to be continued.

	Returns:

	'''

	# Temporarily remove dots outside cells
	nan_cond = np.isnan(t.loc[:, 'cell_ID'])
	vcomp_table = pd.DataFrame()

	subt = t.loc[np.logical_not(nan_cond), :].copy()
	if 0 == subt.shape[0]:
		print("!WARNING! All dots in FoV#%d are outside cells." % (
			t['File'].values[0],))
		return((t, vcomp_table, msg))

	fid = subt['File'].values[0]
	
	# Create empty table to host compartment volume data
	vcomp_table = pd.DataFrame(index = range(1, int(subt['cell_ID'].max()) + 1))
	vcomp_table['File'] = fid
	vcomp_table['cell_ID'] = range(1, int(subt['cell_ID'].max()) + 1)
	vcomp_table['center_bot'] = np.nan
	vcomp_table['center_top'] = np.nan
	vcomp_table['poles'] = np.nan
	vcomp_table['ndots_center_bot'] = np.nan
	vcomp_table['ndots_center_top'] = np.nan
	vcomp_table['ndots_poles'] = np.nan

	for cid in range(int(subt['cell_ID'].max()) + 1):
		if cid in nuclei.keys():
			msg += "    >>> Working on cell #%d...\n" % (cid,)
			cell_cond = cid == subt['cell_ID']

			# Extract dots coordinates -----------------------------------------
			dot_coords = np.vstack([
				subt.loc[cell_cond, 'z'] - nuclei[cid].box_origin[0],
				subt.loc[cell_cond, 'x'] - nuclei[cid].box_origin[1],
				subt.loc[cell_cond, 'y'] - nuclei[cid].box_origin[2]
			])

			# Center coordinates -----------------------------------------------
			x, y, z, xd, yd, zd = stt.centered_coords_3d(
				nuclei[cid].mask, dot_coords)
			coords = np.vstack([x, y, z])
			dot_coords = np.vstack([xd, yd, zd])

			# Rotate data ------------------------------------------------------
			
			# First axis
			xv, yv, zv = stt.extract_3ev(coords)
			theta1 = stt.calc_theta(xv[0], yv[0])
			xt, yt, zt = stt.rotate3d(coords, theta1, 2)
			tcoords = np.vstack([xt, yt, zt])

			# Third axis
			xv, yv, zv = stt.extract_3ev(tcoords)
			theta3 = stt.calc_theta(xv[2], zv[2])
			if np.abs(theta3) > np.pi / 2.:
				if theta3 > 0:
					theta3 = -np.abs(theta3 - np.pi / 2.)
				else:
					theta3 = np.abs(theta3 + np.pi / 2.)
			else:
				theta3 = -np.abs(theta3 + np.pi / 2.)
			xt, yt, zt = stt.rotate3d(tcoords, theta3, 1)
			tcoords = np.vstack([xt, yt, zt])

			# Second axis
			xv, yv, zv = stt.extract_3ev(tcoords)
			theta2 = stt.calc_theta(yv[1], zv[1])
			xt, yt, zt = stt.rotate3d(tcoords, theta2, 0)
			tcoords = np.vstack([xt, yt, zt])

			# Fit ellipsoid ----------------------------------------------------

			# Round up rotated coordinates
			trcoords = tcoords.astype('i')

			# Convert to rotated image
			icoords = np.transpose(trcoords) + abs(trcoords.min(1))
			trbin = np.zeros((icoords.max(0) + 1).tolist()[::-1])
			trbin[icoords[:, 2], icoords[:, 1], icoords[:, 0]] = 1

			# Calculate axes size
			zax_size, xax_size, yax_size = trbin.shape

			el = draw.ellipsoid(zax_size / 2., xax_size / 2., yax_size / 2.)
			el = el[2:(zax_size + 2), 1:(xax_size + 1), 1:(yax_size + 1)]

			# Calculate intersection with fitting ellipsoid
			inter_size = np.logical_and(trbin, el).sum()

			# Log intersection
			comments = []
			comments.append("%s%%%s [%s.%s]." % (
				round(inter_size / float(trbin.sum()) * 100, 2,),
				" of the nucleus is in the ellipsoid", fid, cid,))
			comments.append("%s%%%s [%s.%s]." % (
				round(inter_size / float(el.sum()) * 100, 2,),
				" of the ellipsoid is in the nucleus", fid, cid,))
			msg += "".join(["   >>>> GoF_ellipse: %s\n" % (s,)
				for s in comments])

			# Identify ellipsoid foci
			b = xax_size / 2.
			a = yax_size / 2.
			c = np.sqrt(a**2 - b**2)
			#ecc = c / a

			# Rotate dots ------------------------------------------------------

			dot_coords_t = np.vstack(stt.rotate3d(dot_coords, theta1, 2))
			dot_coords_t = np.vstack(stt.rotate3d(dot_coords_t, theta2, 0))
			dot_coords_t = np.vstack(stt.rotate3d(dot_coords_t, theta3, 1))

			# Assign compartments ----------------------------------------------
			# Compartment code:
			# 0 = center-top
			# 1 = center-bottom
			# 2 = pole
			cf = 1 - 2 * pole_fraction
			status = np.zeros(dot_coords.shape[1])
			status[dot_coords_t[2] < 0] = 1
			status[dot_coords_t[0] > cf * a] = 2
			status[dot_coords_t[0] < -(cf * a)] = 2
			subt.loc[cell_cond, 'compartment'] = status

			# Calculate compartment volume -------------------------------------

			# Round up coordinates
			xt = xt.astype('i')
			zt = zt.astype('i')

			# Count voxels in compartments
			vpole = sum(xt > c) + sum(xt < -c)
			centr_cond = np.logical_and(xt < c, xt > -c)
			vctop = np.logical_and(centr_cond, zt >= 0).sum()
			vcbot = np.logical_and(centr_cond, zt < 0).sum()

			vcomp_table.loc[cid, 'center_top'] = vctop
			vcomp_table.loc[cid, 'center_bot'] = vcbot
			vcomp_table.loc[cid, 'poles'] = vpole
			vcomp_table.loc[cid, 'ndots_center_top'] = (status == 0).sum()
			vcomp_table.loc[cid, 'ndots_center_bot'] = (status == 1).sum()
			vcomp_table.loc[cid, 'ndots_poles'] = (status == 2).sum()

			# Assign volume information
			volume = np.zeros(dot_coords.shape[1])
			volume[:] = vctop
			volume[dot_coords_t[2] < 0] = vcbot
			volume[dot_coords_t[1] > 2 * c - a] = vpole
			volume[dot_coords_t[1] < -(2 * c - a)] = vpole
			subt.loc[cell_cond, 'compartment_volume'] = volume

			# Generate compartment plot with dots ------------------------------
			
			if not type(None) == type(outdir):
				outpng = open(os.path.join(outdir,
					"%s.%s.png" % (fid, cid,)), "wb")
				plt.close("all")
				plot.ortho_3d(tcoords, dot_coords = dot_coords_t, c = a * 0.6)
				plt.suptitle("\n".join(comments))
				plt.savefig(outpng, format = "png")
				plt.close("all")
				outpng.close()

			t.loc[np.logical_not(nan_cond), :] = subt

	return((t, vcomp_table, msg))

def build_nuclei(msg, L, dilate_factor, series_id, thr, dna_bg, sig_bg,
	aspect, offset, logpath, i):
	'''
	Build nuclei objects
	
	Args:
	  msg (string): log message, to be continued.
	  L (np.ndarray): labeled mask.
	  dilate_factor (int): dilation factor.
	  series_id (int): series ID.
	  thr (float): global threshold value.
	  dna_bg (float): DNA channel background.
	  sig_bg (float): signal channel background.
	  aspect (tuple): Z,Y,X voxel sides in real units.
	  offset (tuple): tuple with pixel offset for bounding box.
	  logpath (string): path to log file.
	  i (np.array): image.
	
	Returns:
	  (string, list): log message and list of Nucleus objects.
	'''
	
	# Prepare input for Nucleus class
	kwargs = {
		'series_id' : series_id, 'thr' : thr,
		'dna_bg' : dna_bg, 'sig_bg' : sig_bg,
		'aspect' : aspect, 'offset' : offset,
		'logpath' : logpath, 'i' : i
	}

	istruct = imt.mkIsoStruct(dilate_factor, aspect)

	# Default nuclear ID list and empty dictionary
	seq = range(1, L.max() + 1)
	curnuclei = {}

	# Log operation
	if 0 != dilate_factor:
		msg += "   - Saving %d nuclei with dilation [%d]...\n" % (
			L.max(), dilate_factor)
	else:
		msg += "   - Saving %d nuclei...\n" % (L.max(),)

	# Iterate through nuclei
	for n in seq:
		# Make nucleus
		if 0 != dilate_factor:
			# With dilated mask
			mask = dilation(L == n, istruct)
			nucleus = Nucleus(n = n, mask = mask, **kwargs)
		else:
			mask = L == n
			nucleus = Nucleus(n = n, mask = mask, **kwargs)

		# Apply box
		msg += "    > Applying nuclear box [%d]...\n" % (n,)
		mask = imt.apply_box(mask, nucleus.box)

		# Store nucleus
		nucleus.mask = mask
		nucleus.box_origin = np.array([c[0] + 1 for c in nucleus.box])
		nucleus.box_sides = np.array([np.diff(c) for c in nucleus.box])
		nucleus.box_mass_center = center_of_mass(mask)
		nucleus.dilate_factor = dilate_factor
		curnuclei[n] = nucleus

	return((msg, curnuclei))

def flag_G1_cells(t, nuclei, outdir, dilate_factor, dot_file_name):
	'''
	Assign a binary flag identifying the predominant cell population
	based on flatten size and intensity sum
	
	Args:
	  t (pd.DataFrame): DOTTER output table.
	  nuclei (list(gp.Nucleus)): identified nuclei.
	  outdir (string): path to output folder.
	  dilate_factor (int): number of dilation operations.
	  dot_file_name (string): output file name.
	
	Returns:
	  pd.DataFrame:.
	'''
	
	print("> Flagging G1 cells...")

	# Retrieve nuclei summaries ------------------------------------------------
	print('   > Retrieving nuclear summary...')
	summary = np.zeros(len(nuclei),
		dtype = const.DTYPE_NUCLEAR_SUMMARY)
	for i in range(len(nuclei)):
		summary[i] = nuclei[i].get_summary()

	# Filter nuclei ------------------------------------------------------------
	print('   > Filtering nuclei based on flatten size and intensity...')
	cond_name = 'none'
	sigma = .1
	nsf = (const.NSEL_FLAT_SIZE, const.NSEL_SUMI)
	out_dir = '.'

	# Filter features
	sel_data = {}
	ranges = {}
	plot_counter = 1
	for nsfi in nsf:
		# Identify Nuclear Selection Feature
		nsf_field = const.NSEL_FIELDS[nsfi]
		nsf_name = const.NSEL_NAMES[nsfi]
		print('   >> Filtering %s...' % (nsf_name,))

		# Start building output
		d = {'data' : summary[nsf_field]}

		# Calculate density
		d['density'] = stt.calc_density(d['data'], sigma = sigma)

		# Identify range
		args = [d['density']['x'], d['density']['y']]
		d['fwhm_range'] = stt.get_fwhm(*args)
		ranges[nsf_name] = d['fwhm_range']

		# Plot
		sel_data[nsf_field] = d

	# Select based on range
	f = lambda x, r: x >= r[0] and x <= r[1]
	for nsfi in nsf:
		nsf_field = const.NSEL_FIELDS[nsfi]
		nsf_name = const.NSEL_NAMES[nsfi]
		print("   > Selecting range for %s ..." % (nsf_name,))

		# Identify nuclei in the FWHM range
		nsf_data = sel_data[nsf_field]
		nsf_data['sel'] = [f(i, nsf_data['fwhm_range'])
			for i in nsf_data['data']]
		sel_data[nsf_field] = nsf_data

	# Select those in every FWHM range
	print("   > Applying selection criteria")
	nsfields = [const.NSEL_FIELDS[nsfi] for nsfi in nsf]
	selected = [sel_data[f]['sel'] for f in nsfields]
	g = lambda i: all([sel[i] for sel in selected])
	selected = [i for i in range(len(selected[0])) if g(i)]
	sub_data = np.array(summary[selected])

	# Identify selected nuclei objects
	sel_nuclei_labels = ["_%d.%d_" % (n, s)
		for (n, s) in sub_data[['s', 'n']]]
	sel_nucl = [n for n in nuclei
		if "_%d.%d_" % (n.s, n.n) in sel_nuclei_labels]

	# Check which dots are in which nucleus and update flag --------------------
	print("   > Matching DOTTER cells with GPSeq cells...")
	t['G1'] = 0
	t.loc[np.where(np.isnan(t['cell_ID']))[0], 'G1'] = np.nan
	t['universalID'] =  ["_%s.%s_" % x for x in zip(
		t['File'].values, t['cell_ID'].values
	)]
	g1ids = [i for i in range(t.shape[0])
		if t.loc[i, 'universalID'] in sel_nuclei_labels]
	t.loc[g1ids, 'G1'] = 1
	t = t.drop('universalID', 1)

	# Add G1 status to summary -------------------------------------------------
	summary = pd.DataFrame(summary)
	summary['G1'] = np.zeros((summary.shape[0],))
	summary['universalID'] =  ["_%s.%s_" % x
		for x in zip(summary['s'].values, summary['n'].astype("f").values)]
	g1ids = [i for i in range(summary.shape[0])
		if summary.loc[i, 'universalID'] in sel_nuclei_labels]
	summary.loc[g1ids, 'G1'] = 1
	summary = summary.drop('universalID', 1)

	# Estimate radius ----------------------------------------------------------
	summary['sphere_radius'] = summary['size'].values * 3 / (4 * math.pi)
	summary['sphere_radius'] = (summary['sphere_radius'])**(1/3.)

	# Export -------------------------------------------------------------------

	# Export feature ranges
	s = ""
	for (k, v) in ranges.items():
		s += "%s\t%f\t%f\n" % (k, v[0], v[1])
	f = open("%s/feature_ranges.txt" % (outdir,), "w+")
	f.write(s)
	f.close()

	# Export summary
	outname = "%s/nuclei.out.dilate%d.%s" % (
		outdir, dilate_factor, dot_file_name)
	summary.to_csv(outname, sep = '\t', index = False)

	# Output -------------------------------------------------------------------
	print("> Flagged G1 cells...")
	return(t)

# END ==========================================================================

################################################################################

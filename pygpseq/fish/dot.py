# -*- coding: utf-8 -*-

'''
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
@description: methods for FISH dots manipulation.
'''

# DEPENDENCIES =================================================================

import numpy as np
import pandas as pd

from scipy.ndimage.morphology import distance_transform_edt

from ..tools import image as imt

# FUNCTIONS ====================================================================

def add_allele(data):
	'''
	Add allele labels to DOTTER-based table with GPSeq-like centrality.
	
	Labels:
	  NaN : dot outside of cells.
	  -1  : more than 2 dots per cell.
	  0   : less than 2 dots per cell.
	  1   : central dot.
	  2   : peripheral dot.
	
	Args:
	  data (pd.DataFrame): DOTTER-based table with GPSeq-like centrality.
						   Required columns:
							  cell_ID, lamin_dist_norm, File, Channel
		Returns:
	  pd.DataFrame: input data table with added Allele column (label).
	'''

	# Initial checks -----------------------------------------------------------

	# Check that the format corresponds
	if not type(data) == type(pd.DataFrame()):
		print("Input should be a DataFrame from the pandas library.")
		return(data)

	# Check that required columns are present
	req_cols = ['cell_ID', 'lamin_dist_norm', 'File', 'Channel']
	check_cols = [True for c in req_cols if c in data.columns.tolist()]
	if not all(check_cols):
		miss_cols = [req_cols[i]
			for i in range(len(req_cols)) if not check_cols[i]]
		print("Some required columns are missing: %s" % (", ".join(miss_cols),))
		return(data)

	# Universal index and dots in cells ----------------------------------------

	# Default value of np.nan for dots outside of nuclei
	data['Allele'] = np.nan

	# Identify dots within cells
	validIdx = np.where(np.logical_not(np.isnan(data['cell_ID'])))[0]
	subt = data.loc[validIdx, :]

	# Assemble universal index
	subt['universalID'] =  ["%s_%s_%s" % t for t in zip(
		subt['File'].values, subt['Channel'].values, subt['cell_ID'].values
	)]

	# Count dots per universalID
	uID,  uCount = np.unique(subt.loc[validIdx, 'universalID'],
		return_index = False, return_counts = True)
	IDmap = zip(subt.loc[validIdx, 'universalID'],
		[dict(zip(uID, uCount))[ID]
		for ID in subt.loc[validIdx, 'universalID']])
	IDmap = np.array(list(IDmap))
	
	# Stop if no dots are inside a cell
	if 0 == sum(IDmap.shape):
		return(data)

	# Fill Allele column -------------------------------------------------------

	# -1 if more than 2 dots
	cond = IDmap[:,1].astype('i') > 2
	if 0 != sum(cond):
		subt.loc[validIdx[cond], 'Allele'] = -1

	#  0 if less than 2 dots
	cond = IDmap[:,1].astype('i') == 1
	if 0 != sum(cond):
		subt.loc[validIdx[cond], 'Allele'] = 0

	# Iterate over 2-dots cases
	cond = IDmap[:,1].astype('i') == 2
	if 0 != sum(cond):
		uID = np.unique(IDmap[cond, 0]).tolist()
		for ID in uID:
			dotPair = subt.loc[subt['universalID'] == ID, :]
			ldn = dotPair['lamin_dist_norm'].tolist()
			if ldn[0] == ldn[1]:
				# Same centrality
				subt.loc[dotPair.index[0], 'Allele'] = 1 # Central
				subt.loc[dotPair.index[1], 'Allele'] = 2 # Peripheral
			else: # Different centrality
				# Peripheral
				subt.loc[dotPair['lamin_dist_norm'].argmin(), 'Allele'] = 2
				# Central
				subt.loc[dotPair['lamin_dist_norm'].argmax(), 'Allele'] = 1

	# Output -------------------------------------------------------------------
	data.loc[validIdx, 'Allele'] = subt['Allele']
	return(data)

def calc_dot_distances(msg, t, nuclei, aspect):
	'''
	Calculate distance of dots from lamina and central area
	
	Args:
	  msg (string): log message, to be continued.
	  t (pd.DataFrame): DOTTER output table.
	  nuclei (list(gp.Nucleus)): identified nuclei.
	  aspect (tuple): Z,Y,X voxel sides in real units.
	
	Returns:
	  pd.DataFrame: updated dotter table.
	  str: message log..
	'''

	# Skip if no cells are present
	if ( np.all(np.isnan(t['cell_ID'])) ):
		return((t, msg))

	# Calculate distances ------------------------------------------------------
	for cid in range(int(np.nanmax(t['cell_ID'])) + 1):
		if cid in nuclei.keys():
				msg += "    >>> Working on cell #%d...\n" % (cid,)
				cell_cond = cid == t['cell_ID']

				# Distance from lamina and center
				laminD = distance_transform_edt(nuclei[cid].mask, aspect)
				centrD = distance_transform_edt(laminD != laminD.max(), aspect)

				t.loc[cell_cond, 'lamin_dist'] = laminD[
					t.loc[cell_cond, 'z'] - nuclei[cid].box_origin[0],
					t.loc[cell_cond, 'x'] - nuclei[cid].box_origin[1],
					t.loc[cell_cond, 'y'] - nuclei[cid].box_origin[2]
				]

				t.loc[cell_cond, 'centr_dist'] = centrD[
					t.loc[cell_cond, 'z'] - nuclei[cid].box_origin[0],
					t.loc[cell_cond, 'x'] - nuclei[cid].box_origin[1],
					t.loc[cell_cond, 'y'] - nuclei[cid].box_origin[2]
				]

	# Normalize distances ------------------------------------------------------

	# Max distance for each dot
	fnorm = t.loc[:, 'lamin_dist'] + t.loc[:, 'centr_dist']
	t.loc[:, 'centr_dist_norm'] = t.loc[:, 'centr_dist'] / fnorm
	t.loc[:, 'lamin_dist_norm'] = t.loc[:, 'lamin_dist'] / fnorm

	# Output
	return((t, msg))

def dots2cells(t, nuclei, dilate_factor):
	'''
	Assign dots to cells
	
	Args:
	  t (pd.DataFrame): DOTTER output subset.
	  nuclei (list(gp.Nucleus)): identified nuclei.
	  dilate_factor (int): number of dilation operations.
	
	Returns:
	  pd.DataFrame: updated DOTTER output.
	'''
	
	for idx in t.index:
		coords = ( t.loc[idx, 'z'], t.loc[idx, 'x'], t.loc[idx, 'y'] )
		for (nid, n) in nuclei.items():
			if imt.in_mask(coords - n.box_origin, n.mask):
				t.loc[idx, 'cell_ID'] = nid
				continue

	# Output
	return(t)

# END ==========================================================================

################################################################################

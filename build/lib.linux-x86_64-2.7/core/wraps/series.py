# -*- coding: utf-8 -*-

import math
import os

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import numpy as np
from scipy import ndimage as ndi
from skimage.filters import threshold_otsu
from skimage.io import imread, imsave
from skimage.measure import label

from .. import const
from ..tools import io as iot, image as imt, plot, stat as stt, string as st
from ..tools import vector as vt

from .nucleus import Nucleus

class Series(iot.ioInterface):
	"""docstring for series"""

	# Class-bound attributes
	__version__ = const.VERSION

	# Series id
	n = 0

	# Wrapped nuclei
	nuclei = []

	# Series name
	name = ''

	# Condition path (series folder)
	basedir = '.'

	# Series files info
	filist = []

	def __init__(self, ds, condition = None, **kwargs):
		"""
		@param:
		 - ds <dict> series information list
		 - condition <pyGPSeq.wrap.condition> condition wrapper (opt)
		"""

		# If required, inherit from `condition` wrap
		if None != condition:
			logpath = condition.logpath
			super(Series, self).__init__(path = logpath, append = True)
			self.basedir = condition.path
		else:
			super(Series, self).__init__()
		
		# Save input parameters
		self.name = ds[0]
		self.filist = ds[1]
		self.n = ds[2]

	def get_c(self):
		""" Returns number of channels in the series. """
		return(len(self.filist))

	def get_channel_names(self, channel_field = None):
		""" Returns the names of the channels in the series. """
		if None == channel_field:
			channel_field = const.REG_CHANNEL_NAME
		return([c[channel_field] for c in self.filist.values()])

	def find_channel(self, channel_names):
		""" Return the first channel to correspond to channel_names. """

		# Fix the param type
		if type(str()) == type(channel_names):
			channel_names = [channel_names]

		# Cycle through the available channels
		for cname in channel_names:

			# Identify the requested channel
			idx = self.find_channel_id(cname)

			# Return the channel
			if -1 != idx:
				return([i for i in self.filist.items()][idx])

		# Return empty dictionary if no matching channel is found
		return({})

	def find_channel_id(self, channel_name):
		""" Return the id of the channel file with the specified name. """

		# Retrieve available channel names
		names = self.get_channel_names()

		if 0 != names.count(channel_name):
			# Return matching channel id
			return(names.index(channel_name))
		else:
			# Return -1 if no matching channel is found
			return(-1)

	def find_nuclei(self, **kwargs):
		"""
		Segment current series.

		@param:
		 - dna_names <tuple[string]> dna channel names
		 - cond_name <string> condition wrapper name
		 - seg_type <pyGPSeq.const> segmentation type
		 - rm_z_tips <bool> remove nuclei touching the tips of the stack
		 - radius_interval <tuple[float]> allowed nuclear radius interval
		 - offset <tuple[int]> dimensions box/square offset
		 - aspect <tuple[float]> pixel/voxel dimension proportion
		"""

		# Set output suffix
		if not 'suffix' in kwargs.keys():
			suffix = ''
		else:
			suffix = st.add_leading_dot(kwargs['suffix'])

		# Check plotting
		if not 'plotting' in kwargs.keys():
			kwargs['plotting'] = True

		log = ""
		log += self.printout('Current series: "' + self.name + '"...', 1)

		# Read images
		kwargs, alog = self.adjust_options(read_only_dna = False, **kwargs)
		log += alog

		# Extract from kwargs
		seg_type = kwargs['seg_type']
		dna_ch = kwargs['dna_ch']
		sig_ch = kwargs['sig_ch']

		# Make new channel copy
		i = dna_ch.copy()

		# Make Z-projection
		if const.SEG_3D != seg_type and 2 != len(i.shape):
			msg = 'Generating Z-projection [' + seg_type + ']...'
			log += self.printout(msg, 2)
			i = imt.mk_z_projection(i, seg_type)

		# Binarize
		thr = threshold_otsu(i)
		log += self.printout('Thresholding image, global thr: ' + str(thr), 2)
		mask = imt.binarize(i, thr)

		# Perform adaptive threshold if necessary
		if 1 < kwargs['adp_thr']:
			msg = 'Applying adaptive threshold to neighbourhood: '
			msg += str(kwargs['adp_thr'])
			log += self.printout(msg, 2)
			lmask = imt.threshold_adaptive(i, kwargs['adp_thr'])

			# Combine adaptive and global masks
			log += self.printout('Combining global and local thresholds...', 2)
			mask = np.logical_and(mask, lmask)

		# Remove objects touching the borders
		log += self.printout('Removing objects touching the image border...', 2)
		mask = imt.clear_borders(mask, kwargs['rm_z_tips'])

		# Fill holes
		log += self.printout('Filling holes...', 2)
		mask = ndi.binary_fill_holes(mask)
		if 3 == len(mask.shape):
			# Single slice filling
			for sliceid in range(mask.shape[0]):
				mask[sliceid, :, :] = ndi.binary_fill_holes(mask[sliceid, :, :])

		# Estimate background 
		kwargs['dna_bg'] = imt.estimate_background(dna_ch, mask, seg_type)
		kwargs['sig_bg'] = imt.estimate_background(sig_ch, mask, seg_type)
		log += self.printout('Estimating background:', 2)
		log += self.printout('DNA channel: ' + str(kwargs['dna_bg']), 3)
		log += self.printout('Signal channel: ' + str(kwargs['sig_bg']), 3)

		# Filter objects XY size
		log += self.printout('Filtering objects XY size...', 2)
		sinter = stt.r_to_size(kwargs['radius_interval'], seg_type)
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

		# Filter objects Z size
		doFilterZsize = 0 != int(math.ceil(kwargs['min_z_size']))
		doFilterZsize = doFilterZsize and kwargs['an_type'] == const.AN_3D
		if doFilterZsize:
			log += self.printout('Filtering objects XY size...', 2)
			if kwargs['min_z_size'] > 1:
				kwargs['min_z_size'] = int(math.ceil(kwargs['min_z_size']))
			else:
				kwargs['min_z_size'] = kwargs['min_z_size'] * mask.shape[0]
				kwargs['min_z_size'] = int(math.ceil(kwargs['min_z_size']))
			msg = 'Minimum ' + str(kwargs['min_z_size']) + ' slices.'
			log += self.printout(msg, 3)

			# Identify objects Z size
			log += self.printout('Retrieving objects Z size...', 3)
			L = label(mask)
			zsizes = imt.get_objects_zsize(L)
			log += self.printout('Found ' + str(L.max()) + ' objects.', 3)
			
			# Select objects to be discarded
			torm = np.array(zsizes) < kwargs['min_z_size']
			torm = [ii for ii, x in enumerate(torm) if x]
			log += self.printout('Discarding ' + str(len(torm)) + ' objects.', 3)

			# Remove objects lower than minimum size
			mask = vt.rm_from_mask(L, torm)
			L = label(mask)

		# Save mask
		log += self.printout('Saving series object mask...', 2)

		# Plot
		fig = plt.figure()
		if 3 == len(mask.shape):
			plt.imshow(L.max(0).astype('u4'))
		else:
			plt.imshow(L.astype('u4'))
		plt.gca().get_xaxis().set_visible(False)
		plt.gca().get_yaxis().set_visible(False)
		plot.set_font_size(kwargs['font_size'])

		title = 'Nuclei in "' + kwargs['cond_name'] + '", ' + str(self.name)
		title += ' [' + str(L.max()) + ' objects]'
		plt.title(title)

		# Export as png
		fname = kwargs['out_dir'] + const.OUTDIR_MASK + kwargs['cond_name']
		fname += '.' + self.name + '.mask' + suffix + '.png'
		if kwargs['plotting']: plot.export(fname, 'png')

		# Close plot figure
		plt.close(fig)
		
		# Initialize nuclei
		log += self.printout('Bounding ' + str(L.max()) + ' nuclei...', 2)
		kwargs['logpath'] = self.logpath
		kwargs['i'] = i
		kwargs['thr'] = thr
		kwargs['series_id'] = self.n
		seq = range(1, L.max() + 1)
		self.nuclei = [Nucleus(n = n, mask = L == n, **kwargs) for n in seq]

		return((self, log))

	def adjust_options(self, read_only_dna = None, **kwargs):
		"""
		Adjusts options to be passed to the Nucleus class.
		Adds the following kwargs:
		 - series_name <string> series wrap name
		 - basedir <string> series wrap base directory
		 - dna_ch <numpy.array> image (dimensionality based on an_type)
		 - sig_ch <numpy.array> image (dimensionality based on an_type)

		@param:
		 - dna_names <tuple[string]> dna channel names
		 - sig_names <tuple[string]> signal channel names
		 - an_type <pyGPSeq.const> analysis type
		"""

		# Only work on dna channel
		if None == read_only_dna:
			read_only_dna = False

		# Add necessary options
		kwargs['series_name'] = self.name
		kwargs['basedir'] = self.basedir

		# Start log (used when verbosity is off)
		log = ""
		log += self.printout('Reading channels...', 2)

		# Read DNA channel
		dna_file = self.find_channel(kwargs['dna_names'])
		dna_ch = imread(os.path.join(self.basedir, dna_file[0]))
		dna_ch = imt.autoselect_time_frame(dna_ch)
		dna_ch = imt.slice_k_d_img(dna_ch, 3)

		if not read_only_dna:
			# Read Signal channel
			sig_file = self.find_channel(kwargs['sig_names'])
			sig_ch = imread(os.path.join(self.basedir, sig_file[0]))
			sig_ch = imt.autoselect_time_frame(sig_ch)
			sig_ch = imt.slice_k_d_img(sig_ch, 3)

		# Deconvolved images correction
		if 'rescale_deconvolved' in kwargs.keys():
			if kwargs['rescale_deconvolved']:
				# Get DNA scaling factor and rescale
				sf = imt.get_rescaling_factor(dna_file, **kwargs)
				dna_ch = (dna_ch / sf).astype('float')
				msg = 'Rescaled "' + dna_file[0] + '" [' + str(sf) + ']...'
				log += self.printout(msg, 3)

				if not read_only_dna:
					# Get Signal scaling factor and rescale
					sf = imt.get_rescaling_factor(sig_file, **kwargs)
					sig_ch = (sig_ch / sf).astype('float')
					msg = 'Rescaled "' + sig_file[0] + '" [' + str(sf) + ']...'
					log += self.printout(msg, 3)

		# Make Z-projection
		if kwargs['an_type'] in [const.AN_SUM_PROJ, const.AN_MAX_PROJ]:
			msg = 'Generating Z-projection [' + str(kwargs['an_type']) + ']...'
			log += self.printout(msg, 2)
			if 2 != len(dna_ch.shape):
				dna_ch = imt.mk_z_projection(dna_ch, kwargs['an_type'])
			if 2 != len(sig_ch.shape):
				sig_ch = imt.mk_z_projection(sig_ch, kwargs['an_type'])

		# Prepare output
		kwargs['dna_ch'] = dna_ch
		kwargs['sig_ch'] = sig_ch

		# Output
		return((kwargs, log))

	def get_nuclei_data(self, nuclei_ids, **kwargs):
		""" Retrieve a single nucleus from the current series. """

		# Read channel images
		kwargs, log = self.adjust_options(**kwargs)

		# Re-build mask with adaptive threshold
		if 1 < kwargs['adp_thr']:
			log += self.printout('Re-building adaptive threshold mask...', 2)

			# Make mask
			i = kwargs['dna_ch'].copy()

			# Make Z-projection
			if const.SEG_3D != kwargs['seg_type'] and 2 != len(i.shape):
				msg = 'Generating Z-projection [' + kwargs['seg_type'] + ']...'
				log += self.printout(msg, 3)
				i = imt.mk_z_projection(i, kwargs['seg_type'])

			# Global mask
			thr = threshold_otsu(i)
			msg = 'Thresholding image, global thr: ' + str(thr)
			log += self.printout(msg, 3)
			mask = imt.binarize(i, thr)

			# Local mask
			msg = 'Applying adaptive threshold to neighbourhood: '
			msg += str(kwargs['adp_thr'])
			log += self.printout(msg, 3)
			lmask = imt.threshold_adaptive(i, kwargs['adp_thr'])

			# Combine global and local masks
			log += self.printout('Combining global and local thresholds...', 3)
			mask = np.logical_and(mask, lmask)

		# Empty nuclear data array
		data = []
		for nucleus_id in nuclei_ids:
			# Select nucleus
			n = self.nuclei[nucleus_id -1]

			# Setup nucleus instance verbosity
			if not self.verbose:
				n.verbose = False

			# Retrieve nuclear data
			ndata, nlog = n.get_data(mask = mask, **kwargs)

			# Update log and save nuclear data
			log += nlog
			data.append(ndata)

		return((data, log))

	def export_nuclei(self, **kwargs):
		""" Export current series nuclei. """

		# Set output suffix
		if not 'suffix' in kwargs.keys():
			suffix = ''
		else:
			suffix = st.add_leading_dot(kwargs['suffix'])

		# Add necessary options
		self.printout('Current series: "' + self.name + '"...', 1)
		kwargs, log = self.adjust_options(**kwargs)

		# Export nuclei
		[n.export(**kwargs) for n in self.nuclei]
		
		# Produce log
		log = np.zeros(len(self.nuclei), dtype = const.DTYPE_NUCLEAR_SUMMARY)
		for l in [n.get_summary(**kwargs) for n in self.nuclei]:
			# Append nuclear data to the series log
			summary = [self.n]
			summary.extend(l)
			log[i, :] = summary

		# Export series log
		np.savetxt(kwargs['out_dir'] + self.name + '.summary' + suffix + '.csv',
			log, delimiter = ',', comments = '',
			header = ",".join([h for h in log.dtype.names]))

		return(log)

	def propagate_attr(self, key):
		""" Propagate attribute current value to every nucleus. """
		for i in range(len(self.nuclei)):
			self.nuclei[i][key] = self[key]

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

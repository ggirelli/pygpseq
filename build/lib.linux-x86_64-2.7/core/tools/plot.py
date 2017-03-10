# -*- coding: utf-8 -*-

""" Functions for the management of plots """

from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import rc
from matplotlib import ticker
from matplotlib.ticker import ScalarFormatter
import matplotlib.pyplot as plt
import numpy as np
import scipy

from .. import const

from . import path as pt, stat as stt

"""

 NOTES:
 
 - If the font cache is re-built every time that matplotlib is imported,
   try removing the ~/.cache/matplotlib/ folder.

"""

def set_font_size(size):
	rc('font', **{'size'   : size})

def density_with_range(density, fwhm_range,
	new_figure = None, show = None, close = None, **kwargs):
	"""
	Plot density curve and FWHM range.
	Density should be the output of .stat.calc_density
	"""

	# Default values
	if None == new_figure:
		new_figure = True
	if None == show:
		show = True
	if None == close:
		clos = False

	# New figure
	if new_figure:
		fig = plt.figure()

	# Plot
	plt.plot(density['x'], density['y'], 'k')
	plt.hold(True)
	plt.plot(fwhm_range, density['f'](fwhm_range), 'b.')
	plt.axvline(x = fwhm_range[0], color = 'red', linestyle = '-.')
	plt.axvline(x = fwhm_range[1], color = 'red', linestyle = '-.')
	plt.ticklabel_format(style = 'sci', axis = 'x', scilimits = (0, 0))
	plt.ticklabel_format(style = 'sci', axis = 'y', scilimits = (0, 0))

	# Adjust labels
	if 'xlab' in kwargs.keys():
		plt.xlabel(kwargs['xlab'])
	if 'ylab' in kwargs.keys():
		plt.ylabel(kwargs['ylab'])
	if 'title' in kwargs.keys():
		plt.title(kwargs['title'])

	# Show figure
	if new_figure and show:
		plt.show()

	# Close
	if new_figure and close:
		plt.close(fig)

def single_pixel_study(xs, ys, field_label, profile, nbins = None, **kwargs):
	"""
	Plots a study of single-pixel and general behaviour
	of the provided profile data.

	@param:
	 - xs <np.array> pixel (relative) distance from nuclear lamina
	 - ys <np.array> pixel channel intensity
	 - field_label <string> y-axis label
	 - nbins <int> study precision (opt, def 200)
	"""

	# CHECK PARAMS =============================================================

	# Default values
	if None == nbins:
		nbins = 200

	# PREPARE DATA =============================================================

	# XY range
	xyrange = [[xs.min(), xs.max()], [ys.min(), ys.max()]]

	# Bin data
	assigned_bins = np.digitize(xs,
		np.linspace(xs.min(), xs.max(), nbins))
	binned_ydat = [[] for i in range(nbins)]
	for bin_id in range(nbins):
		binned_ydat[bin_id] = ys[np.where(assigned_bins == bin_id)]

	# Calculate local density
	hh, locx, locy = scipy.histogram2d(xs, ys,
		range = xyrange, bins = [nbins, nbins])
	posx = np.digitize(xs, locx)
	posy = np.digitize(ys, locy)

	# Hide empty bins
	pp = (hh.T / hh.T.max(0))
	pp[:, hh.T.sum(0) == 0] = 0

	# MULTIPLOT ================================================================

	# Init figure --------------------------------------------------------------
	fig = plt.figure(figsize = [40, 20])
	rc('font', **{'size' : 15})

	# Boxplot ------------------------------------------------------------------
	ax = plt.subplot2grid((4, 2), (0, 0))
	plt.boxplot(binned_ydat, sym = 'k.')
	plt.ylabel(field_label)
	ax.get_xaxis().set_visible(False)
	plt.ticklabel_format(style = 'sci', axis = 'y', scilimits = (0, 0))

	# Mean plot ----------------------------------------------------------------
	plt.subplot2grid((4, 2), (1, 0))
	plt.plot(profile['mean_raw'], 'b')
	plt.plot(profile['mean_raw'], 'b.')
	plt.ylabel('mean(' + field_label + ')')
	plt.hold(True)
	plt.plot(profile['mean'], 'y')

	# Adjust ticks
	ax = plt.gca()
	ticks = ax.get_xticks() / float(nbins) * xyrange[0][1]
	ticks = [const.SCI_FORMAT % n for n in ticks]
	ax.set_xticklabels(ticks)
	plt.ticklabel_format(style = 'sci', axis = 'y', scilimits = (0, 0))

	# Median plot --------------------------------------------------------------
	plt.subplot2grid((4, 2), (2, 0))
	plt.plot(profile['median_raw'], 'r')
	plt.plot(profile['median_raw'], 'r.')
	plt.ylabel('median(' + field_label + ')')
	plt.hold(True)
	plt.plot(profile['median'], 'g')

	# Adjust ticks
	ax = plt.gca()
	ax.set_xticklabels(ticks)
	plt.ticklabel_format(style = 'sci', axis = 'y', scilimits = (0, 0))

	# Mode plot ----------------------------------------------------------------
	plt.subplot2grid((4, 2), (3, 0))
	plt.plot(profile['mode_raw'], 'c')
	plt.plot(profile['mode_raw'], 'c.')
	plt.ylabel('mode(' + field_label + ')')
	plt.hold(True)
	plt.plot(profile['mode'], 'm')
	plt.xlabel(kwargs['dlabel'])

	# Adjust ticks
	ax = plt.gca()
	ax.set_xticklabels(ticks)
	plt.ticklabel_format(style = 'sci', axis = 'y', scilimits = (0, 0))

	# Number of pixels plot ----------------------------------------------------
	plt.subplot2grid((4, 2), (0, 1))
	plt.plot(range(nbins), hh.sum(1))
	plt.ylabel('Number of pixels')
	plt.xlabel(kwargs['dlabel'])

	# Adjust ticks
	ax = plt.gca()
	ax.set_xticklabels(ticks)
	plt.ticklabel_format(style = 'sci', axis = 'y', scilimits = (0, 0))

	# Local density scatterplot ------------------------------------------------
	plt.subplot2grid((4, 2), (1, 1), rowspan = 2)
	profile_density_scatterplot(pp, field_label, nbins, xyrange, profile,
		dlabel = kwargs['dlabel'], new_figure = False)

	# STD plot -----------------------------------------------------------------
	plt.subplot2grid((4, 2), (3, 1))
	plt.plot(profile['std_raw'], 'r')
	plt.plot(profile['std_raw'], 'r.')
	plt.ylabel('std(' + field_label + ')')
	plt.hold(True)
	plt.plot(profile['std'], 'g')
	plt.xlabel(kwargs['dlabel'])

	# Adjust ticks
	ax = plt.gca()
	ax.set_xticklabels(ticks)
	plt.ticklabel_format(style = 'sci', axis = 'y', scilimits = (0, 0))

	# SINGLE PLOTS =============================================================

	if kwargs['plotting']:

		# Density scatterplot --------------------------------------------------
		
		# Retrieve tmp figure
		tmp_fig = profile_density_scatterplot(pp, field_label, nbins, xyrange,
			profile, dlabel = kwargs['dlabel'])

		# Get filename suffix
		suffix = field_label.lower().split(' ')[0].replace('/', '_')

		# Get condition description
		if kwargs['cond_name'] in kwargs['cdescr'].keys():
			suptitle = kwargs['cdescr'][kwargs['cond_name']]
		else:
			suptitle = kwargs['cond_name']
		plt.suptitle(suptitle)

		# Set output name
		fname = kwargs['outdir'] + const.OUTDIR_PNG_REPORT + kwargs['cond_name']
		fname += '.' + suffix + '.density_profile'

		# Add partial suffix
		if 'partial' in kwargs.keys():
			if kwargs['partial']:
				fname += '.part'

		# Add general suffix
		fname += kwargs['suffix'] + '.png'

		# Export
		export(fname, 'png')

		# Close tmp figure
		plt.close(tmp_fig)

		# Number of pixels plot  -----------------------------------------------
		
		# Init tmp figure
		tmp_fig = plt.figure()

		# Plot
		plt.plot(range(nbins), hh.sum(1))
		plt.ylabel('Number of pixels')
		plt.xlabel(kwargs['dlabel'])
		ticks = ax.get_xticks() / float(nbins) * xyrange[0][1]
		ticks = [const.SCI_FORMAT % n for n in ticks]
		plt.gca().set_xticklabels(ticks)
		plt.ticklabel_format(style = 'sci', axis = 'y', scilimits = (0, 0))

		# Set suptitle
		if 'partial' in kwargs.keys():
			if kwargs['partial']:
				suptitle += ' [partial volume '
				suptitle += str(kwargs['part_n_erosion']) + ']'
		plt.suptitle(suptitle)

		# Set output name
		fname = kwargs['outdir'] + const.OUTDIR_PNG_REPORT + kwargs['cond_name']
		fname += '.' + suffix + '.npx'

		# Add partial suffix
		if 'partial' in kwargs.keys():
			if kwargs['partial']:
				fname += '.part'

		# Add general suffix
		fname += kwargs['suffix'] + '.png'

		# Export
		export(fname, 'png')

		# Close tmp figure
		plt.close(tmp_fig)

		# # CV plot
		# tmp_fig = plt.figure()
		# plt.plot(np.array(std_dat) / np.array(mean_dat), 'r')
		# plt.ylabel('CV(' + field_label + ')')
		# plt.hold(True)
		# plt.plot(np.array(smooth_std) / np.array(smooth_mean), 'g')
		# fname = kwargs['outdir'] + const.OUTDIR_PNG_REPORT
		# fname += kwargs['cond_name'] + '.' + suffix + '.cv'
		# if 'partial' in kwargs.keys():
		# 	if kwargs['partial']:
		# 		fname += '.part'
		# fname += kwargs['suffix'] + '.png'
		# export(fname, 'png')
		# plt.close(tmp_fig)

	# Return MULTIPLOT canvas
	return(fig)

def profile_density_scatterplot(pp, field_label, nbins, xyrange,
	profile, new_figure = None, dlabel = None):
	"""
	Plot profile pixel density scatterplot.

	@param:
	 - pp <np.array> bin-normalized densi
	 - field_label <string> y-axis label
	 - nbins <int> number of bins
	 - xyrange <list> x/y-axis range [min, max]
	 - profile <list(np.array)> mean, median, mode, std (smoothed and raw)
	 - new_figure <boolean> whether to plot in a new figure (and return it)
	"""

	# CHECK PARAMS =============================================================

	if None == new_figure:
		new_figure = True

	if new_figure:
		fig = plt.figure()

	# PLOT =====================================================================
	
	plt.imshow(pp, origin = 'lower', interpolation = 'none',
		cmap = 'gray')

	cbar = plt.colorbar()
	cbar.set_label('Bin-relative local density', rotation = 270, labelpad = 20)

	plt.xlabel(dlabel)
	plt.ylabel(field_label)
	
	plt.hold(True)

	plt.plot(range(nbins), nbins *
		(profile['median'] - xyrange[1][0]) / (xyrange[1][1] - xyrange[1][0]),
		'g', linewidth = 2)
	plt.plot(range(nbins), nbins *
		(profile['mean'] - xyrange[1][0]) / (xyrange[1][1] - xyrange[1][0]),
		'y', linewidth = 2)
	plt.plot(range(nbins), nbins *
		(profile['mode'] - xyrange[1][0]) / (xyrange[1][1] - xyrange[1][0]),
		'm', linewidth = 2)

	plt.xlim([0, nbins - 1])
	plt.ylim([0, nbins - 1])

	# Adjust ticks
	ax = plt.gca()

	ticks = ax.get_xticks() / float(nbins) * xyrange[0][1]
	ticks = [const.SCI_FORMAT % n for n in ticks]
	ax.set_xticklabels(ticks)

	ticks = ax.get_yticks() / float(nbins) * xyrange[1][1]
	ticks = [const.SCI_FORMAT % n for n in ticks]
	ax.set_yticklabels(ticks)

	if new_figure:
		return(fig)

def single_condition_profiles(profiles, font_size = None, color = None,
	hspace = None, wspace = None, n_nuclei = None, new_figure = None,
	cdescr = None, yfield = None, **kwargs):
	"""
	Saves three single condition profiles.

	@param:
	 - profiles <{dtype:{x:float, y:float}, n: int, condition:string}>
	 - font_size <int>
	 - hspace <float> horizontal intra-plot spacing
	 - wspace <float> vertical intra-plot spacing
	 - n_nuclei <int> number of nuclei
	"""
	
	# SET PARAMS ===============================================================

	if None == font_size:
		font_size = 12

	if None == n_nuclei:
		n_nuclei = '~unknown~'

	if None == new_figure:
		new_figure = True

	if None != cdescr:
		if profiles['condition'] in cdescr.keys():
			profiles['condition'] = cdescr[profiles['condition']]

	if None == yfield:
		yfield = 'mean'

	if None == wspace:
		wspace = .4

	if None == hspace:
		hspace = .4

	# Plot options
	opt = {}
	if None != color:
		opt['c'] = color

	# PLOT =====================================================================

	# Init canvas --------------------------------------------------------------
	if new_figure:
		fig = plt.figure(figsize = (12, 8))

		# Main title
		title = 'GPSeq profiles for condition "'
		title += profiles['condition'] + '" (n.nuclei = ' + str(n_nuclei) + ')'

		# Add user-defined comment
		if 'title_comment' in kwargs.keys():
			title += ' [' + kwargs['title_comment'] + ']'

		# Additional information
		title += ' [sigma = ' + str(kwargs['sigma']) + ']'
		title += ' [nbins = ' + str(kwargs['nbins']) + ']'

		# Add suptitle
		plt.suptitle(title)

	# Setup subplots spacing
	plt.subplots_adjust(wspace = wspace, hspace = hspace)
	
	# Number of steps (precision)
	nsteps = len(profiles['ratio']['x'])

	# DNA profile --------------------------------------------------------------
	plt.subplot(3, 2, 1)
	plt.hold(True)
	plt.plot(profiles['dna']['x'], profiles['dna'][yfield],
		label = profiles['condition'], **opt)
	plt.xlabel(kwargs['dlabel'])
	plt.ylabel('DNA [a.u.]')
	plt.ticklabel_format(style = 'sci', axis = 'y', scilimits = (0, 0))
	set_font_size(font_size)

	# DNA first derivative -----------------------------------------------------
	plt.subplot(3, 2, 2)
	plt.hold(True)
	plt.plot(profiles['dna']['x'][1:nsteps], np.diff(profiles['dna'][yfield]),
		label = profiles['condition'], **opt)
	plt.axhline(0, color = 'k', linestyle = '-.')
	plt.xlabel(kwargs['dlabel'])
	plt.ylabel('diff(DNA [a.u.])')
	plt.ticklabel_format(style = 'sci', axis = 'y', scilimits = (0, 0))
	set_font_size(font_size)

	# Signal profile -----------------------------------------------------------
	plt.subplot(3, 2, 3)
	plt.hold(True)
	plt.plot(profiles['sig']['x'], profiles['sig'][yfield],
		label = profiles['condition'], **opt)
	plt.xlabel(kwargs['dlabel'])
	plt.ylabel('Signal [a.u.]')
	plt.ticklabel_format(style = 'sci', axis = 'y', scilimits = (0, 0))
	set_font_size(font_size)

	# Signal first derivative --------------------------------------------------
	plt.subplot(3, 2, 4)
	plt.hold(True)
	plt.plot(profiles['sig']['x'][1:nsteps], np.diff(profiles['sig'][yfield]),
		label = profiles['condition'], **opt)
	plt.axhline(0, color = 'k', linestyle = '-.')
	plt.xlabel(kwargs['dlabel'])
	plt.ylabel('diff(Signal [a.u.])')
	plt.ticklabel_format(style = 'sci', axis = 'y', scilimits = (0, 0))
	set_font_size(font_size)

	# Ratio profile ------------------------------------------------------------
	plt.subplot(3, 2, 5)
	plt.hold(True)
	plt.plot(profiles['ratio']['x'], profiles['ratio'][yfield],
		label = profiles['condition'], **opt)
	plt.xlabel(kwargs['dlabel'])
	plt.ylabel('Signal/DNA')
	plt.ticklabel_format(style = 'sci', axis = 'y', scilimits = (0, 0))
	set_font_size(font_size)

	# Ratio first derivative ---------------------------------------------------
	plt.subplot(3, 2, 6)
	plt.hold(True)
	plt.plot(profiles['ratio']['x'][1:nsteps],
		np.diff(profiles['ratio'][yfield]),
		label = profiles['condition'], **opt)
	plt.axhline(0, color = 'k', linestyle = '-.')
	plt.xlabel(kwargs['dlabel'])
	plt.ylabel('dif(Signal/DNA)')
	plt.ticklabel_format(style = 'sci', axis = 'y', scilimits = (0, 0))
	set_font_size(font_size)

	# Close figure -------------------------------------------------------------
	if new_figure:
		return(fig)

def multi_condition_profiles(profiles, font_size = None,
	hspace = None, wspace = None, yfield = None, **kwargs):
	"""
	Saves three single condition profiles.

	@param:
	 - profiles <[{dtype:{x:float, y:float}, n:int, condition:string}]>
	 - font_size <int>
	 - hspace <float> horizontal intra-plot spacing
	 - wspace <float> vertical intra-plot spacing
	"""

	if None == yfield:
		yfield = 'mean'

	fig = plt.figure(figsize = (12, 10))
	
	# Setup subplots spacing
	for i in range(len(profiles)):
		profile = profiles[i]
		n_nuclei = profile['n']

		single_condition_profiles(profile, color = const.PALETTE[i],
			n_nuclei = n_nuclei, new_figure = False, yfield = yfield, **kwargs)

	# Add legend
	plt.subplot(3, 2, 1)
	plt.legend(bbox_to_anchor = (0., 1.12, 1., .102), loc = 3,
		ncol = 2, mode = "expand", borderaxespad = 0.,
		fontsize = 'xx-small')

	# Add title
	title = 'GPSeq profiles [' + yfield + ']'
	if 'title_comment' in kwargs.keys():
		title += ' [' + kwargs['title_comment'] + ']'
	title += ' [sigma = ' + str(kwargs['sigma']) + ']'
	title += ' [nbins = ' + str(kwargs['nbins']) + ']'
	plt.suptitle(title)

	return(fig)

def boxplot(fig, subplot, data, labels, widths, ylabel, nnuclei,
	xaxis1 = None, xaxis2 = None, ylab_pos = None, forceSciNotation = None):
	""" Generate single condition boxplot """

	# Set default xaxis visibility
	if None == xaxis1:
		xaxis1 = True
	if None == xaxis2:
		xaxis2 = True
	if None == ylab_pos:
		ylab_pos = 'left'
	if None == forceSciNotation:
		forceSciNotation = False

	# Create double xaxis
	ax1 = fig.add_subplot(subplot)
	ax2 = ax1.twiny()

	# Produce boxplot
	bp = ax1.boxplot(data, labels = labels,
		sym = '.', widths = widths, patch_artist = True)

	# Change boxplot appearance
	[b.set(color = '#000000') for b in bp['boxes']]
	[b.set(facecolor = '#ffffff') for b in bp['boxes']]
	[w.set(color = '#000000') for w in bp['whiskers']]
	[c.set(color = '#000000') for c in bp['caps']]
	[m.set(color = 'red') for m in bp['medians']]
	[m.set(linewidth = 2) for m in bp['medians']]
	[f.set(color = '#000000') for f in bp['fliers']]

	# Change axes labels
	ax1.set_xlabel('Condition')
	ax1.set_ylabel(ylabel)
	ax2.set_xlim(ax1.get_xlim())
	ax2.set_xticks(np.asarray(range(len(widths))) + 1)
	ax2.set_xticklabels(nnuclei)
	ax2.set_xlabel("Number of nuclei")

	# Force scientific notation
	if forceSciNotation:
		ax1.ticklabel_format(style = 'sci', axis = 'y', scilimits = (0, 0))

	# Move y axis
	ax1.yaxis.set_label_position(ylab_pos)

	# Hide unnecessary axes
	if not xaxis1:
		ax1.get_xaxis().set_visible(False)
	if not xaxis2:
		ax2.get_xaxis().set_visible(False)

def multi_condition_boxplot(profiles, summaries, merged, **kwargs):
	""" Generate multi condition boxplots """

	# Create new plot figure
	fig = plt.figure(figsize = (5 + 3 * len(profiles), 12))
	
	# Calculate boxplot relative widths
	bp_widths = np.asarray([p['n'] for p in profiles], dtype='float')
	bp_widths = bp_widths / bp_widths.max() * 0.9

	if const.AN_3D == kwargs['an_type']:
		size_ylab = 'Volume [vx]'
		shape_ylab = 'Sphericity'
	else:
		size_ylab = 'Area [px]'
		shape_ylab = 'Circularity'

	# Plot size
	boxplot(fig, 321,
		[s['size'] for s in summaries],
		[p['condition'] for p in profiles], 
		bp_widths.tolist(), size_ylab + '  [per nucleus]',
		[p['n'] for p in profiles],
		xaxis1 = False,
		forceSciNotation = True
	)

	# Plot shape
	boxplot(fig, 322,
		[s['shape'] for s in summaries],
		[p['condition'] for p in profiles], 
		bp_widths.tolist(), shape_ylab + ' [per nucleus]',
		[p['n'] for p in profiles],
		xaxis1 = False, ylab_pos = 'right',
		forceSciNotation = True
	)

	# Plot intensity average
	boxplot(fig, 323,
		[s['meanI'] for s in summaries],
		[p['condition'] for p in profiles], 
		bp_widths.tolist(), 'mean(DNA [a.u.]) [per nucleus]',
		[p['n'] for p in profiles],
		xaxis1 = False, xaxis2 = False,
		forceSciNotation = True
	)

	# Plot intensity sum
	boxplot(fig, 324,
		[s['sumI'] for s in summaries],
		[p['condition'] for p in profiles], 
		bp_widths.tolist(), 'sum(DNA [a.u.]) [per nucleus]',
		[p['n'] for p in profiles],
		xaxis1 = False, xaxis2 = False,
		ylab_pos = 'right',
		forceSciNotation = True
	)

	# Calculate boxplot relative widths
	bp_widths = np.asarray([m.shape[0] for m in merged], dtype='float')
	bp_widths = bp_widths / bp_widths.max() * 0.9

	# Plot single-pixel DNA intensity
	boxplot(fig, 325,
		[m['dna'] for m in merged],
		[p['condition'] for p in profiles], 
		bp_widths.tolist(), 'DNA [a.u.] [per px]',
		[m.shape[0] for m in merged],
		xaxis2 = False,
		forceSciNotation = True
	)

	# Plot single-pixel Signal intensity
	boxplot(fig, 326,
		[m['sig'] for m in merged],
		[p['condition'] for p in profiles], 
		bp_widths.tolist(), 'Signal [a.u.] [per px]',
		[m.shape[0] for m in merged],
		xaxis2 = False, ylab_pos = 'right',
		forceSciNotation = True
	)

	return(fig)

def multi_condition_single_boxplot(profiles, data, data_field,
	bp_widths, ylab, out_png, **kwargs):
	fig = plt.figure(figsize = (5 + 3 * len(profiles), 4))
	boxplot(fig, 111, [d[data_field] for d in data],
		[p['condition'] for p in profiles], bp_widths.tolist(), ylab,
		[p['n'] for p in profiles], forceSciNotation = True
	)
	export(out_png + 'boxplot.' + data_field + kwargs['suffix'] + '.png', 'png')
	plt.close(fig)

def export(path, exp_format = None, plt_obj = None):
	"""
	Exports the current plot to the provided path.
	The formats should be specified.
	Also, a plot object can be exported instead of the current one.
	"""

	if None == exp_format:
		exp_format = 'pdf'

	if exp_format in ['pdf', 'png', 'jpg']:

		# Check presence of the extension
		path = pt.add_extension(path, '.' + exp_format)

		# Save plot
		if exp_format == 'pdf':
			pp = PdfPages(path)
			plt.savefig(pp, format = exp_format)
			pp.close()
		else:
			plt.savefig(path, format = exp_format)

def get_nsf_label(nsfi, seg_type):
	"""
	Get the proper plot label for the selected nuclear selection feature,
	based on the current segmentation type
	"""

	# Check that the provided nuclear selection feature index exists
	if not nsfi in range(len(const.NSEL_FIELDS)):
		return(None)

	# Retrieve label
	label = const.NSEL_LABELS[nsfi]

	# Change label if required
	if 'auto' == label:
		if nsfi == const.NSEL_SIZE:
			if seg_type == const.SEG_3D:
				label = 'Volume [vx]'
			else:
				label = 'Area [px]'
		elif nsfi == const.NSEL_SHAPE:
			if seg_type == const.SEG_3D:
				label = 'Sphericity'
			else:
				label = 'Circularity'

	# Output
	return(label)

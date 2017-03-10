# -*- coding: utf-8 -*-

import pickle as cp
import datetime
import multiprocessing
import os
import pkg_resources
import time
import warnings

from jinja2 import Environment, FileSystemLoader
import logging
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from weasyprint import HTML, CSS

from . import const

from .tools import path as pt, io as iot, plot, string as st

from .wraps.analysis import Analyzer
from .wraps.condition import Condition

class Main(Analyzer):
	""" GPSeq data manager (extends Analyzer) """

	# Class-bound attributes
	__version__ = const.VERSION
	
	# Base directory path with data
	basedir = './'

	# Output directory
	outdir = './output/'

	# Steps to be skipped
	# 1	  : skip instantiation (need gpi.inst.cpickle, otherwise unskips)
	# 2	  : skip segmentation (need gpi.seg.cpickle, otherwise unskips)
	# 3	  : skip analysis (need gpi.an.cpickle, otherwise unskips)
	# 3.5 : skip single-nuclei boxplot (end of step 3)
	# 4	  : final plots
	# 5	  : final report
	skip = []

	# Wrapped conditions
	conds = []

	# Series file attributes
	ext = '.tif'

	# Plot features
	font_size = 8

	# Parallelization
	ncores = 1

	# Series regexp
	reg = '^(?P<' + const.REG_CHANNEL_NAME + '>[^/]*)'
	reg += '\.(?P<' + const.REG_CHANNEL_ID + '>channel[0-9]+)'
	reg += '\.(?P<' + const.REG_SERIES_ID + '>series[0-9]+)'
	reg += '(?P<' + const.REG_EXT + '>\.tif)$'

	# Debug mode, effects:
	#  - Save intermediate nuclear tif images of DNA, signal, distance and mask.
	debugging = False

	# Whether to save plots or not
	plotting = True

	# Output suffix (a leading dot is automatically added)
	suffix = ''

	# User notes
	notes = '...'

	def __init__(self, **kwargs):
		"""
		@param:
		 - ncores <int> number of cores (opt, def: 1)
		 - font_size <int>
		"""
		super(Main, self).__init__()

		# Read additional parameters
		if 'ncores' in kwargs.keys():
			self.ncores = kwargs['ncores']
		if 'font_size' in kwargs.keys():
			self.font_size = kwargs['font_size']

		self.printout('', 0)

	def run(self, **kwargs):
		""" Runs the GPSeq manager """

		# INIT =================================================================

		# Suppress RuntimeWarning(s)
		warnings.simplefilter("ignore", category = FutureWarning)
		warnings.simplefilter("ignore", category = RuntimeWarning)
		warnings.simplefilter("ignore", category = UserWarning)

		# Check number of cores
		if self.ncores > multiprocessing.cpu_count():
			self.ncores = multiprocessing.cpu_count()
			msg = 'Decreased core number to maximum allowed: %i' % self.ncores
			msg += '\nPlease, don\'t ask for the impossible... ಠ_ಠ'
			self.printout(msg, -1)
			self.printout('', 0)

		# Warn for log freezing if parallelization is on
		if self.ncores > 1:
			msg = 'pyGPSeq log might freeze due to parallelization.\n'
			msg += 'But do NOT dispair. Everything will work out, eventually...'
			msg += ' Just be patient.\n┬─┬﻿ ノ( ゜-゜ノ)'
			self.printout(msg, -1)
			self.printout('', 0)

		# Check parameters
		start_time = time.time()
		kwargs['suffix'] = st.add_leading_dot(self.suffix)
		self.basedir = pt.add_trailing_slash(self.basedir)
		self.outdir = pt.add_trailing_slash(self.outdir)

		# Create basedir and outdir if missing
		for d in [self.basedir, self.outdir]:
			if not os.path.isdir(d):
				os.makedirs(d)

		# Make profile output folder
		for d in [const.OUTDIR_PDF, const.OUTDIR_PNG,
			const.OUTDIR_PNG_REPORT, const.OUTDIR_MASK]:
			if not os.path.exists(self.outdir + d):
				os.makedirs(self.outdir + d)
		debug_dir = self.outdir + const.OUTDIR_DEBUG
		if self.debugging and not os.path.exists(debug_dir):
			os.makedirs(debug_dir)

		# Plot font-size
		if not 'font_size' in kwargs.keys():
			kwargs['font_size'] = self.font_size

		# Number of cores for parallelization
		if not 'ncores' in kwargs.keys():
			kwargs['ncores'] = self.ncores

		# Output directory
		kwargs['out_dir'] = self.outdir

		# Distance field for profiles
		if self.normalize_distance:
			kwargs['dfield'] = 'dnorm'
			kwargs['dlabel'] = 'Relative distance from nuclear lamina'
			kwargs['dlabel'] += ' [a.u.]'
		else:
			kwargs['dfield'] = 'd'
			kwargs['dlabel'] = 'Distance from nuclear lamina [a.u.]'

		# Update kwargs with self attributes
		kwargs.update([(n, getattr(self, n)) for n in dir(self)
			if type(getattr(self, n)) in const.KWARGS_TYPELIST and
			not n.startswith('__') and not n in const.KWARGS_AVOIDLIST])

		# SINGLE STEPS =========================================================

		# INSTANTIATION --------------------------------------------------------
		# Check whether to skip instantiation
		if self.is_skipped(1):
			fname = self.outdir + 'gpi.inst' + kwargs['suffix'] + '.cpickle'
			if os.path.exists(fname):
				self.printout('Skipping instantiation...', 0)
				self.printout('Loading dumped instance...\n', 0)

				# Load
				f = open(fname, 'rb')
				keys = list(const.PARAM_SEG)
				keys.extend(list(const.PARAM_AN))
				self = self.load(f, keys)
				f.close()

				# Dump
				f = open(fname, 'wb')
				cp.dump(self, f)
				f.close()
			else:
				self.printout('Unskipping instantiation...', 0)
				self.unskip(1)
		
		if not self.is_skipped(1):
			# Run instantiation if not skipped
			self.run_initialization(**kwargs)

			# Dump
			fname = self.outdir + 'gpi.inst' + kwargs['suffix'] + '.cpickle'
			f = open(fname, 'wb')
			cp.dump(self, f)
			f.close()

		# SEGMENTATION ---------------------------------------------------------
		# Check whether to skip segmentation
		if self.is_skipped(2):
			fname = self.outdir + 'gpi.seg' + kwargs['suffix'] + '.cpickle'
			if os.path.exists(fname):
				self.printout('Skipping segmentation...', 0)
				self.printout('Loading dumped instance...\n', 0)

				# Load
				f = open(fname, 'rb')
				self = self.load(f, const.PARAM_AN)
				f.close()

				# Dump
				f = open(fname, 'wb')
				cp.dump(self, f)
				f.close()
			else:
				self.printout('Unskipping segmentation...', 0)
				self.unskip(2)
		
		if not self.is_skipped(2):
			# Run segmentation if not skipped
			self.run_segmentation(**kwargs)

			# Dump
			fname = self.outdir + 'gpi.seg'+ kwargs['suffix'] +'.cpickle'
			f = open(fname, 'wb')
			cp.dump(self, f)
			f.close()

		# ANALYSIS -------------------------------------------------------------
		# Check whether to skip analysis
		if self.is_skipped(3):
			fname = self.outdir + 'gpi.an' + kwargs['suffix'] + '.cpickle'
			if os.path.exists(fname):
				self.printout('Skipping analysis...', 0)
				self.printout('Loading dumped instance...\n', 0)

				# Load
				f = open(fname, 'rb')
				self, profiles, sumd = self.load(f)
				f.close()

				# Dump
				f = open(fname, 'wb')
				cp.dump((self, profiles, sumd), f)
				f.close()
			else:
				self.printout('Unskipping analysis...', 0)
				self.unskip(3)
		
		if not self.is_skipped(3):
			# Run analysis if not skipped
			profiles, sumd, md = self.run_analysis(**kwargs)

			# Dump
			fname = self.outdir + 'gpi.an' + kwargs['suffix'] + '.cpickle'
			f = open(fname, 'wb')
			cp.dump((self, profiles, sumd), f)
			f.close()

			# Generate general boxplots if not skipped
			if not self.is_skipped(3.5):
				self.mk_general_boxplots(profiles, sumd, md, **kwargs)

		# FINAL PLOTS ----------------------------------------------------------
		# Produce final plots if not skipped
		if not self.is_skipped(4):
			self.mk_general_plots(profiles, **kwargs)
		else:
			self.printout('Skipping final plots...', 0)

		# FINAL REPORT ---------------------------------------------------------
		# Execution end
		end_time = time.time()

		# Generate final report if not skipped
		if not self.is_skipped(5):
			self.printout('Generating final report...', 0)
			self.mk_report(start_time, end_time)
		else:
			self.printout('Skipping final report...', 0)

		# CONCLUSION ===========================================================

		# Final remarks
		self.printout('', 0)
		self.printout('~ took %s s ~' % (end_time - start_time), 0)
		self.printout('\n ~ FIN ~', 0)
		self.printout('\n└[∵┌] └[ ∵ ]┘ [┐∵]┘\n\n', 0)

		return(self)

	def run_initialization(self, **kwargs):
		""" Initialize run """

		self.printout('Starting GPSeq manager...', 0)

		self.printout('Looking into folder: "' + self.basedir + '"', 0)
		self.printout('', 0)

		# Instantiate OOP architecture
		self.printout('* Building OOP Architecture *', 0)

		# Select condition folders
		self.conds = pt.select_folders(self.basedir, self.ext)
		msg = 'Found ' + str(len(self.conds)) + ' condition folder(s)...'
		self.printout(msg, 0)
		self.printout('', 0)

		# Instantiate conditions
		self.conds = [Condition(c, main = self) for c in self.conds]
		self.printout('', 0)

	def run_segmentation(self, **kwargs):
		""" Run segmentation. """

		# Check analysis/segmentation types
		self.check_anseg_types()

		# Identify nuclei
		self.printout('* Looking for nuclei *', 0)
		self.printout('', 0)
		[c.find_nuclei(**kwargs) for c in self.conds]
		self.printout('', 0)

	def run_analysis(self, **kwargs):
		""" Run analysis. """

		# Retrieve and save nuclear data
		self.printout('* Retrieving nuclear data *', 0)
		self.printout('', 0)

		# [{dtype:{x:float, y:float}, n:int, condition:string}]
		data = [c.analyze_nuclei(**kwargs) for c in self.conds]
		profiles = [d[0] for d in data]
		sumd = [d[1] for d in data]
		md = [d[2] for d in data]

		return((profiles, sumd, md))

	def mk_general_boxplots(self, profiles, sumd, md, **kwargs):
		""" Generate general boxplots. """
		
		# Common destinations
		out_pdf = self.outdir + 'out_pdf/'
		out_png = self.outdir + 'out_png/'
		
		# Multi-condition single-nucleus boxplots
		self.printout('Preparing general single-nuclei boxplot...', 0)
		fig = plot.multi_condition_boxplot(profiles, sumd, md, **kwargs)

		# Export PDF
		fname = out_pdf + 'boxplots' + kwargs['suffix'] + '.pdf'
		if self.plotting: plot.export(fname, 'pdf')

		# Export PNG
		fname = out_png + 'boxplots' + kwargs['suffix'] + '.png'
		if self.plotting: plot.export(fname, 'png')

		# Close figure
		plt.close(fig)

		# Calculate boxplot relative widths
		bp_widths = [p['n'] for p in profiles]
		bp_widths = np.asarray(bp_widths, dtype = 'float')
		bp_widths = bp_widths / bp_widths.max() * 0.9

		# Per nucleus boxplots
		bpitems = [
			# Per nucleus size boxplots
			(plot.get_nsf_label(const.NSEL_SIZE, kwargs['seg_type'])
				+ ' [per nucleus]', 'size'),

			# Per nucleus shape boxplots
			(plot.get_nsf_label(const.NSEL_SHAPE, kwargs['seg_type'])
				+ ' [per nucleus]', 'shape'),

			# Per nucleus intensity average boxplots
			('mean(DNA [a.u.]) [per nucleus]', 'meanI'),

			# Per nucleus intensity sum boxplots
			('sum(DNA [a.u.]) [per nucleus]', 'sumI')
		]
		for (ylab, field) in bpitems:
			plot.multi_condition_single_boxplot(profiles, sumd, field,
				bp_widths, ylab, self.outdir + const.OUTDIR_PNG_REPORT,
				**kwargs)

		# Calculate boxplot relative widths
		bp_widths = np.asarray([m.shape[0] for m in md], dtype='float')
		bp_widths = bp_widths / bp_widths.max() * 0.9

		# Per pixel boxplots
		bpitems = [
			# Per pixel DNA intensity
			('DNA [a.u.] [per pixel]', 'dna'),

			# Per pixel Signal intensity
			('Signal [a.u.] [per pixel]', 'sig')
		]
		for (ylab, field) in bpitems:
			plot.multi_condition_single_boxplot(profiles, md, field,
				bp_widths, ylab, self.outdir + const.OUTDIR_PNG_REPORT,
				**kwargs)

	def mk_general_plots(self, profiles, **kwargs):
		""" Generate final plots """

		# Common destinations
		out_pdf = self.outdir + 'out_pdf/'
		out_png = self.outdir + 'out_png/'

		for yfield in ['mean', 'median', 'mode']:
			msg = 'Preparing multi-condition profiles plot [' + yfield + ']...'
			self.printout(msg, 0)

			# Plot
			fig = plot.multi_condition_profiles(profiles,
				yfield = yfield, **kwargs)

			# Export PDF
			common_name = 'profiles.' + yfield + kwargs['suffix']
			fname = out_pdf + common_name + '.pdf'
			if self.plotting: plot.export(fname, 'pdf')

			# Export PNG
			fname = out_png + common_name + '.png'
			if self.plotting: plot.export(fname, 'png')

			# Close figure
			plt.close(fig)

		# Export profiles to CSV
		fname = self.outdir + 'profiles' + kwargs['suffix'] + '.csv'
		if self.plotting: iot.export_profiles(profiles, fname)

		# Partial volume profile
		if const.AN_3D == self.an_type:
			self.printout('Selecting partial-volume pixels...', 0)
			part_profiles = [p['part'] for p in profiles]

			for yfield in ['mean', 'median', 'mode']:
				msg = 'Preparing partial volume multi-condition profiles plot ['
				msg += yfield + ']...'
				self.printout(msg, 0)

				# Plot
				fig = plot.multi_condition_profiles(part_profiles,
					yfield = yfield, title_comment = 'partial_volume '
					+ str(kwargs['part_n_erosion']), **kwargs)

				# Export PDF
				common_name = 'profiles.part.' + yfield + kwargs['suffix']
				fname = out_pdf + common_name + '.pdf'
				if self.plotting: plot.export(fname,'pdf')

				# Export PNG
				fname = out_png + common_name + '.png'
				if self.plotting: plot.export(fname,'png')

				# Close
				plt.close(fig)

			# Export profiles to CSV
			fname = self.outdir + 'profiles.part' + kwargs['suffix'] + '.csv'
			iot.export_profiles(part_profiles, fname)

	def mk_report(self, start_time, end_time, path = None):
		""" Produce PDF report. """

		# From time to timestamp
		start_time = datetime.datetime.fromtimestamp(start_time)
		start_time = start_time.strftime('%Y-%m-%d %H:%M:%S')
		end_time = datetime.datetime.fromtimestamp(end_time)
		end_time = end_time.strftime('%Y-%m-%d %H:%M:%S')

		# Load template
		gpdir = pkg_resources.resource_filename(const.PACK_NAME, '') + '/static/'
		env = Environment(loader = FileSystemLoader(gpdir))
		env.filters['type'] = type
		env.filters['has_key'] = lambda a, b: b in a.keys()
		template = env.get_template("report_template.html")

		# Prepare variables to fill the template
		tempv = {
			'starttime' : start_time,
			'endtime' : end_time,
			'basedir' : self.basedir,
			'outdir' : self.outdir,
			'logpath' : self.logpath,
			'reg' : self.reg,
			'ext' : self.ext,
			'skip' : [const.STEP_DESCR[i - 1] for i in self.skip
				if i in range(len(const.STEP_DESCR))],
			'ncores' : self.ncores,
			#'correctca' : self.correctCA,
			'verbose' : self.verbose,
			'debugging' : self.debugging,
			'suffix' : self.suffix,
			'dna_names' : self.dna_names,
			'sig_names' : self.sig_names,
			'plotting' : self.plotting,
			'fontsize' : self.font_size,
			'nbins' : self.nbins,
			'sigma' : self.sigma,
			'seg_type' : const.SEG_LABELS[self.seg_type],
			'adp_thr' : self.adp_thr,
			'radius_interval' : self.radius_interval,
			'min_z_size' : self.min_z_size,
			'offset' : self.offset,
			'rm_z_tips' : self.rm_z_tips,
			'an_type' : const.AN_LABELS[self.an_type],
			'aspect' : self.aspect,
			'nsf' : [const.NSEL_NAMES[i] for i in self.nsf
				if i in range(len(const.NSEL_NAMES))],
			'part_n_erosion' : self.part_n_erosion,
			'norm_d' : self.normalize_distance,
			'rescale_deconvolved' : self.rescale_deconvolved,
			'notes' : self.notes,
			'conds' : self.conds,
			'cnuclei' : [sum([len(s.nuclei) for s in c.series])
				for c in self.conds],
			'cdescr' : self.cdescr
		}

		# Escape characters
		for (k,v) in tempv.items():
			if type('') == type(v):
				tempv[k] = v.replace('<', '&lt;').replace('>', '&gt;')

		# Fill template
		html_out = template.render(tempv)

		# Hide CSS warnings
		logger = logging.getLogger('weasyprint')
		logger.handlers = []
		logger.addHandler(logging.FileHandler('/tmp/weasyprint.log'))

		# Output
		suffix = datetime.datetime.fromtimestamp(time.time())
		suffix = suffix.strftime('%Y-%m-%d %H:%M:%S')
		fname = self.outdir + 'report.' + suffix + '.pdf'
		HTML(string = html_out).write_pdf(fname)

		# f = open(self.outdir + 'report.' + suffix + '.htm', 'w')
		# f.write(html_out)
		# f.close()		

	def is_skipped(self, step):
		return(step in self.skip)

	def unskip(self, step):
		self.skip = [i for i in self.skip if step != i]

	def load(self, f, more_keys = None):
		"""
		Re-loads a dumped pyGPSeq instance, keeping some parameters
		from the current one.
		"""

		# Parameters to keep from current instance
		keys = list(const.PARAM_STATIC)
		if not None == more_keys:
			keys.extend(list(more_keys))
		vals = [self[k] for k in keys]

		# Load
		loaded = cp.load(f)

		# If the loaded cpickle bundle contained a Main instance
		if type(loaded) == type(self):
			# Set parameters back
			for i in range(len(keys)):
				loaded[keys[i]] = vals[i]

			# Propagate parameters
			for param in const.PARAM_PROPAGATE:
				loaded.propagate_attr(param)
		elif not 0 == len(loaded):
			# Identify the Main instance
			i = filter(lambda i: type(loaded[i]) == type(self),
				range(len(loaded)))[0]
			tmp = loaded[i]

			# Set parameters back
			for j in range(len(keys)):
				tmp[keys[j]] = vals[j]

			# Propagate parameters
			for param in const.PARAM_PROPAGATE:
				tmp.propagate_attr(param)

			# Prepare output
			loaded = [tmp if j == i else loaded[j] for j in range(len(loaded))]
			loaded = tuple(loaded)

		return(loaded)

	def propagate_attr(self, key):
		""" Propagate attribute current value to every condition. """
		for i in range(len(self.conds)):
			self.conds[i][key] = self[key]

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
			return(super(Main, self).__setattr__(name, value))
		else:
			# Don't set the attribute
			return(None)

	def check_attr(self, name, value):
		""" Check attribute format and value. """

		# Default answer
		checked = True

		if name in ['basedir', 'outdir']:
			# Require a string
			if not type('') == type(value):
				checked = False
			# # Require an existing folder
			# elif not os.path.isdir(value):
			# 	checked = False

			if not checked:
				msg = '"' + name + '" must be an existing folder\'s path.\n'
				msg += 'Keeping previous value.'
				self.printout(msg, -1)
		elif 'ncores' == name:
			# Require an integer
			if not type(0) == type(value):
				checked = False
			# Require a positive non-zero integer
			elif 0 >= value:
					checekd = False

			if not checked:
				msg = '"' + name + '" must be a positive non-zero integer.\n'
				msg += 'Keeping previous value.'
				self.printout(msg, -1)
		elif 'ext' == name:
			# Require a string
			if not type('') == type(value):
				checked = False
			# Require a file extension (\..+)
			elif not value.startswith('.') or 1 >= len(value):
				checked = False

			if not checked:
				msg = '"' + name + '" must be a file extension'
				msg += ' (start with ".").\n'
				msg += 'Keeping previous value.'
				self.printout(msg, -1)
		elif 'skip' == name:
			# Require a list
			if not type([]) == type(value):
				checked = False
			elif not 0 == len(value):
				# Require integer list
				types = [-1 for i in value if type(i) != type(0)]
				types = 0 == len(types)
				if not types:
					checked = False

			if not checked:
				msg = '"' + name + '" must be an integer list (can be empty).\n'
				msg += 'Keeping previous value.'
				self.printout(msg, -1)
		elif 'font_size' == name:
			# Require positive number
			if not type(value) in [type(0), type(0.)]:
				checked = False
			elif 0 >= value:
				checked = False

			if not checked:
				msg = '"' + name + '" must be a positive non-zero number.\n'
				msg += 'Keeping previous value.'
				self.printout(msg, -1)
		elif 'reg' == name:
			# Require a string
			if not type('') == type(value):
				checked = False

				msg = '"' + name + '" must be a regular expression (string).\n'
				msg += 'Keeping previous value.'
				self.printout(msg, -1)
		
		# Output
		return(checked)

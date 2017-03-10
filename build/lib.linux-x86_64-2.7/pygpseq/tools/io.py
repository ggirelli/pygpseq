# -*- coding: utf-8 -*-

""" Functions for the management of IO operations """

from datetime import datetime
import numpy as np
import os, sys, tempfile
from time import time

from .. import const

from . import path as pt, string as st

class IOinterface(object):
	"""Interface with input-output methods

	Attributes:
		logpath (string): path to log file.
		verbose (bool): True to be verbose.
	"""

	logpath = ''
	verbose = True

	def __init__(self, **kwargs):
		"""
		Args:
			path (string): path to the log file
			append (bool): whether to append to existing log (opt)
		"""

		super(IOinterface, self).__init__()

		# Append to existing log?
		if not 'append' in kwargs.keys():
			append = False
		else:
			append = kwargs['append']

		# Select path for log file
		if not 'path' in kwargs.keys():
			curpath = pt.add_trailing_slash(tempfile.gettempdir())
			curpath += 'pyGPSeq/log/'
			now = datetime.fromtimestamp(time()).strftime('%Y-%m-%d_%H:%M:%S')
			curpath += now + '.log'
		else:
			curpath = kwargs['path']

		# Create log file directory if missing
		self.check_log_dir(curpath)

		if not append:
			# Edit path if logfile exists already
			c = 1
			while os.path.exists(curpath):
				fname, fext = os.path.splitext(curpath)
				fname = fname + ' (' + str(c) + ')'
				curpath = fname + fext
				c += 1

		self.logpath = curpath

	def printout(self, s, lvl):
		"""Output for user.

		Args:
			s (string): message string.
			lvl (int): message level.

		Returns:
			string: formatted string.
		"""

		if self.verbose:
			s = printout(s, lvl)
			self.log(s)
		return(printout(s, lvl, verbose = False))

	def log(self, s):
		"""Write to logfile.

		Args:
			s (string): formatted message string.

		Returns:
			None: writes to logfile.
		"""

		# Create log file directory if missing
		self.check_log_dir()

		# Write to logfile
		f = open(self.logpath, 'a')
		f.write(s)
		f.close()
	
	def check_log_dir(self, path = None):
		"""Create log file directory if missing.

		Args:
			path (string): optional log path.

		Returns:
			None: creates the folder that will contain the log, if missing.
		"""

		if None == path:
			path = self.logpath

		dpath = os.path.dirname(path)
		if not os.path.isdir(dpath):
			os.makedirs(dpath)

	def gen_log_name(self):
		"""Generate logfile name. """
		now = datetime.fromtimestamp(time()).strftime('%Y-%m-%d_%H:%M:%S')
		name = now + '_pyGPSeq.log'
		return(name)

def printout(s, lvl, verbose = True):
	"""Log to shell.

	Args:
		s (string): message string.
		lvl (int): message level.
		verbose (bool): True to display formatted message.
	"""
	
	# Only a string can be logged
	if type(str()) != type(s):
		return()

	# Add level-based prefix
	if -2 == lvl:
		s = '\n~~ ERROR ~~ ლ(ಠ益ಠლ)\n' + s + '\nTerminated.\n'
		print(s)
		sys.exit()
	elif -1 == lvl:
		s = '\n~~ WARNING ~~ (ノ ゜Д゜)ノ ︵ ┻━┻\n' + s
	elif 0 == lvl:
		s = ' ' + s
	elif 1 == lvl:
		s = '  · ' + s
	elif 2 == lvl:
		s = '    > ' + s
	elif 3 == lvl:
		s = '    >> ' + s
	elif 4 <= lvl:
		s = '    >>> ' + s

	# Log
	if verbose:
		print(s)
	return(st.add_trailing_new_line(s))

def merge_profiles(profiles):
	"""Formats the profiles into a single table.

	Args:
		profiles (list): list of profile dictionaries.

	Returns:
		np.array: merged profiles.
	"""

	nsteps = len(profiles[0]['ratio']['x'])
	nrows = nsteps * len(profiles)

	out = np.zeros((nrows,), dtype = const.DTYPE_PROFILE_EXPORT)

	c = 0
	for profile in profiles:
		subcsv = np.zeros((nsteps,), dtype = const.DTYPE_PROFILE_EXPORT)
		subcsv['condition'] = [profile['condition'] for i in range(nsteps)]
		subcsv['x'] = profile['ratio']['x']
		subcsv['dna_mean'] = profile['dna']['mean']
		subcsv['dna_median'] = profile['dna']['median']
		subcsv['dna_mode'] = profile['dna']['mode']
		subcsv['dna_std'] = profile['dna']['std']
		subcsv['sig_mean'] = profile['sig']['mean']
		subcsv['sig_median'] = profile['sig']['median']
		subcsv['sig_mode'] = profile['sig']['mode']
		subcsv['sig_std'] = profile['sig']['std']
		subcsv['ratio_mean'] = profile['ratio']['mean']
		subcsv['ratio_median'] = profile['ratio']['median']
		subcsv['ratio_mode'] = profile['ratio']['mode']
		subcsv['ratio_std'] = profile['ratio']['std']
		subcsv['dna_mean_raw'] = profile['dna']['mean_raw']
		subcsv['dna_median_raw'] = profile['dna']['median_raw']
		subcsv['dna_mode_raw'] = profile['dna']['mode_raw']
		subcsv['dna_std_raw'] = profile['dna']['std_raw']
		subcsv['sig_mean_raw'] = profile['sig']['mean_raw']
		subcsv['sig_median_raw'] = profile['sig']['median_raw']
		subcsv['sig_mode_raw'] = profile['sig']['mode_raw']
		subcsv['sig_std_raw'] = profile['sig']['std_raw']
		subcsv['ratio_mean_raw'] = profile['ratio']['mean_raw']
		subcsv['ratio_median_raw'] = profile['ratio']['median_raw']
		subcsv['ratio_mode_raw'] = profile['ratio']['mode_raw']
		subcsv['ratio_std_raw'] = profile['ratio']['std_raw']
		subcsv['n'] = [profile['n'] for i in range(nsteps)]
		out[c:(c+nsteps),] = subcsv
		c += nsteps

	return(out)

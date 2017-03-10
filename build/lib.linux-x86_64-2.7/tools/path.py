# -*- coding: utf-8 -*-

""" Functions for the management of paths """

import os
import re

from .. import const

def add_extension(path, ext):
	""" If the provided path is missing the extension, add it. """

	# Add lading dot to extension, if needed
	ext = add_leading_dot(ext)

	# Check presence of the extension
	if not path.endswith(ext):
		path += ext

	# Output
	return(path)

def add_trailing_slash(path):
	""" Add the trailing slash to a path """
	if '/' != path[len(path) - 1]:
		path = path + '/'
	return(path)

def add_leading_dot(ext):
	""" Add the leading dot to an extension """
	if '.' != ext[0]:
		ext = '.' + ext
	return(ext)

def select_folders(path, ext):
	""" Select subdirectories containing files with the given extension """

	# Check input params
	path = add_trailing_slash(os.path.abspath(path))
	ext = add_leading_dot(ext)

	# Return an empty list if the provided path is not a directory
	if not os.path.isdir(path):
		return([])
	else:
		# Retreive subdirectory list
		sdirs = [add_trailing_slash(path + x) for x in next(os.walk(path))[1]]

		# Will contain the selected subdirectories
		selected = []

		for sdir in sdirs:

			# Retreive file list
			flist = os.listdir(sdir)

			# Look for files with the proper extension
			if 0 != len(select_files(sdir, ext)):
				selected.append(sdir)

		selected.sort()

		return(selected)

def select_files(path, ext):
	""" Select the files with the proper .ext in the provided path """
	
	# Check input params
	path = add_trailing_slash(os.path.abspath(path))
	ext = add_leading_dot(ext)

	# Return an empty list if the provided path is not a directory
	if not os.path.isdir(path):
		return([])
	else:
		# Retreive file list
		flist = os.listdir(path)

		# Will contain the selected files
		selected = []

		for f in flist:
			if os.path.isfile(path + f):
				# Retreive file extension
				fname, fext = os.path.splitext(f)

				# Compare file extension with .ext
				if fext == ext:
					selected.append(f)

		selected.sort()

		return(selected)

def select_series(flist, reg, series_field = None):
	""" Group a list of files by series based on the provided regexp """
	
	# Default series string field name
	if None == series_field:
		series_field = const.REG_SERIES_ID

	if type(str()) == type(reg):
		# Compile regexp
		reg = re.compile(reg)

	# Output dictionary
	d = {}

	# Select files matching the regexp
	for f in flist:

		# Search for matches
		m = re.search(reg, f)
		
		if None != m:
			# Identify the series string
			scur = m.groupdict()[series_field]

			# Save in the output dictionary
			if scur in d.keys():
				# Add channel to existing series
				d[scur][f] = m.groupdict()
			else:
				# Add new series
				d[scur] = {f : m.groupdict()}
	
	return(d)

#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
./run_pyGPSeq_from_table.py paramTable i

Runs pyGPSeq by reading the required parameter from the paramTable row i.
"""

# Load library
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import pygpseq as gp
import sys


# CHECK PARAMETERS
# =================================

# Check number of provided parameters
if 3 != len(sys.argv):
	sys.exit('Proper usage: ./run_pyGPSeq_from_table paramTable i.\nSTOPPED')

# Check that paramTable exists
paramTable = sys.argv[1]
if not os.path.isfile(paramTable):
	msg = 'Proper usage: ./run_pyGPSeq_from_table paramTable i.\n'
	msg += 'paramTable not found at "' + paramTable + '".\nSTOPPED'
	sys.exit(msg)

# Read paramTable
t = pd.read_table(paramTable, delimiter = '\t', dtype = 'string')

# Check that i is a number
i = sys.argv[2]
if not all([d.isdigit() for d in i]):
	msg = 'Proper usage: ./run_pyGPSeq_from_table paramTable i.\n'
	msg += 'i should be an integer.\nSTOPPED'
	sys.exit(msg)
i = int(i)

# Check that i is a row of paramTable
if not i in range(t.shape[0] + 1):
	msg = 'Proper usage: ./run_pyGPSeq_from_table paramTable i.\n'
	msg += 'Cannot find row i [' + str(i) + '] in paramTable [max ' + str(t.shape[0]) + '].\nSTOPPED'
	sys.exit(msg)
i -= 1

# FUNCTIONS
# ==============================

def runParam(t, i):
	if 'ncores' in t.columns:
		gpi = gp.Main(ncores = int(t['ncores'][i]))
	else:
		gpi = gp.Main()

	# Select only proper parameters
	proper = lambda pname: pname in t.columns and not '' == t[pname][i]
	plist = [pname for pname in dir(gpi) if proper(pname)]

	# Store parameters
	for pname in plist:
		if not type('') == type(gpi[pname]):
			gpi[pname] = eval(t[pname][i])
		else:
			gpi[pname] = t[pname][i]

	gpi.reg = '^(?P<' + gp.const.REG_CHANNEL_NAME + '>[^/]*)'
	gpi.reg += '\.(?P<' + gp.const.REG_CHANNEL_ID + '>channel[0-9]+)'
	gpi.reg += '\.(?P<' + gp.const.REG_SERIES_ID + '>series[0-9]+)'
	gpi.reg += '(?P<' + gp.const.REG_EXT + '>(_cmle)?\.tif)$'

	gpi.run()

	plt.close('all')

	return(None)

# RUN
# ==============================

if 'todo' in t.columns:
	if t['todo'][i] in ['False', 'True']:
		if eval(t['todo'][i]):
			runParam(t, i)
		else:
			sys.exit('Skipped row ' + str(i+1) + ' [todo = False].')
	else:
		sys.exit('Skipped row ' + str(i+1) + ' due to unrecognized todo value [todo = ' + t['todo'][i] + '].')
else:
	runParam(t, i)
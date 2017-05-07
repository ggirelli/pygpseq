#!/usr/bin/python
# -*- coding: utf-8 -*-

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 1.0.0
# Description: script to interface with the pygpseq package.
# ------------------------------------------------------------------------------



# DEPENDENCIES =================================================================

import argparse
import numpy as np
import os
import pandas as pd
import pygpseq as gp
import sys

# PARAMETERS ===================================================================

# Add script description
parser = argparse.ArgumentParser(
	description = """Run pygpseq-based analysis."""
)

'''

Mandatory
- basedir
- outdir
- aspect 300,216.6,216.6

Default
- dna_names dapi
- sig_names cy5,tmr
- seg_type 1|2|3
- an_type 1|2|3|4
- min_z_size .25
- nuclear_selection 1|2|3|4|5|6
- logpath
- rescale_deconvolved
- normalize_distance

Optional
- ncores 1
- skip 1,2,3,4
- description_file
- notes

Advanced
- regexp

'''

# Add params
parser.add_argument('inDir', type = str, nargs = 1,
	help = """Path to input directory, containing single-condition directories
	with TIF files.""")

parser.add_argument('chrlen', metavar = 'chrlen', type = str, nargs = 1,
	help = 'Path to file with chromosome lengths (chr, length)' +
	' or chromosome length.')
parser.add_argument('chr', metavar = 'chr', type = str, nargs = 1,
	help = 'The chromosome to bin. E.g., chr1')

# Add flags
parser.add_argument('--binsize', metavar = 'bsi', type = int, nargs = 1,
	default = [1e6],
	help = 'Bin size. Default: 1e6')
parser.add_argument('--binstep', metavar = 'bst', type = int, nargs = 1,
	default = [1e6],
	help = 'Bin step. Non-overlapping bins if equal to bin size. Default: 1e5')

# Parse arguments
args = parser.parse_args()

# Retrieve arguments
chrfile = args.chrlen[0]
schr = args.chr[0]
size = int(args.binsize[0])
step = int(args.binstep[0])


# RUN ==========================================================================

# End --------------------------------------------------------------------------

################################################################################

# Import library

# Create pyGPSeq analyzer instance
gpi = gp.Main(ncores = 10)

# Steps to be skipped
# 1	  : skip instantiation (need gpi.inst.cpickle, otherwise unskips)
# 2	  : skip segmentation (need gpi.seg.cpickle, otherwise unskips)
# 3	  : skip analysis (need gpi.an.cpickle, otherwise unskips)
# 3.5 : skip single-nuclei boxplot (end of step 3)
# 4	  : final plots
# 5   : final report
gpi.skip = [1,2,3,4]

gpi.sig_names = ('dapi', 'cy5', 'tmr')
gpi.dna_names = ('dapi',)

# Data directory
gpi.basedir = '/home/gire/Desktop/BiCro-Data/Imaging/iGG043_050_deco'

# Output directory
gpi.outdir = '/home/gire/Desktop/BiCro-Analysis/Imaging/iGG043_050_deco'

# Segmentation type for nuclear identification
# SEG_SUM_PROJ	: through sum Z projection
# SEG_MAX_PROJ	: through max Z projection
# SEG_3D		: 3D segmentation
gpi.seg_type = gp.const.SEG_3D

# Single-nucleus analysis type
# AN_SUM_PROJ	: on nucleus sum Z projection
# AN_MAX_PROJ	: on nucleus max Z projection
# AN_3D			: on whole and partial nuclear volume
# AN_MID		: only on nucleus mid-section (highest sum of DNA staining intensity)
gpi.an_type = gp.const.AN_MID

# Voxel aspect proportions (or sizes, ZYX)
gpi.aspect = (300., 216.6, 216.6)

# Minimum percentage of stack to be occupied by a cell
gpi.min_z_size = .25

# Nuclear selection
# NSEL_SIZE			: nuclear size (either volume or area)
# NSEL_SURF			: nuclear surface (only if AN_3D)
# NSEL_SHAPE		: nuclear shape (either sphericity or circularity)
# NSEL_SUMI			: nuclear DNA intensity sum
# NSEL_MEANI		: nuclear DNA intensity average
# NSEL_FLAT_SIZE	: nuclear Z-projected area
gpi.nsf = (gp.const.NSEL_FLAT_SIZE, gp.const.NSEL_SUMI)

# Regular expression to identify image files
gpi.reg = '^(?P<' + gp.const.REG_CHANNEL_NAME + '>[^/]*)'
gpi.reg += '\.(?P<' + gp.const.REG_CHANNEL_ID + '>channel[0-9]+)'
gpi.reg += '\.(?P<' + gp.const.REG_SERIES_ID + '>series[0-9]+)'
gpi.reg += '(?P<' + gp.const.REG_EXT + '>(_cmle)?\.tif)$'

# Where to save the run log
gpi.logpath = '/home/gire/Desktop/BiCro-Analysis/Imaging/iGG043_050_deco/log'
#gpi.logpath += gpi.gen_log_name()

# Perform deconvolved image rescaling?
gpi.rescale_deconvolved = True

# Normalize distance?
gpi.normalize_distance = True

# Better condition naming
gpi.cdescr['iGG043_20170411_001'] = 'DAPI 1 ug/ml'
gpi.cdescr['iGG044_20170411_001'] = 'DAPI 200 ng/ml'
gpi.cdescr['iGG045_20170412_001'] = 'DAPI 50 ng/ml'
gpi.cdescr['iGG046_20170412_001'] = 'DAPI 20 ng/ml'
gpi.cdescr['iGG047_20170412_001'] = 'Hoechst 1 ug/ml'
gpi.cdescr['iGG048_20170412_001'] = 'Hoechst 200 ng/ml'
gpi.cdescr['iGG049_20170412_001'] = 'Hoechst 50 ng/ml'
gpi.cdescr['iGG050_20170412_001'] = 'Hoechst 20 ng/ml'

# Notes
gpi.notes = 'DNA staining test on HeLa P20/21.'

# Start the analysis
gpi = gpi.run()


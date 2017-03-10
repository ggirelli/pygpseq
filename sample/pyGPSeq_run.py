#!/usr/bin/python
# -*- coding: utf-8 -*-

# Import library
import pygpseq as gp

# Create pyGPSeq analyzer instance
gpi = gp.Main(ncores = 5)

# Steps to be skipped
# 1	  : skip instantiation (need gpi.inst.cpickle, otherwise unskips)
# 2	  : skip segmentation (need gpi.seg.cpickle, otherwise unskips)
# 3	  : skip analysis (need gpi.an.cpickle, otherwise unskips)
# 3.5 : skip single-nuclei boxplot (end of step 3)
# 4	  : final plots
# 5   : final report
gpi.skip = []

# Data directory
gpi.basedir = '/home/gire/Desktop/BiCro-Data/Imaging/iTK_61_66_deconvolved'

# Output directory
gpi.outdir = '/home/gire/Desktop/BiCro-Analysis/Imaging/iTK_61_66_deconvolved'

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
gpi.logpath = '/home/gire/Desktop/BiCro-Analysis/Imaging/iTK_61_66_deconvolved/log'
#gpi.logpath += gpi.gen_log_name()

# Perform deconvolved image rescaling?
gpi.rescale_deconvolved = True

# Normalize distance?
gpi.normalize_distance = True

# Better condition naming
gpi.cdescr['iTK61_031116_001'] = '0N'
gpi.cdescr['iTK62_031116_001'] = '30min'
gpi.cdescr['iTK63_031116_001'] = '15min'
gpi.cdescr['iTK64_031116_001'] = '5min'
gpi.cdescr['iTK65_031116_001'] = '1min'
gpi.cdescr['iTK66_031116_001'] = 'negative'

# Notes
gpi.notes = 'HeLa cells, 200U HindIII, super long Y-FISH linker'

# Start the analysis
gpi = gpi.run()

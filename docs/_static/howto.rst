How to
======

Set up your data for pyGPSeq
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**pyGPSeq** requires your data (images) to be already in `.tif` format. The following steps include also convertion of your images to `.tif` format using the `toTif_fiji.js` script.

1. First of all, create an empty directory which will contain your data. We will refer to this diretory as **basedir**.
2. Your `.tif` images should be moved to subfolders of **basedir**. Note that each subfolder in **basedir** is expected to be a different condition. Thus, group in the same subfolder images from the same condition, and divide in different subfolders images from different conditions.
    * If your images are not in `.tif` format but in a microscope-specific format, you can use the `toTif_fiji.js` script. The script will convert your images to `.tif` format (one image per channel), moving them in a subfolder with the same name as the initial image.
    * Given how the script works, it is convenient to have one microscope-specific formatted image per condition. If that is the case, one subfolder will be created for every condition, containing the `.tif` images in a proper format. These subfolders can then be moved to **basedir** to proceed with the analysis.

Run pyGPSeq
~~~~~~~~~~~

**pyGPSeq** can be run on a single dataset or on multiple datasets (batch) at once. In the latter case, a queue is created and only one dataset at a time is analyzed. The following is an example of a script used to run the pipeline, where most of the available parameters are reported with a description.

```

	#!/usr/bin/python
	# -*- coding: utf-8 -*-

	# Import library
	import pygpseq as gp

	# Create pyGPSeq analyzer instance
	# ncores specifies the number of threads to be used for parallelization!!! (!IMPORTANT!)
	gpi = gp.Main(ncores = 10)

	# Steps to be skipped
	# 1	  : skip instantiation (need gpi.inst.cpickle, otherwise unskips)
	# 2	  : skip segmentation (need gpi.seg.cpickle, otherwise unskips)
	# 3	  : skip analysis (need gpi.an.cpickle, otherwise unskips)
	# 3.5 : skip single-nuclei boxplot (end of step 3)
	# 4	  : final plots
	# 5   : final report
	gpi.skip = [1,2,3,4]

	# Name of the signal channels, a final comma is needed if only a channel is available
	gpi.sig_names = ('cy5', 'tmr')

	# Name of the DNA staining channels, a final comma is needed if only a channel is available
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

	# Better condition naming (convert folder name to condition description)
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
```

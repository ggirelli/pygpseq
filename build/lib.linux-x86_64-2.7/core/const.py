# -*- coding: utf-8 -*-

import pkg_resources

# Taken from:
# http://code.activestate.com/recipes/65207-constants-in-python/?in=user-97991
class _const:
	class ConstError(TypeError): pass
	def __setattr__(self, name, value):
		if name in self.__dict__.keys():
			raise(self.ConstError, "Can't rebind const(%s)" % name)
		self.__dict__[name] = value

# Package version
#_const.VERSION = pkg_resources.get_distribution('pyGPSeq').version
_const.VERSION = '0.1.0'
_const.PACK_NAME = 'pygpseq'

# Series regexp fields
_const.REG_PATH = 'abs_path'
_const.REG_CHANNEL_NAME = 'channel_name'
_const.REG_CHANNEL_ID = 'channel_id'
_const.REG_SERIES_ID = 'series_id'
_const.REG_EXT = 'ext'

# Steps
_const.STEP_DESCR = ('Instantiation', 'Segmentation', 'Analysis',
	'General boxplot', 'Final plots', 'Final report')

# Projection types
_const.SUM_PROJ = 0
_const.MAX_PROJ = 1

# Segmentation types
_const.SEG_SUM_PROJ = _const.SUM_PROJ
_const.SEG_MAX_PROJ = _const.MAX_PROJ
_const.SEG_3D = 2
_const.SEG_LABELS = ('Sum Z projection', 'Max Z projection', '3D')

# Analysis types
_const.AN_SUM_PROJ = 0
_const.AN_MAX_PROJ = 1
_const.AN_3D = 2
_const.AN_MID = 3
_const.AN_LABELS = ('Sum Z projection', 'Max Z projection', '3D', 'Mid-section')

# Colors
_const.PALETTE = ('#a6cee3', '#1f78b4', '#b2df8a', '#33a02c', '#fb9a99',
	'#e31a1c', '#fdbf6f', '#ff7f00', '#cab2d6', '#6a3d9a', '#ffff99', '#b15928')

# DTYPE export
_const.DTYPE_NUCLEAR_SUMMARY = [('s', 'u4'), ('n', 'u4'),
	('flat_size', 'u8'), ('size', 'u8'),
	('surf', 'f8'), ('sumI', 'f8'), ('meanI', 'f8'), ('shape', 'f8')]
_const.DTYPE_NUCLEAR_DATA = [('s', 'u4'), ('n', 'u4'), ('dna', 'u8'),
	('sig', 'u8'), ('d', 'f8'), ('dnorm', 'f8'), ('part', 'b')]
_const.DTYPE_PROFILE_EXPORT = [
	('condition', 'S100'), ('x', 'f'),
	('dna_mean', 'f'), ('sig_mean', 'f'), ('ratio_mean', 'f'),
	('dna_mean_raw', 'f'), ('sig_mean_raw', 'f'), ('ratio_mean_raw', 'f'),
	('dna_median', 'f'), ('sig_median', 'f'), ('ratio_median', 'f'),
	('dna_median_raw', 'f'), ('sig_median_raw', 'f'), ('ratio_median_raw', 'f'),
	('dna_mode', 'f'), ('sig_mode', 'f'), ('ratio_mode', 'f'),
	('dna_mode_raw', 'f'), ('sig_mode_raw', 'f'), ('ratio_mode_raw', 'f'),
	('dna_std', 'f'), ('sig_std', 'f'), ('ratio_std', 'f'),
	('dna_std_raw', 'f'), ('sig_std_raw', 'f'), ('ratio_std_raw', 'f'),
	('n', 'u4')]

# Nuclear selection features
_const.NSEL_SIZE = 0
_const.NSEL_SURF = 1
_const.NSEL_SHAPE = 2
_const.NSEL_SUMI = 3
_const.NSEL_MEANI = 4
_const.NSEL_FLAT_SIZE = 5
_const.NSEL_FIELDS = ('size', 'surf', 'shape', 'sumI', 'meanI', 'flat_size')
_const.NSEL_NAMES = ('Size', 'Surface', 'Shape', 'Intensity Sum',
	'Mean Intensity', 'Area')
_const.NSEL_LABELS = ('auto', 'Surface [a.u.]', 'auto',
	'Intensity Sum [a.u.]', 'Mean Intensity [a.u.]', 'Area [px]')

# Output folders (MUST have trailing slash)
_const.OUTDIR_PDF = 'out_pdf/'
_const.OUTDIR_PNG = 'out_png/'
_const.OUTDIR_PNG_REPORT = _const.OUTDIR_PNG + 'report/'
_const.OUTDIR_TIF = 'out_tif/'
_const.OUTDIR_MASK = 'out_masks/'
_const.OUTDIR_DEBUG = 'debugging/'

# Plot constants
_const.SCI_FORMAT = '%2.e'

# kwargs automatic update
_const.KWARGS_TYPELIST = (type([]), type(()), type(0), type(''),
	type(True), type({}), type(0.0))
_const.KWARGS_AVOIDLIST = ('conds')

# Step-related main() class parameters
_const.PARAM_STATIC = ('basedir', 'cdescr', 'debugging', 'font_size', 'logpath',
	'ncores', 'notes', 'outdir', 'plotting', 'skip', 'suffix', 'verbose')
_const.PARAM_SEG = ('adp_thr', 'calc_n_surface', 'dna_names', 'ext',
	'min_z_size', 'seg_type', 'sig_names', 'offset', 'radius_interval', 'reg',
	'rescale_deconvolved', 'rm_z_tips', 'seg_type', 'sig_names', 'sigma')
_const.PARAM_AN = ('an_type', 'aspect', 'nbins', 'normalize_distance', 'nsf',
	'part_n_erosion')
_const.PARAM_PROPAGATE = ('logpath',)

# Save constants
import sys
sys.modules[__name__] = _const()

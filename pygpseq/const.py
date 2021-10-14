# -*- coding: utf-8 -*-

# ------------------------------------------------------------------------------
#
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Description: keeps constant values for the pygpseq package.
#
# ------------------------------------------------------------------------------

# DEPENDENCIES =================================================================

import pkg_resources

# FUNCTIONS ====================================================================

# Taken from:
# http://code.activestate.com/recipes/65207-constants-in-python/?in=user-97991
class _const:
    class ConstError(TypeError):
        pass

    def __setattr__(self, name, value):
        if name in self.__dict__.keys():
            raise (self.ConstError, "Can't rebind const(%s)" % name)
        self.__dict__[name] = value


# CONSTANTS ====================================================================

# Package-related --------------------------------------------------------------

# Package version
_const.VERSION = pkg_resources.get_distribution("pygpseq").version
_const.PACK_NAME = "pygpseq"

# Parameter-related ------------------------------------------------------------

# kwargs automatic update
_const.KWARGS_TYPELIST = (
    type([]),
    type(()),
    type(0),
    type(""),
    type(True),
    type({}),
    type(0.0),
    type(None),
)
_const.KWARGS_AVOIDLIST = "conds"

# Step-related main() class parameters
_const.PARAM_STATIC = (
    "basedir",
    "cdescr",
    "debugging",
    "font_size",
    "logpath",
    "ncores",
    "notes",
    "outdir",
    "plotting",
    "skip",
    "suffix",
    "verbose",
)
_const.PARAM_SEG = (
    "adp_thr",
    "calc_n_surface",
    "dna_names",
    "ext",
    "min_z_size",
    "seg_type",
    "sig_names",
    "offset",
    "radius_interval",
    "reg",
    "rescale_deconvolved",
    "rm_z_tips",
    "seg_type",
    "sig_names",
)
_const.PARAM_AN = (
    "an_type",
    "aspect",
    "nbins",
    "normalize_distance",
    "nsf",
    "part_n_erosion",
    "sigma_smooth",
    "sigma_density",
)
_const.PARAM_PROPAGATE = ("logpath",)

# Analysis-related -------------------------------------------------------------

# Series regexp fields
_const.REG_PATH = "abs_path"
_const.REG_CHANNEL_NAME = "channel_name"
_const.REG_CHANNEL_ID = "channel_id"
_const.REG_SERIES_ID = "series_id"
_const.REG_EXT = "ext"

# Steps
_const.STEP_DESCR = (
    "Instantiation",
    "Segmentation",
    "Analysis",
    "General boxplot",
    "Final plots",
    "Final report",
)

# Projection types
_const.SUM_PROJ = 0
_const.MAX_PROJ = 1

# Segmentation types
_const.SEG_SUM_PROJ = _const.SUM_PROJ
_const.SEG_MAX_PROJ = _const.MAX_PROJ
_const.SEG_3D = 2
_const.SEG_DEFAULT = _const.SEG_3D
_const.SEG_LABELS = ("Sum Z projection", "Max Z projection", "3D")
_const.SEG_ARG_LABELS = ("sum_proj", "max_proj", "3d")

# Analysis types
_const.AN_SUM_PROJ = 0
_const.AN_MAX_PROJ = 1
_const.AN_3D = 2
_const.AN_MID = 3
_const.AN_DEFAULT = _const.AN_MID
_const.AN_LABELS = ("Sum Z projection", "Max Z projection", "3D", "Mid-section")
_const.AN_ARG_LABELS = ("sum_proj", "max_proj", "3d", "mid")

# Mid-section selection mode
_const.MID_SEC_CENTRAL = 0
_const.MID_SEC_LARGEST = 1
_const.MID_SEC_MAXSUMI = 2
_const.MID_SEC_DEFAULT = _const.MID_SEC_LARGEST
_const.MID_SEC_LABELS = ("central", "largest", "max intensity sum")
_const.MID_SEC_ARG_LABELS = ("central", "largest", "maxIsum")

# Lamina distance mode
_const.LD_CENTER_MAX = 0
_const.LD_CENTER_PERC = 1
_const.LD_DIFFUSION = 2
_const.LD_DEFAULT = _const.LD_CENTER_PERC
_const.LD_ARG_LABELS = ("center_max", "center_percentile", "diffusion")

# Nuclear selection features
_const.NSEL_SIZE = 0
_const.NSEL_SURF = 1
_const.NSEL_SHAPE = 2
_const.NSEL_SUMI = 3
_const.NSEL_MEANI = 4
_const.NSEL_FLAT_SIZE = 5
_const.NSEL_FLAT_SUMI = 5
_const.NSEL_FIELDS = (
    "size",
    "surf",
    "shape",
    "sumI",
    "meanI",
    "flat_size",
    "flat_sumI",
)
_const.NSEL_NAMES = (
    "Size",
    "Surface",
    "Shape",
    "Intensity Sum",
    "Mean Intensity",
    "Area",
    "Intensity Sum in Sum projection",
)
_const.NSEL_LABELS = (
    "auto",
    "Surface [a.u.]",
    "auto",
    "Intensity Sum [a.u.]",
    "Mean Intensity [a.u.]",
    "Area [px]",
    "Intensity Sum in Sum projection [a.u.]",
)

# Output-related ---------------------------------------------------------------

_const.DLAMIN_LABEL = "lamin_d"
_const.DLAMIN_NORM_LABEL = "lamin_dnorm"

# Nuclear summary
_const.DTYPE_NUCLEAR_SUMMARY_3D = [
    ("s", "u4"),
    ("n", "u4"),
    ("flat_size", "u8"),
    ("size", "u8"),
    ("surf", "f8"),
    ("sumI", "f8"),
    ("flat_sumI", "f8"),
    ("meanI", "f8"),
    ("shape", "f8"),
    ("box_start_slice", "u4"),
    ("box_start_row", "u4"),
    ("box_start_col", "u4"),
    ("box_end_slice", "u4"),
    ("box_end_row", "u4"),
    ("box_end_col", "u4"),
    ("com_slice", "f8"),
    ("com_row", "f8"),
    ("com_col", "f8"),
]
_const.DTYPE_NUCLEAR_SUMMARY_2D = [
    ("s", "u4"),
    ("n", "u4"),
    ("flat_size", "u8"),
    ("size", "u8"),
    ("surf", "f8"),
    ("sumI", "f8"),
    ("flat_sumI", "f8"),
    ("meanI", "f8"),
    ("shape", "f8"),
    ("box_start_row", "u4"),
    ("box_start_col", "u4"),
    ("box_end_row", "u4"),
    ("box_end_col", "u4"),
    ("com_row", "f8"),
    ("com_col", "f8"),
]

# Single pixel nuclear data
_const.DTYPE_NUCLEAR_DATA = [
    ("s", "u4"),
    ("n", "u4"),
    ("dna", "u8"),
    ("sig", "u8"),
    (_const.DLAMIN_LABEL, "f8"),
    ("centr_d", "f8"),
    (_const.DLAMIN_NORM_LABEL, "f8"),
    ("part", "b"),
]

# Nuclear data export
_const.DTYPE_NDATA_EXPORT_2D = [("c", "u4"), *_const.DTYPE_NUCLEAR_SUMMARY_2D]
_const.DTYPE_NDATA_EXPORT_3D = [("c", "u4"), *_const.DTYPE_NUCLEAR_SUMMARY_3D]
# Profile table
_const.PROFILE_TYPES = ("mean", "median", "mode", "std", "max")
_const.DATA_FOCUS = ("dna", "sig", "ratio")
_const.DTYPE_PROFILE_EXPORT = [("condition", "S100"), ("x", "f"), ("n", "u4")]
for df in _const.DATA_FOCUS:
    for pt in _const.PROFILE_TYPES:
        _const.DTYPE_PROFILE_EXPORT.append(("%s_%s" % (df, pt), "f"))
        _const.DTYPE_PROFILE_EXPORT.append(("%s_%s_raw" % (df, pt), "f"))

# Output folders (MUST have trailing slash)
_const.OUTDIR_PDF = "out_pdf/"
_const.OUTDIR_PNG = "out_png/"
_const.OUTDIR_CSV = "out_csv/"
_const.OUTDIR_MASK = "out_masks/"
_const.OUTDIR_PNG_REPORT = _const.OUTDIR_PNG + "report/"
_const.OUTDIR_TIF = "out_tif/"
_const.OUTDIR_DEBUG = "debugging/"

# Plot-related------------------------------------------------------------------

# Plot constants
_const.SCI_FORMAT = "%2.e"

# RUN ==========================================================================

# Save constants
import sys

sys.modules[__name__] = _const()

# END ==========================================================================

################################################################################

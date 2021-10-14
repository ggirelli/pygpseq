# ------------------------------------------------------------------------------
#
# MIT License
#
# Copyright (c) 2017 Gabriele Girelli
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
#
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Project: GPSeq
# Description: script to interface with the pygpseq.anim package.
#
# ------------------------------------------------------------------------------


# DEPENDENCIES =================================================================

import argparse
import datetime
import multiprocessing
import numpy as np
import os
import pandas as pd
import pygpseq as gp
import sys

from ggc.prompt import ask
from ggc.args import check_threads, export_settings

from pygpseq.tools.io import printout


# FUNCTION =====================================================================
version = "2.1.1"


def print_settings(gpi, args, readable_cdescr, readable_nsf, clear=True):
    """Show input settings, for confirmation.

    Args:
            args (Namespace): arguments parsed by argparse.
            clear (bool): clear screen before printing.
    """
    s = " # GPSeq image analysis v%s\n" % version

    s += """
---------- SETTING:  VALUE ----------

Input directory :  %s
Output directory :  %s
        Log file :  %s

    Skipped steps :  %s

    DNA channels :  %s
Signal channels :  %s

    Segmentation :  %s
    Mask folder :  %s
    Mask prefix :  %s
        Labeled :  %r
        Compressed :  %r

        Analysis :  %s
    Middle section :  %s
    Distance mode :  %s

Voxel aspect (ZYX) :  %s
    Aspect unit :  %s
Minimum Z portion :  %.2f
    Minimum radius :  %.2f vx
        Fill holes :  %r

    Sigma (smooth) :  %.4f
Sigma (density) :  %.4f
            #bins :  %d

Condition descr. : %s

Nuclear selection :  %s

        Threads :  %d
            Note :  %s

            Regexp :  '%s'

    Correct shift :  %r
Rescale deconv. :  %r
Normalize dist. :  %r
        Debug mod :  %r

    """ % (
        gpi.basedir,
        gpi.outdir,
        gpi.logpath,
        str(args.skip),
        str(gpi.dna_names),
        str(gpi.sig_names),
        args.seg_type,
        args.mask_folder,
        args.mask_prefix,
        args.labeled,
        args.compressed,
        args.an_type,
        args.mid_type,
        args.dist_type,
        str(gpi.aspect),
        gpi.umes,
        gpi.min_z_size,
        gpi.radius_interval[0],
        gpi.do_fill_holes,
        gpi.sigma_smooth,
        gpi.sigma_density,
        gpi.nbins,
        "\n                     ".join(readable_cdescr),
        readable_nsf,
        gpi.ncores,
        gpi.notes,
        gpi.reg,
        gpi.correct_shift,
        gpi.rescale_deconvolved,
        gpi.normalize_distance,
        gpi.debugging,
    )

    if clear:
        print("\033[H\033[J%s" % s)
    else:
        print(s)
    return s


def run():
    # PARAMETERS ===================================================================

    # Add script description
    parser = argparse.ArgumentParser(
        description="""
	Analyse a multi-condition GPSeq image dataset.
	""",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # Positional parameters
    parser.add_argument(
        "inDir",
        type=str,
        help="""Path to input directory, containing single-condition directories
		with TIF files.""",
    )
    parser.add_argument(
        "outDir",
        type=str,
        help="""Path to output directory, must be different from the input
		directory.""",
    )

    # Optional parameters
    parser.add_argument(
        "--skip",
        type=str,
        nargs="*",
        help="""Space-separated phases to be skipped.
		Use -- after the last one.""",
        choices=["inst", "seg", "an", "box", "plot", "report"],
    )
    parser.add_argument(
        "-l",
        "--logpath",
        type=str,
        help="""Path to log file. By default: outDir/log""",
        metavar="log",
    )
    parser.add_argument(
        "-a",
        "--aspect",
        type=float,
        nargs=3,
        help="""Physical size of Z, Y and X voxel sides.
		Default: 300.0 216.6 216.6""",
        metavar=("Z", "Y", "X"),
        default=[300.0, 216.6, 216.6],
    )
    parser.add_argument(
        "-U",
        "--umes",
        type=str,
        help="""Unit of measure for the aspect. Default: nm""",
        metavar="unit",
        default="nm",
    )
    parser.add_argument(
        "-d",
        "--dna-channels",
        type=str,
        nargs="+",
        help="""Space-separated names of DNA staining channels.
		Use -- after the last one.""",
        default=["dapi"],
        metavar="dna_name",
    )
    parser.add_argument(
        "-s",
        "--sig-channels",
        type=str,
        nargs="+",
        help="""Space-separated names of GPSeq signal channels.
		Use -- after the last one.""",
        default=["cy5", "tmr"],
        metavar="sig_name",
    )
    parser.add_argument(
        "-z",
        "--min-z",
        type=float,
        help="""If lower than 1, minimum fraction of stack, if higher than 1,
		minimum number of slices to be occupied by a nucleus. Default: .25""",
        default=0.25,
        metavar="min_z",
    )
    parser.add_argument(
        "-R",
        "--min-radius",
        type=float,
        help="""Minimum radius (in voxels) for segmented objects. Default: 10""",
        default=10.0,
        metavar="min_radius",
    )
    parser.add_argument(
        "--seg-type",
        type=str,
        help="""Segmentation type. Default: '%s'"""
        % (gp.const.SEG_ARG_LABELS[gp.const.SEG_DEFAULT]),
        choices=list(gp.const.SEG_ARG_LABELS),
        default=gp.const.SEG_ARG_LABELS[gp.const.SEG_DEFAULT],
    )
    parser.add_argument(
        "-m",
        "--mask-folder",
        metavar="folder",
        type=str,
        help="""Path to folder containing binarized/labeled images.
		Masks will be saved to this folder if missing.""",
        default=None,
    )
    parser.add_argument(
        "-M",
        "--mask-prefix",
        metavar="prefix",
        type=str,
        help="""Prefix for mask selection. Default: 'mask_'.""",
        default="mask_",
    )
    parser.add_argument(
        "-2",
        "--manual-2d-masks",
        type=str,
        metavar="MAN2DDIR",
        help="""Path to folder with 2D masks with matching name,
		to combine with 3D masks. Used only if no mask was provided through -m.""",
    )
    parser.add_argument(
        "--an-type",
        type=str,
        help="""Analysis type. Default: '%s'"""
        % (gp.const.AN_ARG_LABELS[gp.const.AN_DEFAULT]),
        choices=list(gp.const.AN_ARG_LABELS),
        default=gp.const.AN_ARG_LABELS[gp.const.AN_DEFAULT],
    )
    parser.add_argument(
        "--mid-type",
        type=str,
        help="""Method for mid-section selection. Default: '%s'"""
        % (gp.const.MID_SEC_ARG_LABELS[gp.const.MID_SEC_DEFAULT]),
        choices=list(gp.const.MID_SEC_ARG_LABELS),
        default=gp.const.MID_SEC_ARG_LABELS[gp.const.MID_SEC_DEFAULT],
    )
    parser.add_argument(
        "--dist-type",
        type=str,
        help="""Method for lamina distance calculation/normalization.
		Default: '%s'"""
        % (gp.const.LD_ARG_LABELS[gp.const.LD_DEFAULT]),
        choices=list(gp.const.LD_ARG_LABELS),
        default=gp.const.LD_ARG_LABELS[gp.const.LD_DEFAULT],
    )
    parser.add_argument(
        "--nuclear-sel",
        type=str,
        nargs="*",
        help="""Space-separated features for nuclear selection.
		Use -- after the last one. Default: flat_size sumI""",
        choices=list(gp.const.NSEL_FIELDS),
        default=["flat_size", "sumI"],
    )
    parser.add_argument(
        "--sigma-smooth",
        metavar="sigmaValue",
        type=float,
        default=0.1,
        help="""Sigma value for sparse gaussian smoothing of profiles.""",
    )
    parser.add_argument(
        "--sigma-density",
        metavar="sigmaValue",
        type=float,
        default=0.1,
        help="""Sigma value for calculation of nuclear density distributions.""",
    )
    parser.add_argument(
        "--description",
        type=str,
        nargs="*",
        help="""Space separated condition:description couples.
		'condition' are the name of condition folders.
		'description' are descriptive labels used in plots instead of folder names.
		Use -- after the last one.""",
    )
    parser.add_argument(
        "-t",
        "--threads",
        metavar="nthreads",
        type=int,
        default=1,
        help="""Number of threads to be used for parallelization. Increasing the
		number of threads might increase the required amount of RAM.""",
    )
    parser.add_argument(
        "--note", type=str, help="""Dataset/Analysis description. Use double quotes."""
    )
    regexp = "^(?P<channel_name>[^/]*)\\.(?P<channel_id>channel[0-9]+)"
    regexp += "\\.(?P<series_id>series[0-9]+)(?P<ext>(_cmle)?\\.tif)$"
    parser.add_argument(
        "--regexp",
        type=str,
        help="""Advanced. Regular expression to identify tif images.""",
        default=regexp,
    )
    parser.add_argument(
        "--nbins",
        type=int,
        default=200,
        help="Number of bins for profile calculation. Default: 200",
    )

    # Flag parameters
    parser.add_argument(
        "--shift-correct",
        action="store_const",
        help="""Perform shift correction.""",
        const=True,
        default=False,
    )
    parser.add_argument(
        "--no-hole-filling",
        action="store_const",
        help="""Do not fill holes in segmented masks.""",
        const=True,
        default=False,
    )
    parser.add_argument(
        "--labeled",
        action="store_const",
        dest="labeled",
        const=True,
        default=False,
        help="Export labeled masks instead of binary.",
    )
    parser.add_argument(
        "--compressed",
        action="store_const",
        dest="compressed",
        const=True,
        default=False,
        help="Generate compressed TIF binary masks.",
    )
    parser.add_argument(
        "-r",
        "--rescale-deconvolved",
        action="store_const",
        help="""Perform rescaling of deconvolved images. Requires Huygens
		Professional v4.5 log file for an image to be rescaled.""",
        const=True,
        default=False,
    )
    parser.add_argument(
        "-n",
        "--normalize-distance",
        action="store_const",
        help="""Perform distance normalization. Necessary to compare nuclei
		with different radius.""",
        const=True,
        default=False,
    )
    parser.add_argument(
        "-u",
        "--DEBUG-MODE",
        action="store_const",
        help="""Debugging mode.""",
        const=True,
        default=False,
    )
    parser.add_argument(
        "-y",
        "--do-all",
        action="store_const",
        help="""Do not ask for settings confirmation and proceed.""",
        const=True,
        default=False,
    )

    # Version flag
    parser.add_argument(
        "--version",
        action="version",
        version="%s v%s"
        % (
            sys.argv[0],
            version,
        ),
    )

    # Parse arguments
    args = parser.parse_args()

    # Additional checks
    args.threads = check_threads(args.threads)

    # RUN ==========================================================================

    # Create pyGPSeq analyzer instance
    gpi = gp.anim.Main(ncores=args.threads)

    # Steps to be skipped
    dskip = {"inst": 1, "seg": 2, "an": 3, "box": 3.5, "plot": 4, "report": 5}
    if None is not args.skip:
        gpi.skip = [dskip[e] for e in args.skip]

    # Channel names
    gpi.sig_names = tuple(args.sig_channels)
    gpi.dna_names = tuple(args.dna_channels)

    # Data directory
    gpi.basedir = args.inDir

    # Output directory
    gpi.outdir = args.outDir
    if not os.path.isdir(gpi.outdir):
        os.mkdir(gpi.outdir)

    assert args.seg_type in gp.const.SEG_ARG_LABELS
    gpi.seg_type = gp.const.SEG_ARG_LABELS.index(args.seg_type)
    gpi.mask_folder = args.mask_folder
    gpi.mask_prefix = args.mask_prefix
    gpi.labeled = args.labeled
    gpi.compressed = args.compressed

    assert args.an_type in gp.const.AN_ARG_LABELS
    gpi.an_type = gp.const.AN_ARG_LABELS.index(args.an_type)
    assert args.mid_type in gp.const.MID_SEC_ARG_LABELS
    gpi.mid_type = gp.const.MID_SEC_ARG_LABELS.index(args.mid_type)
    assert args.dist_type in gp.const.LD_ARG_LABELS
    gpi.dist_type = gp.const.LD_ARG_LABELS.index(args.dist_type)

    # Voxel aspect proportions (or sizes, ZYX)
    gpi.aspect = tuple(args.aspect)
    gpi.umes = args.umes

    # Minimum percentage of stack to be occupied by a cell
    gpi.min_z_size = args.min_z
    gpi.radius_interval = (args.min_radius, float("inf"))

    # Do hole filling
    gpi.do_fill_holes = not args.no_hole_filling

    # Nuclear selection
    dnsel = {"size": 0, "surf": 1, "shape": 2, "sumI": 3, "meanI": 4, "flat_size": 5}
    arsel = ["size", "surf", "shape", "sumI", "meanI", "flat_size"]
    gpi.nsf = tuple(dnsel[e] for e in args.nuclear_sel)
    readable_nsf = "*NONE*"
    if len(gpi.nsf) != 0:
        readable_nsf = " ".join([str(arsel[i]) for i in gpi.nsf])

    # Regular expression to identify image files
    gpi.reg = args.regexp

    # Where to save the run log
    if None is args.logpath:
        gpi.logpath = args.outDir + "/" + gpi.gen_log_name()
    else:
        gpi.logpath = args.logpath

    # Perform deconvolved image rescaling?
    gpi.rescale_deconvolved = args.rescale_deconvolved
    gpi.correct_shift = args.shift_correct

    # Normalize distance?
    gpi.normalize_distance = args.normalize_distance

    # Sigma
    gpi.sigma_smooth = args.sigma_smooth
    gpi.sigma_density = args.sigma_density
    gpi.nbins = args.nbins

    # Better condition naming
    if None is not args.description:
        for descr in args.description:
            c, d = descr.split(":")
            gpi.cdescr[c] = d
    cdescr_k = list(gpi.cdescr.keys())
    cdescr_k.sort()
    readable_cdescr = [str(k) + " => " + str(gpi.cdescr[k]) for k in cdescr_k]
    if not readable_cdescr:
        readable_cdescr = ["*NONE*"]

    # Notes
    if None is not args.note:
        gpi.notes = args.note

    # Debugging mode
    gpi.debugging = args.DEBUG_MODE

    # Show current settings
    ssettings = print_settings(gpi, args, readable_cdescr, readable_nsf)
    if not args.do_all:
        ask("Confirm settings and proceed?")

    # Save settings
    with open("%s.settings.txt" % os.path.splitext(gpi.logpath)[0], "w+") as OH:
        OH.write("%s\n" % " ".join(sys.argv))
        OH.write("@%s\n\n" % datetime.datetime.now())
        OH.write("%s\n" % sys.version)
        export_settings(OH, ssettings)
    print(args.regexp)
    # Start the analysis
    gpi = gpi.run()

    # End --------------------------------------------------------------------------

    ################################################################################

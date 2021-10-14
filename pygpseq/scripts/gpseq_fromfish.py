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
# Date: 20170718
# Project: GPSeq/FISH
# Description: Calculate radial position of dots in cells
#
# ------------------------------------------------------------------------------


# DEPENDENCIES =================================================================

import matplotlib.pyplot as plt

import argparse
import datetime
from joblib import Parallel, delayed
import numpy as np
import os
import pandas as pd
import pickle
import re
import sys

from ggc.prompt import ask
from ggc.args import check_threads, export_settings

import pygpseq as gp
from pygpseq.tools import image as imt
from pygpseq.fish.dot import add_allele, add_allele_polarity
from pygpseq.fish.nucleus import flag_G1_cells, plot_nuclei_aggregated
from pygpseq.fish.image import analyze_field_of_view

# FUNCTIONS ====================================================================
version = "7.0.2"


def print_settings(args, clear=True):
    """Show input settings, for confirmation.

    Args:
        args (Namespace): arguments parsed by argparse.
        clear (bool): clear screen before printing.
    """
    s = " # GPSeq analysis of FISH signals v%s\n" % version

    s += """
    --- INPUT =======

        FISH table : %s
    Image folder : %s
    Output folder : %s

        Mask folder : %s
        Mask prefix : '%s'

    --- ANALYSIS ====

        Dilation : %d
    Aspect (Z Y X) : %s
    Distance type : %s
Skipped channels : %s
    Pole fraction : %.3f
            #bins : %d

    --- FLAGS =======

        2D masks : '%s'
            Labeled : %r
        Compressed : %r
    Dilate over Z : %r
Use dilation only
for dot assignment : %r

        Plot all : %r
Plot compartments : %r

    --- ADVANCED ====

    Input regexp : %s
            Delim : '%s'
            Threads : %d
        Debug mode : %r
    """ % (
        args.dotCoords,
        args.imdir,
        args.outdir,
        args.mask_folder,
        args.mask_prefix,
        args.dilate,
        str(args.aspect),
        args.dist_type,
        str(args.skip_channels),
        args.pole,
        args.nbins,
        args.manual_2d_masks,
        args.labeled,
        args.compressed,
        args.doZdilation,
        args.dilate_for_assignment_only,
        not args.noplot,
        not args.no_compartment_plot,
        args.inreg,
        args.delim,
        args.threads,
        args.DEBUG_MODE,
    )

    if clear:
        print("\033[H\033[J%s" % s)
    else:
        print(s)
    return s


def run():

    # PARAMETERS ===================================================================

    if any([x in sys.argv for x in ["--help2", "-H"]]):
        print(
            """
    "gpseq_fromfish" can be used to calculate radial position of FISH signals (dots)
    in nuclei. Dots coordinates are expected to be 1-indexed. If 0-indexed, use the
    '-0' flag. The image dimensions are expected to be: Z (slices), Y (columns) and
    X (rows).

    Images are expected to follow DOTTER filename notation: "channel_series.tif".
    Identified images are first re-scaled (if deconvolved with Huygens software),
    then a global (Otsu) and local (median) thresholds are combined to binarize the
    image in 3D. Holes are filled in 3D and a closing operation is performed
    to remove small objects. Objects are filtered based on volume and Z size, and
    those touching the XY contour of the image are discarded. If the --labeled
    option is used, the generated images have identified objects labeled with
    different intensity values. Both compressed and normal TIFF files are
    compatible as input/output. To produce compressed TIFF files as output use the
    --compressed option.

    If you previously segmented your images (i.e., produced masks), provide the path
    to the folder containing the masks using -m, and the prefix for the mask name
    with -M. For example, with '-m /ex/mask_dir/ -M mask_', in the case of the image
    file '1.tif', the script will look for the mask at "/ex/mask_dir/mask_1.tif".
    If the mask can be found, it will be used, otherwise it will be generated and
    then saved. This can be used also to export masks as tifs to the folder
    specified with -m.

    If needed, 2D masks can be 'enforced' to the 3D masks. This is performed by
    applying a logical AND operation between each slice of the 3D mask and the 2D
    mask. To achieve this, use the -2 option, specifying the path to a folder
    containing the 2D mask tiff files (filename should match input image). The
    enforcing is performed during the segmentation process. Thus, if segmentation is
    skipped (e.g., if a mask is provided as input with -m), this is skipped too.

    The G1 selection is actually a selection of the most represented cell sub-
    -population based on flatten area and integral of DNA stain intensity. In other
    words, it will selected the most represented cell cycle phase in your cell
    population (generally, G1). Basically, the major peak is identified and its FWHM
    is used as a range for the filter on both nuclear features.

    Then, each dot is assigned to the poles (2), bottom center (1) or top center (0)
    compartments. The pole compartments are defined as the 25%% of the whole volume,
    by cutting perpendicularly to the major axis. Such definition can be changed
    with the -P (--pole) option specifying a different axis fraction. Information on
    goodnes of fit is reported and a plot is provided with each nucleus rotated and
    centered.

    Also, each dot is assigned an `Allele` label: -1 to dots from cells with more
    than 2 dots in that channel, 0 to dots from cells with less than 2 dots in that
    channel, and in the case of cells with exactly 2 dots in that channel the most
    peripheral dot is labeled with 1, and the most central with 2. This is useful
    when expecting 2 dots per cell-channel, i.e., in the case of one probe per
    channel in a diploid cell line.

    Plotting can be turned off generally with the --noplot flag, and specifically
    with the --no-compartment-plot flag.

    More details can be found in the online documentation, at:
    https://github.com/ggirelli/pygpseq/wiki/
    """
        )
        sys.exit()

    # Add script description
    parser = argparse.ArgumentParser(
        description="", formatter_class=argparse.RawDescriptionHelpFormatter
    )

    # Add mandatory arguments
    parser.add_argument(
        "dotCoords", type=str, help="Dot coordinates table generated by DOTTER."
    )
    parser.add_argument(
        "imdir", type=str, help="Path to folder containing deconvolved tiff images."
    )
    parser.add_argument(
        "outdir", type=str, help="Path to output folder (created if does not exist)."
    )

    # Optional parameters
    parser.add_argument(
        "-H",
        "--help2",
        action="store_const",
        help="""Shows additional information in a readable format and quit.""",
        const=True,
        default=False,
    )
    parser.add_argument(
        "-a",
        "--aspect",
        type=float,
        nargs=3,
        help="""Physical size of Z, Y and X voxel sides.
        Default: 300.0 130.0 130.0""",
        metavar=("Z", "Y", "X"),
        default=[300.0, 130.0, 130.0],
    )
    parser.add_argument(
        "-d",
        "--delim",
        metavar="sep",
        type=str,
        help="""Input table delimiter. Default: ','""",
        default=",",
    )
    parser.add_argument(
        "-D",
        "--dilate",
        metavar="npx",
        type=int,
        help="""Number of pixels for nuclear mask dilation. It is automatically
        scaled based on the specified aspect to be isotropic in 3D. Default: 0""",
        default=0,
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
        "-t",
        "--threads",
        metavar="nthreads",
        type=int,
        help="""Number of threads for parallelization. Default: 1""",
        default=1,
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
        help="""Path to folder with 2D masks as tiff files with matching name,
        to combine with 3D masks. Used only if no mask was provided through -m.""",
    )
    parser.add_argument(
        "-P",
        "--pole",
        metavar="axis_fraction",
        type=float,
        help="""Fraction of the major nuclear axis to identify a pole.
        Should be in the [0, .5] interval. Default: .25.""",
        default=0.25,
    )
    parser.add_argument(
        "-C",
        "--skip-channels",
        metavar="channel",
        type=str,
        help="""List of space-separated channels to skip. Example: 'a700 ir800'.
        The channel names are forced to lower-case.""",
        nargs="+",
    )
    default_inreg = "^.*\.tiff?$"
    parser.add_argument(
        "--inreg",
        type=str,
        help="""regular expression to identify images from imdir.
        Default: '%s'"""
        % (default_inreg,),
        default=default_inreg,
    )
    parser.add_argument(
        "--nbins",
        type=int,
        default=200,
        help="Number of bins for density profile calculation. Default: 200",
    )

    # Add flags
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
    parser.add_argument(
        "-0",
        "--zero-indexed",
        action="store_const",
        help="""Coordinates are 0-indexed instead of 1-indexed (default).""",
        const=True,
        default=False,
    )
    parser.add_argument(
        "--labeled",
        action="store_const",
        dest="labeled",
        const=True,
        default=False,
        help="Import/Export masks as labeled instead of binary.",
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
        "--dots-assignment-only",
        action="store_const",
        dest="dilate_for_assignment_only",
        const=True,
        default=False,
        help="""Discard dilated masks after dots assignment. In this way, dots
        located outside the non-dilated mask have a lamina distance of 0.""",
    )
    parser.add_argument(
        "--dilate-Z",
        action="store_const",
        dest="doZdilation",
        const=True,
        default=False,
        help="Turn on dilation over Z.",
    )
    parser.add_argument(
        "--noplot",
        action="store_const",
        dest="noplot",
        const=True,
        default=False,
        help="Do not produce any plots.",
    )
    parser.add_argument(
        "--no-compart-plot",
        action="store_const",
        dest="no_compartment_plot",
        const=True,
        default=False,
        help="Do not produce compartments-related plots.",
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

    # Manipulate input -------------------------------------------------------------

    # Assign to in-script variables
    delim = args.delim
    (az, ay, ax) = args.aspect

    # Assign default value
    if not type(None) == type(args.skip_channels):
        args.skip_channels = [c.lower() for c in args.skip_channels]

    # Rename output directory if it exists already
    # And add trailing slashes
    if not args.outdir[-1] == "/":
        while not os.path.isdir(args.outdir) and os.path.exists(args.outdir):
            args.outdir += "_"
        args.outdir += "/"
    if not args.imdir[-1] in ["/\\"]:
        args.imdir += "/"

    # Adjust number of threads
    args.threads = check_threads(args.threads)

    # Limit pole fraction
    if 0 >= args.pole:
        args.pole = 0
    if 0.5 < args.pole:
        args.pole = 0.5

    if not type(None) == type(args.manual_2d_masks):
        assert os.path.isdir(args.manual_2d_masks), "folder not found."

    # Params -----------------------------------------------------------------------

    seg_type = gp.const.SEG_3D
    an_type = gp.const.AN_3D

    # Additional checks ------------------------------------------------------------

    if 0 != args.dilate:
        assert_msg = "cannot apply dilation on images with different X/Y aspect."
        assert ax == ay, assert_msg

    # RUN ==========================================================================

    ssettings = print_settings(args)
    if not args.do_all:
        ask("Confirm settings and proceed?")

    # Create output folder
    if not os.path.isdir(args.outdir):
        os.mkdir(args.outdir)

    # Save settings
    with open("%s/settings.txt" % args.outdir, "w+") as OH:
        export_settings(OH, ssettings)

    ddir = ""
    if args.DEBUG_MODE:
        ddir = "%s/debugging" % args.outdir
        if not os.path.isdir(ddir) and not os.path.isfile(ddir):
            os.mkdir(ddir)

    # Build 3D isotropic structuring element for dilation
    if not args.doZdilation:
        istruct = imt.mkIsoStruct(args.dilate, (0, ay, ax))
    else:
        istruct = imt.mkIsoStruct(args.dilate, (az, ay, ax))

    # Check if isotropic dilation is required
    if az != ax or not args.doZdilation:
        t = np.array(istruct.shape) / 2
        msg = "  Isotropic dilation of "
        msg += "(%d, %d, %d) px in ZYX, respectively." % tuple(t.tolist())
        print(msg)

    # Input ------------------------------------------------------------------------

    # Read table
    t = pd.read_csv(args.dotCoords, delim)

    # Check that all required columns are present
    req_col = ["File", "Channel", "x", "y", "z"]
    misscols = [cn for cn in req_col if cn not in t.columns]
    assert 0 == len(misscols), "missing required column(s): %s" % str(misscols)

    # Remove channels to be skipped from the table
    old_size = t.shape[0]
    if type(None) != type(args.skip_channels):
        print("  Skipping channels: %s" % str(args.skip_channels))
        for i in t.index:
            if t.loc[i, "Channel"] in args.skip_channels:
                t = t.ix[t.index != i, :]
        t.index = range(t.shape[0])
        print("  >>> Reduced table from %d to %d rows." % (old_size, t.shape[0]))

    # Add new empty columns
    t["dilation"] = args.dilate
    t["version"] = version

    # Identify images --------------------------------------------------------------

    # Extract FoV number
    t["File"] = [int(os.path.splitext(os.path.basename(f))[0]) for f in t["File"]]

    # Round up float coordinates and convert to integer to get matricial coordinates
    t["xi"] = np.round(t["x"]).astype("i")
    t["yi"] = np.round(t["y"]).astype("i")
    t["zi"] = np.round(t["z"]).astype("i")

    # Switch to 0-indexing (e.g., Python) if 1-indexing (e.g., MATLAB)
    if not args.zero_indexed:
        t["x"] = t["x"] - 1.0
        t["y"] = t["y"] - 1.0
        t["z"] = t["z"] - 1.0
        t["xi"] = t["xi"] - 1
        t["yi"] = t["yi"] - 1
        t["zi"] = t["zi"] - 1

    # Identify tiff images
    imlist = []
    for (dirpath, dirnames, filenames) in os.walk(args.imdir):
        imlist.extend([f for f in filenames if 0 != len(re.findall(args.inreg, f))])
        break

    # Assign field of views to images
    im2fov = {}
    for i in set(t["File"]):
        imsel = [im for im in imlist if "%03d" % (i,) in im]
        if not 0 == len(imsel):
            im2fov[i] = os.path.join(args.imdir, imsel[0])
        else:
            t = t.ix[t["File"] != i, :]
            print("  Missing image for field #%d, skipped." % (i,))
    t.index = range(t.shape[0])

    # Iterate ----------------------------------------------------------------------

    print("  > Analyzing fields of view... [n.threads=%d]" % (args.threads,))
    kwargs = {
        "data": t,
        "im2fov": im2fov,
        "dilate_factor": args.dilate,
        "istruct": istruct,
        "discard_dilation_mode": args.dilate_for_assignment_only,
        "aspect": args.aspect,
        "mask_dir": args.mask_folder,
        "mask_prefix": args.mask_prefix,
        "plotCompartments": not args.no_compartment_plot,
        "pole_fraction": args.pole,
        "outdir": args.outdir,
        "noplot": args.noplot,
        "labeled": args.labeled,
        "compressed": args.compressed,
        "mask2d_dir": args.manual_2d_masks,
        "an_type": an_type,
        "seg_type": seg_type,
        "dist_type": gp.const.LD_ARG_LABELS.index(args.dist_type),
        "nbins": args.nbins,
        "debug": args.DEBUG_MODE,
        "debug_dir": ddir,
    }
    if 1 != args.threads:
        anData = Parallel(n_jobs=args.threads, verbose=11)(
            delayed(analyze_field_of_view)(ii, **kwargs) for ii in im2fov.keys()
        )
    else:
        anData = []
        for k in im2fov.keys():
            anData.append(analyze_field_of_view(k, verbose=True, **kwargs))

    # Parse output and store log report --------------------------------------------

    with open(os.path.join(args.outdir, "fov_analysis.log"), "w+") as hlog:
        nuclei = []
        tvdata = []
        t = []
        t_columns = []
        density_profile = []
        volume_profile = []
        for data in anData:
            if type(None) != type(data):
                curnuclei, subt, tvcomp, dp, nv, msg = data
                hlog.write(msg)

                nuclei.extend(curnuclei.values())
                tvdata.append(tvcomp)
                t.append(subt)
                density_profile.append(dp)
                volume_profile.append(nv)

                for c in subt.columns:
                    if not c in t_columns:
                        t_columns.append(c)

        assert 0 != len(t)
        t = pd.concat(t).loc[:, t_columns]
        tvdata = pd.concat(tvdata)
        tvdata.index = range(tvdata.shape[0])

    # Plot aggregated FISH view ----------------------------------------------------
    if not args.no_compartment_plot:
        aggdir = os.path.join(args.outdir, "agg_vis/")
        if not os.path.isdir(aggdir):
            os.mkdir(aggdir)
        plot_nuclei_aggregated(t, tvdata, args.aspect, aggdir)

    # Prepare density profile ------------------------------------------------------
    dp = pd.concat(density_profile)
    dp["c"] = args.dotCoords
    outname = "%s/density_profile.tsv" % (args.outdir,)
    dp.to_csv(outname, sep="\t", index=False)

    # Prepare volume profile -------------------------------------------------------
    nv = pd.concat(volume_profile)
    nv["c"] = args.dotCoords
    outname = "%s/volume_profile.tsv" % (args.outdir,)
    nv.to_csv(outname, sep="\t", index=False)

    # Identify G1 cells ------------------------------------------------------------
    bname = os.path.basename(args.dotCoords).split(".")
    bname[-1] = "tsv"
    bname = ".".join(bname)
    t = flag_G1_cells(t, nuclei, args.outdir, args.dilate, bname)

    # Export -----------------------------------------------------------------------

    # Export compartment volume data
    outname = "%s/nuclear_compartment.out.dilate%d.%s" % (
        args.outdir,
        args.dilate,
        bname,
    )
    tvdata.to_csv(outname, sep="\t", index=False)

    # Export nuclei object vector
    with open("%s/nuclei.pickle" % args.outdir, "wb+") as f:
        pickle.dump(nuclei, f)

    # Export table before allele labeling
    # outname = "%s/wCentr.out.noAllele.dilate%d.%s" % (
    #   args.outdir, args.dilate, bname)
    # t.to_csv(outname, sep = '\t', index = False)

    # Add allele information -------------------------------------------------------
    print("  - Adding allele information...")
    t = add_allele(t)

    # Calculate angle on nucleus centroid between alleles --------------------------
    print("  - Adding allele polarity information...")
    t = add_allele_polarity(t, nuclei, args.aspect)

    # Write output -----------------------------------------------------------------
    outname = "%s/wCentr.out.dilate%d.%s" % (args.outdir, args.dilate, bname)
    t.to_csv(outname, sep="\t", index=False)

    # END ==========================================================================

    ################################################################################

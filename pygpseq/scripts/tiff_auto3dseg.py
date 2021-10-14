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
# Date: 20171205
# Project: bioimaging
# Description: automatic 3D segmentation of nuclear staining.
#
# ------------------------------------------------------------------------------


# DEPENDENCIES =================================================================

import argparse
from joblib import Parallel, delayed
import math
import multiprocessing
import numpy as np
import os
import re
from scipy import ndimage as ndi
import sys
import tifffile
from tqdm import tqdm
import warnings

from ggc.prompt import ask
from ggc.args import check_threads, export_settings

from skimage.filters import threshold_local, threshold_otsu
import skimage.io as io
from skimage.measure import label
from skimage.morphology import closing, cube, square
from skimage.segmentation import clear_border

from pygpseq import const
from pygpseq.tools import Binarize
from pygpseq.tools import image as imt
from pygpseq.tools import path as pt
from pygpseq.tools import plot
from pygpseq.tools import stat as stt
from pygpseq.tools import vector as vt


# FUNCTIONS ====================================================================
version = "3.1.1"


def run_segmentation(args, imgpath, imgdir, radius_interval):
    # Perform 3D segmentation of nuclear staining image.
    #
    # Args:
    #   imgpath (string): input image file name.
    #   imgdir (string): input image folder.
    #
    # Returns:
    #   string: path to output image.

    # Preparation --------------------------------------------------------------

    # Read image
    irf = imt.get_rescaling_factor(os.path.join(imgdir, imgpath))
    img = imt.read_tiff(os.path.join(imgdir, imgpath), 3, rescale=irf)

    # binarize -----------------------------------------------------------------

    binarization = Binarize(
        an_type=const.AN_3D,
        seg_type=const.SEG_3D,
        verbose=False,
        radius_interval=radius_interval,
        min_z_size=args.min_Z,
        do_clear_Z_borders=args.do_clear_Z,
        adp_neigh=args.neighbour,
    )

    mask2d = None
    if args.combineWith2D:
        mask2d_path = os.path.join(args.manual_2d_masks, imgpath)
        if os.path.isfile(mask2d_path):
            mask2d = imt.read_tiff(mask2d_path)
        else:
            print("Warning: 2D mask not found at '%s'" % mask2d_path)

    if type(None) != type(mask2d):
        (mask, thr, log) = binarization.run(img, mask2d, args.labeled)
    else:
        (mask, thr, log) = binarization.run(img)

    # Filter based on object size
    mask, tmp = binarization.filter_obj_XY_size(mask)
    mask, tmp = binarization.filter_obj_Z_size(mask)

    # Perform dilate-fill-erode operation
    if args.dilate_fill_erode != 0:
        strel = args.dilate_fill_erode
        strel = cube(strel) if len(mask.shape) == 3 else square(strel)
        mask = imt.dilate_fill_erode(mask, strel)

    # Label nuclei if not done already
    if not (args.combineWith2D and args.labeled):
        L = label(mask)
    else:
        L = mask
        if type(None) != type(mask2d):
            # Re-assign extra-mask labels
            L = binarization.combine_2d_mask(L, mask2d, args.labeled)

    # Output -------------------------------------------------------------------
    outpath = "%s%s" % (args.outFolder, args.outprefix + imgpath)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        if not args.labeled:
            L[np.nonzero(L)] = 255
        plot.save_tif(outpath, L, "uint8", args.compressed, "ZYX")
    print("Segmentation output written to %s" % (outpath,))

    return outpath


def print_settings(args, radius_interval, clear=True):
    """Show input settings, for confirmation.

    Args:
        args (Namespace): arguments parsed by argparse.
        clear (bool): clear screen before printing.
    """
    s = " # Automatic 3D segmentation v%s\n" % version

    s += """
---------- SETTING :  VALUE ----------

Input directory :  '%s'
Output directory :  '%s'

    Mask prefix :  '%s'
    Neighbourhood :  %d
        2D masks : '%s'
        Labeled :  %r
        Compressed :  %r

Dilate-fill-erode :  %d
Minimum Z portion :  %.2f
    Minimum radius :  [%.2f, %.2f] vx
        Clear Z :  %r

        Threads :  %d
            Regexp :  '%s'

    """ % (
        args.imgFolder,
        args.outFolder,
        args.outprefix,
        args.neighbour,
        args.manual_2d_masks,
        args.labeled,
        args.compressed,
        args.dilate_fill_erode,
        args.min_Z,
        radius_interval[0],
        radius_interval[1],
        args.do_clear_Z,
        args.threads,
        args.inreg,
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
    Perform automatic 3D segmentation of DNA staining. Images are first identified
    based on a regular expression matching the file name. Then, they are first
    re-scaled if deconvolved with Huygens software, then a global (Otsu) and
    local (median) thresholds are combined to binarize the image in 3D. Then, holes
    are filled in 3D and a closing operation is performed to remove small objects.
    Objects are filtered based on volume and Z size, and those touching the XY
    contour of the image are discarded. The generated images have identified objects
    labeled with different intensity levels.
    """,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # Add mandatory arguments
    parser.add_argument(
        "imgFolder", type=str, help="Path to folder containing deconvolved tiff images."
    )
    parser.add_argument(
        "outFolder",
        type=str,
        help="""Path to output folder where imt.binarized images will be stored
        (created if does not exist).""",
    )

    # Optional parameters
    default_inreg = "^.*\.tiff?$"
    parser.add_argument(
        "--inreg",
        type=str,
        help="""regular expression to identify images from imgFolder.
        Default: '%s'"""
        % (default_inreg,),
        default=default_inreg,
    )
    parser.add_argument(
        "--outprefix",
        type=str,
        help="""prefix to add to the name of output imt.binarized images.
        Default: 'mask_', 'cmask_' if --compressed is used.""",
        default="mask_",
    )
    parser.add_argument(
        "--neighbour",
        type=int,
        help="""Side of neighbourhood square/cube. Default: 101""",
        default=101,
    )
    parser.add_argument(
        "--radius",
        type=float,
        nargs=2,
        help="""Range of object radii [vx] to be considered a nucleus.
        Default: [10, Inf]""",
        default=[10.0, float("Inf")],
    )
    parser.add_argument(
        "--min-Z",
        type=float,
        help="""Minimum fraction of stack occupied by an object to be considered a
        nucleus. Default: .25""",
        default=0.25,
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        help="""Number of threads for parallelization. Default: 1""",
        default=1,
    )
    parser.add_argument(
        "-2",
        "--manual-2d-masks",
        type=str,
        metavar="MAN2DDIR",
        help="""Path to folder with 2D masks with matching name,
        to combine with 3D masks.""",
    )
    parser.add_argument(
        "-F",
        "--dilate-fill-erode",
        type=int,
        metavar="DFE",
        help="""Number of pixels for dilation/erosion in a dilate-fill-erode
        operation. Default: 10. Set to 0 to skip.""",
        default=10,
    )

    # Flags
    parser.add_argument(
        "--clear-Z",
        action="store_const",
        dest="do_clear_Z",
        const=True,
        default=False,
        help="""Remove objects touching the bottom/top of the stack.""",
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
        version="%s %s"
        % (
            sys.argv[0],
            version,
        ),
    )

    # Parse arguments
    args = parser.parse_args()

    # Assign to in-script variables
    inreg = re.compile(args.inreg)
    radius_interval = [args.radius[0], args.radius[1]]

    if args.compressed and args.outprefix == "mask_":
        args.outprefix = "cmask_"

    args.combineWith2D = type(None) != type(args.manual_2d_masks)
    if args.combineWith2D:
        assert_msg = "2D mask folder not found, '%s'" % args.manual_2d_masks
        assert os.path.isdir(args.manual_2d_masks), assert_msg

    # Additional checks
    args.threads = check_threads(args.threads)

    # RUN ==========================================================================

    # Show current settings
    ssettings = print_settings(args, radius_interval)
    if not args.do_all:
        ask("Confirm settings and proceed?")

    # Add trailing slashes
    args.imgFolder = pt.add_trailing_slash(args.imgFolder)
    args.outFolder = pt.add_trailing_slash(args.outFolder)

    # Stop if input folder does not exist
    if not os.path.isdir(args.imgFolder):
        sys.exit("!ERROR! Image folder not found: %s" % (args.imgFolder,))

    # Create output folder
    if not os.path.isdir(args.outFolder):
        os.mkdir(args.outFolder)

    # Save confirmed settings
    with open(os.path.join(args.outFolder, "settings.txt"), "w+") as OH:
        export_settings(OH, ssettings)

    # Identify images --------------------------------------------------------------

    # Identify tiff images
    imglist = [
        f
        for f in os.listdir(args.imgFolder)
        if os.path.isfile(os.path.join(args.imgFolder, f))
        and type(None) != type(re.match(inreg, f))
    ]

    print("Found %d image(s) to segment. Starting..." % len(imglist))

    # Start iteration --------------------------------------------------------------

    if args.threads == 1:
        for imgpath in tqdm(imglist):
            run_segmentation(args, imgpath, args.imgFolder, radius_interval)
    else:
        outlist = Parallel(n_jobs=args.threads, verbose=11)(
            delayed(run_segmentation)(args, imgpath, args.imgFolder, radius_interval)
            for imgpath in imglist
        )

    # END ==========================================================================

    ################################################################################

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
# Date: 20180302
# Project: bioimaging
# Description: Uncompress compressed tiff.
#
# ------------------------------------------------------------------------------


# DEPENDENCIES =================================================================

import argparse
from joblib import Parallel, delayed
import multiprocessing
import numpy as np
import os
import re
from scipy import ndimage as ndi
import sys
import tifffile
import warnings

from pygpseq.tools import image as imt
from pygpseq.tools import path as pt
from pygpseq.tools import plot
from pygpseq.tools.io import printout

from skimage.filters import threshold_local, threshold_otsu
import skimage.io as io
from skimage.measure import label
from skimage.morphology import closing, cube, square
from skimage.segmentation import clear_border


# FUNCTIONS ====================================================================


def run(imgpath, imgdir, outdir, outpath=None, compress=None):
    # Perform 3D segmentation of nuclear staining image.
    #
    # Args:
    #   imgpath (string): input image file name.
    #   imgdir (string): input image folder.
    #
    # Returns:
    #   string: path to output image.

    if type(None) == type(compress):
        compress = False

    # Preparation --------------------------------------------------------------

    # Read image
    img = imt.read_tiff(os.path.join(imgdir, imgpath))

    # Write image
    if type(None) == type(outpath):
        outpath = imgpath

    if not compress:
        plot.save_tif(
            os.path.join(outdir, outpath), img, imt.get_dtype(img.max()), False
        )
        label = "Uncompressed"
    else:
        plot.save_tif(
            os.path.join(outdir, outpath), img, imt.get_dtype(img.max()), True
        )
        label = "Compressed"

    print("%s '%s'." % (label, os.path.join(imgdir, imgpath)))
    return os.path.join(outdir, outpath)


def run():

    # PARAMETERS ===================================================================

    # Add script description
    parser = argparse.ArgumentParser(
        description="""
    (Un)compress TIFF images. Provide either a single input/output image path, or 
    input/output folder paths. In case of folder input/output, all tiff files in the
    inout folder with a name matching the specified pattern are (un)compressed and
    saved to the output folder.

    When (un)compressing multiple files, the --threads option allows to parallelize
    on multiple threads. Disk read/write operations become the bottleneck when
    parallelizing, thus working on a SSD is advised.
    """,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # Add mandatory arguments
    parser.add_argument(
        "input",
        type=str,
        nargs=1,
        help="""Path to the TIFF image to uncompress, or to a folder containing
        multiple TIFF images. In the latter case, the --inreg pattern is used to
        identify the images.""",
    )
    parser.add_argument(
        "output",
        type=str,
        nargs=1,
        help="""Path to output TIFF image, or output folder if the input is a
        older.""",
    )

    # Optional parameters
    default_inreg = "^.*\.tiff?$"
    parser.add_argument(
        "--inreg",
        type=str,
        nargs=1,
        help="""Regular expression to identify image files, when the input is a
        folder. Default: '%s'"""
        % (default_inreg,),
        default=[default_inreg],
    )
    parser.add_argument(
        "--threads",
        type=int,
        nargs=1,
        help="""Number of threads for parallelization. Used only to uncompress
        multiple images (i.e., input is a folder). Default: 1""",
        default=[1],
    )

    # Add flags
    parser.add_argument(
        "-u",
        action="store_const",
        dest="doUncompress",
        const=True,
        default=False,
        help="Uncompress TIFF files.",
    )
    parser.add_argument(
        "-c",
        action="store_const",
        dest="doCompress",
        const=True,
        default=False,
        help="Compress TIFF files.",
    )

    # Version flag
    version = "1.0.1"
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
    doCompress = args.doCompress
    doUncompress = args.doUncompress
    inpattern = re.compile(args.inreg[0])
    ncores = args.threads[0]

    # Additional checks
    if not doCompress and not doUncompress:
        printout("""Please, use either -c (compress) or -u (uncompress).""", -2)
    if doCompress and doUncompress:
        printout("""Please, use either -c (compress) or -u (uncompress).""", -2)
    maxncores = multiprocessing.cpu_count()
    if maxncores < ncores:
        print("Lowered number of threads to maximum available: %d" % (maxncores))
        ncores = maxncores

    # Select multitple/single operation style
    do_multiple = False
    if os.path.isdir(args.input[0]):
        do_multiple = True
        imgdir = args.input[0]
        if os.path.isfile(args.output[0]):
            printout(
                """Output must be a folder if the provided input is a folder
                too.""",
                -2,
            )
        else:
            outdir = args.output[0]
    elif os.path.isfile(args.input[0]):
        imgpath = args.input[0]
        if os.path.isdir(args.output[0]):
            printout(
                """Output can be a folder only if the provided input is a
                    folder too.""",
                -2,
            )
        else:
            outpath = args.output[0]
    else:
        printout("Input file not found: %s" % (args.input[0],), -2)

    # RUN ==========================================================================

    if do_multiple:
        # Uncompess multiple images ------------------------------------------------

        # Add trailing slashes
        imgdir = pt.add_trailing_slash(imgdir)
        outdir = pt.add_trailing_slash(outdir)

        # Stop if input folder does not exist
        if not os.path.isdir(imgdir):
            sys.exit("!ERROR! Image folder not found: %s" % (imgdir,))

        # Create output folder
        if not os.path.isdir(outdir):
            os.mkdir(outdir)

        # Identify images
        imglist = [
            f
            for f in os.listdir(imgdir)
            if os.path.isfile(os.path.join(imgdir, f))
            and type(None) != type(re.match(inpattern, f))
        ]

        # Start iteration
        outlist = Parallel(n_jobs=ncores)(
            delayed(run)(imgpath, imgdir, outdir, compress=doCompress)
            for imgpath in imglist
        )
    else:
        # Uncompress a single image ------------------------------------------------
        imgdir = os.path.dirname(imgpath)
        imgpath = os.path.basename(imgpath)
        outdir = os.path.dirname(outpath)
        outpath = os.path.basename(outpath)
        run(imgpath, imgdir, outdir, outpath, doCompress)

    # END ==========================================================================

    ################################################################################

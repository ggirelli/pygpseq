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
# Date: 20171129
# Project: bioimaging
# Description: extract intensity distribution on Z.
#
# ------------------------------------------------------------------------------


# DEPENDENCIES =================================================================

import matplotlib.pyplot as plt

import argparse
from joblib import Parallel, delayed
import math
import numpy as np
import os
import re
import sys
from tqdm import tqdm

from ggc.args import check_threads
from pygpseq.tools import image as imt, plot, stat


# FUNCTIONS ====================================================================


def mkPlot(pdata, path):
    """Generate pdf plot of sum intensity per Z slice.

    Args:
        pdata (dict): for each FoV, a dict with 'x' and 'y' paired coordinates.
        path (string): path to pdf output file.

    Returns:
        None: writes to disk.
    """

    plt.figure(figsize=[12, 8])

    xmax = max(max(f["x"]) for f in pdata.values())
    ymax = max(max(f["y"]) for f in pdata.values())

    for (f, data) in pdata.items():
        plt.plot(data["x"], data["y"], linewidth=0.5)

    plt.xlabel("Z-slice index")
    if args.intensity_sum:
        plt.ylabel("Intensity sum [a.u.]")
    else:
        plt.ylabel("Gradient magnitude [a.u.]")
    plt.title("Out-of-focus study")

    plt.legend(
        list(pdata.keys()),
        bbox_to_anchor=(1.04, 1),
        loc="upper left",
        prop={"size": 6},
    )
    plt.subplots_adjust(right=0.75)

    plt.gca().axvline(x=xmax * args.range[0] / 2, ymax=ymax, linestyle="--", color="k")
    plt.gca().axvline(
        x=xmax - xmax * args.range[0] / 2, ymax=ymax, linestyle="--", color="k"
    )

    plot.export(path)

    plt.show()


def isOOF(args, impath):
    # Read image
    im = imt.read_tiff("%s%s" % (args.imdir[0], impath))

    # Select first time frame
    while len(im.shape) > 3:
        im = im[0]

    # Iterate through slices
    intlist = []

    profile_data = {}
    sout = ""
    for zi in range(im.shape[0]):
        if args.intensity_sum:
            intlist.append(im[zi].sum())
        else:
            dx = stat.gpartial(im[zi, :, :], 1, 1)
            dy = stat.gpartial(im[zi, :, :], 2, 1)
            intlist.append(np.mean(np.mean((dx ** 2 + dy ** 2) ** (1 / 2))))

        # Output string
        sout += "%s\t%d\t%f\n" % (impath, zi + 1, intlist[zi])

        # If plot is required, update profile data
        if args.plot:
            if impath in profile_data:
                profile_data[impath]["x"].append(zi + 1)
                profile_data[impath]["y"].append(intlist[zi])
            else:
                profile_data[impath] = {"x": [zi + 1], "y": [intlist[zi]]}

    # Identify maximum slice
    maxid = intlist.index(max(intlist))
    hrange = im.shape[0] * args.range[0] / 2.0
    hstack = im.shape[0] / 2.0
    if maxid >= (hstack - hrange) and maxid <= (hstack + hrange):
        summary = "%s is in-focus.\n" % (impath,)
    else:
        summary = "%s is out-of-focus.\n" % (impath,)
        if args.rename:
            os.rename(
                "%s%s" % (args.imdir[0], impath),
                "%s%s.old" % (args.imdir[0], impath),
            )

    return (sout, summary, profile_data)


def run():

    # PARAMETERS ===================================================================

    # Add script description
    parser = argparse.ArgumentParser(
        description="""
    Calculate gradient magnitude over Z for every image in the input folder with a
    filename matching the --pattern. Use --range to change the in-focus definition.
    """,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # Add mandatory arguments
    parser.add_argument(
        "imdir", type=str, nargs=1, help="Path to folder with tiff images."
    )
    parser.add_argument("output", type=str, nargs=1, help="Path to output table file.")

    # Add arguments with default value
    parser.add_argument(
        "-r",
        "--range",
        type=float,
        nargs=1,
        metavar="range",
        help="""Specify %% of stack where the maximum of intensity distribution
        over Z is expected for an in-focus field of view. Default: 50%%""",
        default=[0.5],
    )
    parser.add_argument(
        "-p",
        "--pattern",
        type=str,
        nargs=1,
        metavar="regexp",
        help='''Provide a regular expression pattern matching the images in the
        image folder that you want to check. Default: "^.*\.tif$"''',
        default=["^.*\.tif$"],
    )
    parser.add_argument(
        "-t",
        "--threads",
        metavar="nthreads",
        type=int,
        help="""Number of threads for parallelization. Default: 1""",
        default=1,
    )

    # Flag arguments
    parser.add_argument(
        "-P",
        "--plot",
        action="store_const",
        help="""Generate pdf plot of intensity sum per Z-slice.""",
        const=True,
        default=False,
    )
    parser.add_argument(
        "-S",
        "--intensity-sum",
        action="store_const",
        help="""Use intensity sum instead of gradient magnitude.""",
        const=True,
        default=False,
    )
    parser.add_argument(
        "-R",
        "--rename",
        action="store_const",
        help="""Rename out-of-focus images by adding the '.old' suffix.""",
        const=True,
        default=False,
    )
    parser.add_argument(
        "-s",
        "--silent",
        action="store_const",
        help="""Silent run.""",
        const=True,
        default=False,
    )

    # Version flag
    version = "0.3.1"
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

    # Adjust number of threads
    args.threads = check_threads(args.threads)

    # RUN ==========================================================================

    # Add trailing slash to image folder path
    if args.imdir[0][-1] != "/":
        args.imdir[0] += "/"

    # Check that image folder path exists
    if not os.path.isdir(args.imdir[0]):
        sys.exit("!ERROR: specified imdir does not exist.\n%s" % (args.imdir[0],))

    # If plot is required, prepare plot_data
    if args.plot:
        profile_data = {}

    # Identify tiff images
    flist = []
    for (dirpath, dirnames, filenames) in os.walk(args.imdir[0]):
        flist.extend(filenames)
        break
    immatch = lambda f: type(None) != type(re.match(args.pattern[0], f))
    imlist = [f for f in flist if immatch(f)]

    # Iterate through fields of view
    if args.threads == 1:
        if args.silent:
            t = imlist
        else:
            t = tqdm(imlist)
            t.set_description(os.path.dirname(args.imdir[0]))

        data = [isOOF(args, impath) for impath in t]
    else:
        verbosity = 11 if not args.silent else 0
        data = Parallel(n_jobs=args.threads, verbose=verbosity)(
            delayed(isOOF)(args, impath) for impath in imlist
        )

    with open(args.output[0], "w+") as fout:
        with open("%s.log" % os.path.splitext(args.output[0])[0], "w+") as lout:
            for (s, summary, pdata) in data:
                fout.write(s)
                lout.write(summary)
                if args.plot:
                    profile_data.update(pdata)

    # Generate plot
    outpath = os.path.splitext(args.output[0])[0] + ".pdf"
    if args.plot:
        mkPlot(profile_data, outpath)

    # END ==========================================================================

    ################################################################################

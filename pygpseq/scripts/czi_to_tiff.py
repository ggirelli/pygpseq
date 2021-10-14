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
# Date: 20180918
# Project: bioimaging
# Description: convert czi images to tiff.
#
# ------------------------------------------------------------------------------


# DEPENDENCIES =================================================================

import argparse
import czifile
import numpy as np
import os
import sys
import tifffile
from tqdm import tqdm
import warnings
import xml.etree.ElementTree as ET

from pygpseq.tools import image as imt
from pygpseq.tools import plot
from pygpseq.tools.io import printout


# FUNCTIONS ====================================================================


def log_samples(images):
    """Log the number of FoVs and channels in the CZI images."""
    axes = images.axes

    if "S" not in axes:
        nFoVs = 1
        print("Found 1 field of view.")
    else:
        Si = axes.index("S")
        nFoVs = images.shape[Si]
        print("Found %d fields of view." % (nFoVs))

    return nFoVs


def log_axes(pixels, axes):
    """Log image shape."""
    print("; ".join(["%s:%d" % (a, pixels.shape[axes.index(a)]) for a in axes]))


def get_channel_names(images):
    """Extracts channel names from CZI images."""
    channel_path = "Metadata/DisplaySetting/Channels/Channel/DyeName"
    channel_names = [
        x.text for x in ET.fromstring(images.metadata).findall(channel_path)
    ]
    channel_names = [x.replace(" ", "").lower() for x in channel_names]
    return channel_names


def get_resolution(images):
    """Extracts voxel sides from CZI images."""
    res_path = "Metadata/Scaling/Items/Distance"
    resolution = [
        x
        for x in ET.fromstring(images.metadata).findall(res_path)
        if x.attrib["Id"] in ["X", "Y", "Z"]
    ]
    resolution = dict([(x.attrib["Id"], float(x[0].text)) for x in resolution])
    return resolution


def squeeze_axes(pixels, axes, targets=None, skip=None):
    """Squeeze specified single-dimension axes.

    Args:
        pixels (np.ndarray): stacks.
        axes (str): stacks axis labels.
        targets (str): axes to squeeze.
        skip (str): axes not to squeeze.
    """
    assert_msg = "axes expected to be %s, %s instead." % (type(""), type(axes))
    assert type("") == type(axes), assert_msg
    axes = list(axes)

    if type(None) != type(targets):
        for axis in targets:
            if axis not in axes:
                continue
            pixels = np.squeeze(pixels, axes.index(axis))
            axes.pop(axes.index(axis))

    if type(None) != type(skip):
        for axis in axes:
            if axis in skip:
                continue
            pixels = np.squeeze(pixels, axes.index(axis))
            axes.pop(axes.index(axis))

    axes = "".join(axes)
    return (pixels, axes)


def reorder_axes(pixels, axes, target):
    """Reorder stacks axes.

    Args:
        pixels (np.ndarray): stacks.
        axes (str): stacks axis labels.
        target (str): target axis order.
    """

    assert_msg = "[axes] %s expected, %s instead." % (type(""), type(axes))
    assert type("") == type(axes), assert_msg

    assert_msg = "[target] %s expected, %s instead." % (type(""), type(target))
    assert type("") == type(target), assert_msg

    if axes == target:
        return (pixels, target)

    axes = list(axes)
    target = list(target)

    assert len(target) == len(axes)
    for a in axes:
        assert a in target
    for a in target:
        assert a in axes

    target_positions = [target.index(a) for a in axes]
    pixels = np.moveaxis(pixels, range(len(axes)), target_positions)

    return (pixels, "".join(target))


def select_fov(pixels, si, mode="GPSeq"):
    """Generator that yields one channel stack with output path at a time.

    Args:
        pixels (np.ndarray): stacks.
        si (int): ID of the current Field of View (only for output porpuses).
        mode (str): output path notation.

    """

    assert mode in output_modes

    for ci in range(pixels.shape[0]):
        stack = pixels[ci]

        # Identify ouytput file name notation
        if mode == "GPSeq":
            outpath = "%s.channel%03d.series%03d.tif" % (
                channel_names[ci],
                ci + 1,
                si + 1,
            )
        elif mode == "DOTTER":
            outpath = "%s_%03d.tif" % (channel_names[ci], si + 1)

        yield ((stack, outpath))


def run():

    # PARAMETERS ===================================================================

    # Add script description
    parser = argparse.ArgumentParser(
        description="""
    Convert a czi file into single channel tiff images.
    Output file name is either in GPSeq (default) or DOTTER notation.
    Channel names are lower-cased.

    DOTTER:  dapi_001.tif
            <channelName>_<seriesNumber>.tif

    GPSeq:   dapi.channel001.series001.tif
            <channelName>.channel<channelNumber>.series<seriesNumber>.tif
    """,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # Add mandatory arguments
    parser.add_argument("input", type=str, help="""Path to the czi file to convert.""")

    # Add arguments with default value
    parser.add_argument(
        "-o",
        "--outdir",
        metavar="outdir",
        type=str,
        help="""Path to output TIFF folder, created if missing. Default to a
        folder with the input file basename.""",
        default=None,
    )
    output_modes = ("DOTTER", "GPSeq")
    parser.add_argument(
        "-m",
        "--mode",
        type=str,
        choices=output_modes,
        metavar="mode",
        help="""Output filename notation. Default: GPSeq.""",
        default="GPSeq",
    )

    # Add flags
    parser.add_argument(
        "--compressed",
        action="store_const",
        dest="doCompress",
        const=True,
        default=False,
        help="Force compressed TIFF as output.",
    )

    # Version flag
    version = "0.0.1"
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
    if type(None) == type(args.outdir):
        args.outdir = os.path.splitext(os.path.basename(args.input))[0]
        args.outdir = os.path.join(os.path.dirname(args.input), args.outdir)
        print("Output directory: '%s'" % args.outdir)

    # Additional checks
    assert os.path.isfile(args.input), "input file not found: %s" % args.input
    assert_msg = "output directory cannot be a file: %s" % args.outdir
    assert not os.path.isfile(args.outdir), assert_msg
    if not os.path.isdir(args.outdir):
        os.mkdir(args.outdir)

    # RUN ==========================================================================

    # Create buffer pointer to czi image
    images = czifile.CziFile(args.input)
    nFoVs = log_samples(images)

    with warnings.catch_warnings(record=True) as wlist:
        pixels = images.asarray()
    axes = images.axes

    if "T" in axes:
        assert_msg = "time-lapse images not supported."
        assert 1 == pixels.shape[axes.index("T")], assert_msg

    nChannels = pixels.shape[axes.index("C")]
    channel_names = get_channel_names(images)
    assert len(channel_names) == nChannels, "channel mismatch."

    resolution = get_resolution(images)
    resolutionXY = (1e-6 / resolution["X"], 1e-6 / resolution["Y"])

    pixels, axes = squeeze_axes(pixels, axes, skip="SCZYX")

    log_axes(pixels, axes)

    if not "S" in axes:
        # Single field of view
        axes_target = "CZYX"
        pixels, axes = reorder_axes(pixels, axes, axes_target)

        def fovGenerator():
            yield from select_fov(pixels, 1, mode=args.mode)

    else:
        axes_target = "SCZYX"
        pixels, axes = reorder_axes(pixels, axes, axes_target)

        def fovGenerator():
            for si in range(nFoVs):
                yield from select_fov(pixels[si, :], si, mode=args.mode)

    for (stack, outpath) in tqdm(fovGenerator(), total=nFoVs * nChannels):
        plot.save_tif(
            os.path.join(args.outdir, outpath),
            stack,
            imt.get_dtype(stack.max()),
            args.doCompress,
            bundled_axes="ZYX",
            resolution=resolutionXY,
            inMicrons=True,
            ResolutionZ=resolution["Z"] * 1e6,
            forImageJ=True,
        )

    ################################################################################

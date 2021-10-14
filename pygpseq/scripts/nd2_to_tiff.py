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
# Date: 20180316
# Project: bioimaging
# Description: convert nd2 images to tiff.
#
# ------------------------------------------------------------------------------


# DEPENDENCIES =================================================================

import argparse
from nd2reader import ND2Reader
from nd2reader.parser import Parser as ND2Parser
import numpy as np
import os
import sys
import tifffile
from tqdm import tqdm

from pygpseq.tools import image as imt
from pygpseq.tools import plot
from pygpseq.tools.io import printout


# FUNCTIONS ====================================================================


def getOutpath(args, metadata, cid, fid):
    outpath = None
    # Identify ouytput file name notation
    if args.mode == "GPSeq":
        outpath = "%s.channel%03d.series%03d.tif" % (
            metadata["channels"][cid].lower(),
            cid + 1,
            fid + 1,
        )
    elif args.mode == "DOTTER":
        outpath = "%s_%03d.tif" % (metadata["channels"][cid].lower(), fid + 1)
    return outpath


def export_channel(args, ch, outpath, metadata, bundled_axes, resolutionZ=None):
    resolutionXY = (1 / metadata["pixel_microns"], 1 / metadata["pixel_microns"])

    plot.save_tif(
        os.path.join(args.outdir, outpath),
        ch,
        imt.get_dtype(ch.max()),
        args.doCompress,
        bundled_axes=bundled_axes.upper(),
        resolution=resolutionXY,
        inMicrons=True,
        forImageJ=True,
        ResolutionZ=resolutionZ,
    )


def export_fov_3d(fov, metadata, fid, bundled_axes):
    """Export a field of view after bundling the axes to separate channel TIFF.

    Args:
        fov (np.ndarray): input multi-channel field of view array.
        metadata (dict): nd2 metadata dictionary.
        fid (int): field of view 0-based index.
    """

    # Get Z resolution
    if type(None) != type(args.deltaZ):
        resolutionZ = args.deltaZ
    else:
        with open(args.input, "rb") as fh:
            p = ND2Parser(fh)

            Zdata = np.array(p._raw_metadata.z_data)
            Zlevels = np.array(p.metadata["z_levels"]).astype("int")
            Zlevels = Zlevels + len(Zlevels) * fid
            Zdata = Zdata[Zlevels]

            resolutionZ = set(np.round(np.diff(Zdata), 3))

            assert_msg = "Z resolution is not constant: %s" % (str(resolutionZ))
            assert 1 == len(resolutionZ), assert_msg

            resolutionZ = list(resolutionZ)[0]

    if "c" in bundled_axes:
        # Iterate over channels
        for cid in range(fov.shape[3]):
            outpath = getOutpath(args, metadata, cid, fid)
            if type(None) == type(outpath):
                return ()
            ch = fov[:, :, :, cid]
            export_channel(args, ch, outpath, metadata, bundled_axes, resolutionZ)
    else:
        outpath = getOutpath(args, metadata, 0, fid)
        if type(None) == type(outpath):
            return ()
        export_channel(args, fov, outpath, metadata, bundled_axes, resolutionZ)


def export_fov_2d(fov, metadata, fid, bundled_axes):
    """Export a field of view after bundling the axes to separate channel TIFF.

    Args:
        fov (np.ndarray): input multi-channel field of view array.
        metadata (dict): nd2 metadata dictionary.
        fid (int): field of view 0-based index.
    """

    if "c" in bundled_axes:
        # Iterate over channels
        for cid in range(fov.shape[3]):
            outpath = getOutpath(args, metadata, cid, fid)
            if type(None) == type(outpath):
                return ()
            ch = fov[:, :, cid]
            export_channel(args, ch, outpath, metadata, bundled_axes)
    else:
        outpath = getOutpath(args, metadata, 0, fid)
        if type(None) == type(outpath):
            return ()
        export_channel(args, fov, outpath, metadata, bundled_axes)


def run():

    # PARAMETERS ===================================================================

    # Add script description
    parser = argparse.ArgumentParser(
        description="""
    Convert a nd2 file into single channel tiff images.
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
    parser.add_argument("input", type=str, help="""Path to the nd2 file to convert.""")

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
    parser.add_argument(
        "-m",
        "--mode",
        type=str,
        choices=("DOTTER", "GPSeq"),
        metavar="mode",
        help="""Output filename notation. Default: GPSeq.""",
        default="GPSeq",
    )
    parser.add_argument(
        "-Z",
        "--deltaZ",
        type=float,
        metavar="dZ",
        help="""If provided (in um), the script does not check delta Z consistency
        and instead uses the provided one.""",
        default=None,
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
    version = "2.1.0"
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
    if not os.path.isfile(args.input):
        printout("Input file not found: %s" % args.input, -2)
    if os.path.isfile(args.outdir):
        printout(
            "The specified output directory cannot be a file. Path: %s"
            % (args.outdir,),
            -2,
        )
    if not os.path.isdir(args.outdir):
        os.mkdir(args.outdir)

    # RUN ==========================================================================

    if type(None) != type(args.deltaZ):
        print("Enforcing a deltaZ of %.3f um." % args.deltaZ)

    # Create buffer pointer to nd2 image
    images = ND2Reader(args.input)
    print(args)
    print(images)

    if "t" in images.axes:
        assert (
            images.sizes["t"] == 1
        ), "conversion of time-course images not implemented."

    if "z" in images.axes:  # 3D
        axes_for_bundling = "zyxc" if "c" in images.axes else "zyx"

        if "v" not in images.axes:
            if "c" in images.axes:
                print("Found 1 field of view, with %d channels." % (images.sizes["c"]))
            else:
                print("Found 1 field of view, with 1 channel.")

            # Export single field of view
            images.bundle_axes = axes_for_bundling
            export_fov_3d(images[0], images.metadata, 0, axes_for_bundling)
        else:
            if "c" in images.axes:
                print(
                    "Found %d field of view, with %d channels."
                    % (images.sizes["v"], images.sizes["c"])
                )
            else:
                print("Found %d field of view, with 1 channel." % images.sizes["v"])

            # Export multiple fields of view
            images.iter_axes = "v"
            images.bundle_axes = axes_for_bundling
            for fid in tqdm(range(images.sizes["v"])):
                fov = images[fid]
                export_fov_3d(fov, images.metadata, fid, axes_for_bundling)
    else:  # 2D
        axes_for_bundling = "yxc" if "c" in images.axes else "yx"

        if "v" not in images.axes:
            if "c" in images.axes:
                print("Found 1 field of view, with %d channels." % (images.sizes["c"]))
            else:
                print("Found 1 field of view, with 1 channel.")

            # Export single field of view
            images.bundle_axes = axes_for_bundling
            export_fov_2d(images[0], images.metadata, 0, axes_for_bundling)
        else:
            print(
                "Found %d fields of view, with %d channels each."
                % (images.sizes["v"], images.sizes["c"])
            )

            # Export multiple fields of view
            images.iter_axes = "v"
            images.bundle_axes = axes_for_bundling
            for fid in tqdm(range(images.sizes["v"])):
                fov = images[fid]
                export_fov_2d(fov, images.metadata, fid, axes_for_bundling)

    ################################################################################

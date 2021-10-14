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
# Description: split TIFF in smaller images.
#
# Changelog:
#  v1.0.0 - 20180316: converted from MATLAB to Python3 and merged with pygpseq.
#
# ------------------------------------------------------------------------------


# DEPENDENCIES =================================================================

import argparse
import configparser as cp
import numpy as np
import os
import sys
from tqdm import tqdm

from ggc.prompt import ask

from pygpseq.tools import image as imt
from pygpseq.tools import plot
from pygpseq.tools.io import printout


# FUNCTIONS ====================================================================
version = "2.0.0"


def tiff_XY_enlarge(img, increases):
    """Resize image by increasing it on X and Y.

    Args:
        img (np.ndarray): image to resize.
        increases (list): list of int increases in px. (XY)

    Returns:
        np.ndarray: resized image.
    """
    N = len(img.shape)

    new_shape = list(img.shape)
    new_shape[-1] += int(increases[0])
    new_shape[-2] += int(increases[1])

    new_img = img.copy()
    new_img = np.zeros(new_shape)

    new_img[np.ix_(*[range(img.shape[i]) for i in range(len(img.shape))])] = img
    return new_img


def calc_loss(img, sides, step):
    """Calculate how many pixel/voxels would be lost if the image was split
    without enlarging.

    Args:
        img (np.ndarray): image to split.
        sides (list): [column, row] side in px.
        step (float): step fraction

    Returns:
        int: number of pixels/voxels lost during split
    """
    N = len(img.shape)
    assert len(sides) <= N

    if step is None:
        missed = [img.shape[-i - 1] % sides[i] for i in range(len(sides))]
    else:
        assert len(sides) == len(step)
        missed = [img.shape[-i - 1] % sides[i] % step[i] for i in range(len(sides))]

    loss = []
    for i in range(len(sides)):
        otherd = [img.shape[j] for j in range(N) if N - i - 1 != j]
        otherd.append(missed[-i - 1])
        loss.append(np.prod(otherd))
    loss = int(np.sum(loss) - np.prod(img.shape[:-2]) * np.prod(missed))

    return (*missed, loss, loss / np.prod(img.shape) * 100)


def tiff_split(img, sides, step, inverted=False):
    """Split 2D image in sub-images of x_side x y_side x stack_depth.
    Output is saved to outdir with the suffix .subN,
    where N is the sub-image index.

    Args:
        img (np.ndarray): image to split.
        sides (list): [column, row] side in px.
        step (float): step fraction.
        inverted (bool): split top-to-bottom, left-to-right.

    Returns:
        generator of np.ndarray: output split images
    """

    if step is None:
        step = sides

    n = img.shape[-1] // sides[0] * img.shape[-2] // sides[1]
    print("Output %d images." % n)
    if n == 0:
        return

    ys = [y for y in range(0, img.shape[-2], step[1]) if y + sides[1] <= img.shape[-2]]
    xs = [x for x in range(0, img.shape[-1], step[0]) if x + sides[0] <= img.shape[-1]]

    if inverted:
        print("Image split top-to-bottom, left-to-right.")
        xy_gen = ((x, y) for x in xs for y in ys)
    else:
        print("Image split left-to-right, top-to-bottom.")
        xy_gen = ((x, y) for y in ys for x in xs)

    assert len(img.shape) in [2, 3]
    if len(img.shape) == 3:
        tsplit = lambda i, x, y, s: i[:, y : (y + s[1]), x : (x + s[0])]
    elif len(img.shape) == 2:
        tsplit = lambda i, x, y, s: i[y : (y + s[1]), x : (x + s[0])]

    with tqdm(range(n)) as pbar:
        for (x_start, y_start) in xy_gen:
            yield tsplit(img, x_start, y_start, sides)
            pbar.update(1)


def print_settings(args, sides, clear=True):
    """Show input settings, for confirmation.

    Args:
        args (Namespace): arguments parsed by argparse.
        clear (bool): clear screen before printing.
    """

    s = " # TIFF split v%s\n" % version
    s += """
        Input file :  %s
Output directory :  %s

            X side : %d
            Y side : %d

        Overlap : %r
            Step : %r

            Slice : %r
        Enlarge : %r
        Inverted : %r
    """ % (
        args.input,
        args.outdir,
        sides[0],
        sides[1],
        args.overlap,
        args.step,
        args.slice,
        args.enlarge,
        args.inverted,
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
    Split a TIFF image in smaller TIFF images of the specified side(s). If two
    different sides are provided, the smaller images will be rectangular. The
    first side corresponds to the X (columns) and the second to the Y (rows).
    By default, only one side is required, which is used by the script for both
    X and Y sides. In other words, square smaller images are produced.

    If the original dimensions are not multiples of the specified side, a portion
    of the image is lost, unless the --enlarge option is used. In that case,
    the smaller images generated from the image border will contain empty pixels.

    If the input image is a 3D stack, it will be split only on XY and the output
    images will have the same number of slices. Using the --slice option, it is
    possible to specify which slice to split (i.e., the output will be in 2D).
    Defaults to first slice (--slice 0).

    It is also possible to generate overlapping split images. This can be achieved
    by using either the -S or -O options (which cannot be used together). With the
    -S option, you can specify the step used when splitting the image, as a fraction
    of its sides or as an absolute number of pixels. With the -O option, you can
    specify the overlapping region between consecutive split images as a fraction of
    their sides or as absolute pixels. In other words, the
    options -S 0.9 and -O 0.1 yield the same result. It is possible to provide two
    values to -S and -O, to obtain different overlaps in X and Y.

    By default, split images are generated left-to-right, top-to-bottom, e.g.,
    1 2 3
    4 5 6
    7 8 9

    Use the option -I to generate them top-to-bottom, left-to-right, e.g.,
    1 4 7
    2 5 8
    3 6 9

    Examples:

    - Square images of 100x100 px
    tiff_split big_image.tif split_out_dir 100 -e

    - Rectangular images of 125x100 px
    tiff_split big_image.tif split_out_dir 100 125 -e

    - Square images of 100x100 px, overlapping for 10 px in X and Y
    tiff_split big_image.tif split_out_dir 100 -e -S 0.9
    tiff_split big_image.tif split_out_dir 100 -e -S 90
    tiff_split big_image.tif split_out_dir 100 -e -O 0.1
    tiff_split big_image.tif split_out_dir 100 -e -O 10

    - Square images of 100x100 px, overlapping for 10 px in X and 20 px in Y
    tiff_split big_image.tif split_out_dir 100 -e -S 0.9 0.8
    tiff_split big_image.tif split_out_dir 100 -e -S 90 80
    tiff_split big_image.tif split_out_dir 100 -e -O 0.1 0.2
    tiff_split big_image.tif split_out_dir 100 -e -O 10 20
    """,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument("input", type=str, help="""Path to the TIFF image to split.""")
    parser.add_argument(
        "outdir", type=str, help="""Path to output TIFF folder, created if missing."""
    )
    parser.add_argument(
        "side",
        type=int,
        nargs="+",
        help="""One or two (XY) sides,
        used to specify the smaller images dimensions.""",
    )

    parser.add_argument(
        "-S",
        "--step",
        metavar="step",
        type=float,
        nargs="+",
        help="""Step for splitting, defined as a fraction of the
        specified side(s).""",
    )
    parser.add_argument(
        "-O",
        "--overlap",
        metavar="overlap",
        type=float,
        help="""Overlap fraction of splitted images, defined as a fraction of the
        specified side(s).""",
        nargs="+",
    )
    parser.add_argument(
        "-s",
        "--slice",
        metavar="slice",
        type=int,
        help="""ID of slice to be extracted from Z-stacks, 1-indexed.""",
    )

    parser.add_argument(
        "-e",
        "--enlarge",
        action="store_const",
        dest="enlarge",
        const=True,
        default=False,
        help="Expand to avoid pixel loss.",
    )
    parser.add_argument(
        "-I",
        "--invert",
        action="store_const",
        dest="inverted",
        const=True,
        default=False,
        help="""Split top-to-bottom, left-to-right.""",
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
        "--version",
        action="version",
        version="%s %s"
        % (
            sys.argv[0],
            version,
        ),
    )

    args = parser.parse_args()

    assert os.path.isfile(args.input), "input file not found: %s" % args.input
    assert_msg = "output directory cannot be a file: %s" % (args.outdir)
    assert not os.path.isfile(args.outdir), assert_msg

    sides = (args.side[0], args.side[0]) if len(args.side) == 1 else args.side[:2]
    if args.slice is not None:
        assert args.slice > 0
    assert_msg = "-S and -O are incompatible"
    assert args.step is None or args.overlap is None, assert_msg

    relative_steps = True
    if args.step is not None:
        assert all(step > 0 for step in args.step)
        step_is_relative = all(step <= 1 for step in args.step)
        step_is_absolute = all(step > 1 for step in args.step)
        assert step_is_absolute or step_is_relative

        while len(args.step) < len(sides):
            args.step.append(args.step[0])
        args.step = args.step[: len(sides)]

        if step_is_absolute:
            relative_steps = False
            args.overlap = np.array(
                [sides[i] - args.step[i] for i in range(len(args.step))]
            ).astype("int")
        elif step_is_relative:
            args.overlap = [np.round(1 - s, 3) for s in args.step]

    if args.overlap is not None:
        assert all(overlap >= 0 for overlap in args.overlap)
        overlap_is_relative = all(overlap < 1 for overlap in args.overlap)
        overlap_is_absolute = all(overlap > 1 for overlap in args.overlap)
        assert overlap_is_absolute or overlap_is_relative

        while len(args.overlap) < len(sides):
            args.overlap.append(args.overlap[0])
        args.overlap = args.overlap[: len(sides)]

        if overlap_is_absolute:
            relative_steps = False
            args.step = np.array(
                [sides[i] - args.overlap[i] for i in range(len(args.overlap))]
            ).astype("int")
        elif overlap_is_relative:
            args.overlap = [np.round(1 - s, 3) for s in args.overlap]

    if (args.overlap is not None or args.step is not None) and relative_steps:
        args.step = np.array(
            [np.round(sides[i] * args.step[i]) for i in range(len(args.step))]
        ).astype("int")
        args.overlap = np.array(
            [np.round(sides[i] * args.overlap[i]) for i in range(len(args.overlap))]
        ).astype("int")

    if not os.path.isdir(args.outdir):
        os.mkdir(args.outdir)

    # RUN ==========================================================================

    # Show current settings
    settings_string = print_settings(args, sides)
    if not args.do_all:
        ask("Confirm settings and proceed?")

    config = cp.ConfigParser()
    config["MAIN"] = {"input": args.input, "outdir": args.outdir}
    config["ADVANCED"] = {
        "x_side": str(sides[0]),
        "y_side": str(sides[1]),
        "x_step": str(args.step[0]),
        "y_step": str(args.step[1]),
        "slice": str(args.slice),
        "enlarge": str(args.enlarge),
        "inverted": str(args.inverted),
    }
    with open(os.path.join(args.outdir, "config.ini"), "w") as CF:
        config.write(CF)

    # Read input image
    print("Reading input image...")
    img = imt.read_tiff(args.input, k=3)

    # Check image shape and select appropriate analysis method ---------------------

    if len(img.shape) == 3:
        print("3D stack found: %s" % str(img.shape))
        if args.slice is not None:
            print("Enforcing 2D split (extracting slice #%d only)." % args.slice)
            umes = "pixel"
            assert args.slice <= img.shape[0]
            img = img[args.slice - 1, :, :].copy()
        else:
            umes = "voxel"
    elif len(img.shape) == 2:
        print("2D image found: %s" % str(img.shape))
        umes = "pixel"
    else:
        printout("Cannot split a 1D image. File: %s" % args.input, -2)

    # Enlarge or calculate pixel loss ----------------------------------------------

    x_loss, y_loss, loss, perc_loss = calc_loss(img, sides, args.step)
    if args.enlarge:
        img = tiff_XY_enlarge(img, np.array(sides) - np.array((x_loss, y_loss)))
        print("Image enlarged to %s" % str(img.shape))
    else:
        print("%d %ss in X and %d %ss in Y are lost." % (x_loss, umes, y_loss, umes))
        print(
            "In total, %d %ss lost (%.2f%%). Use -e to avoid loss."
            % (loss, umes, perc_loss)
        )

    # Split image ------------------------------------------------------------------

    prefix = os.path.splitext(os.path.basename(args.input))[0]
    ext = os.path.splitext(os.path.basename(args.input))[1]

    for ic, isplit in enumerate(
        tiff_split(img, sides, args.step, args.inverted), start=1
    ):
        opath = os.path.join(args.outdir, "%s.sub%d%s" % (prefix, ic, ext))
        plot.save_tif(opath, isplit, imt.get_dtype(isplit.max()), False)
    # END ==========================================================================

    print("DONE")

    ################################################################################

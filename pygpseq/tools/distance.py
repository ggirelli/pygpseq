# -*- coding: utf-8 -*-

"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
@description: distance calculation library.
"""

# DEPENDENCIES =================================================================

import os
import sys

import numpy as np

from scipy.ndimage.filters import convolve
from scipy.ndimage.morphology import distance_transform_edt
from scipy.stats import norm

import pygpseq as gp
from pygpseq.tools import image as imt

# FUNCTIONS ====================================================================


def calc_lamina_distance(mask, aspect):
    """Calculate lamin distance. If 3D image, add empty slices to top and
    bottom of the stack.

    Args:
        mask (np.ndarray): binary image.
        aspect (tuple[float]): pixel/voxel dimension proportion.
    """

    # Check for binary image
    assert_msg = "binary image expected, instead got [%f, %f]." % (
        mask.min(),
        mask.max(),
    )
    assert mask.max() <= 1 and mask.min() >= 0, assert_msg

    # Calculate lamin distance
    if len(mask.shape) != 3:
        return distance_transform_edt(mask, aspect[3 - len(mask.shape) :])
    # Add top/bottom empty slices
    zero_slice = np.zeros((1, mask.shape[1], mask.shape[2]))
    return distance_transform_edt(
        imt.add_top_bottom_slides(mask), aspect[3 - len(mask.shape) :]
    )[1:-1, :, :]


def calc_center_distance(laminD, aspect, asPercentile=False):
    """Calculate center distance. Center is by default defined as the voxel(s)
    with the highest lamin distance, otherwise as the top percentile.

    Args:
        mask (np.ndarray): binary image.
        aspect (tuple[float]): pixel/voxel dimension proportion.
        asPercentile (bool): define center as percentile.
    """

    return (
        distance_transform_edt(
            laminD < np.percentile(laminD, 99.0), aspect[3 - len(laminD.shape) :]
        )
        if asPercentile
        else distance_transform_edt(
            laminD != laminD.max(), aspect[3 - len(laminD.shape) :]
        )
    )


def mkGaussianKernel(size, sigma):
    """Generate a 1D Gaussian kernel of given size and sigma."""
    import numpy as np
    import scipy.stats as st

    assert 1 == size % 2
    D = np.array(range(size)) - size // 2
    mu = 0
    kernel = norm(mu, sigma).pdf(D)
    kernel /= kernel.sum()
    return kernel


def ndGuassianSmooth(T, sigma, aspect, normalized=False):
    """Perform iterative 1-D Gaussian convolution over each dimension of T."""
    radius = round(4 * sigma + 2)
    if 0 == radius % 2:
        radius += 1
    d = int(2 * radius + 1)
    V = T.copy()
    fill = V.max()
    for di in range(len(V.shape)):
        new_shape = np.ones(len(V.shape)).astype("i")
        new_shape[di] = d
        k = np.reshape(mkGaussianKernel(d, sigma * aspect[di] / aspect[0]), new_shape)
        V = convolve(V, k, mode="constant", cval=fill)
        if normalized:
            V /= convolve(np.ones(V.shape), k, mode="constant", cval=fill)
    return V


def quick_normalize(m):
    m = m.copy() - np.nanmin(m)
    return m / np.nanmax(m)


def simulate_diffusion(mask, sigma, aspect, simthr=0.7):
    """Simulates enzyme diffusion with constant external concentration and
    iterative anisotropic Gaussian blurring."""

    timebox = np.zeros(mask.shape)
    timebox[mask == 1] = -np.inf
    reached = 0 == mask
    simbox = (1 - mask).astype(np.float64)

    iterc = 1
    # outside = np.absolute((mask - 1).sum())
    # to_be_reached = np.prod(mask.shape) - outside
    while np.isinf(timebox.sum()):
        simbox = ndGuassianSmooth(simbox, sigma, aspect, True)
        simbox[mask == 0] = 1
        cond = simbox >= simthr
        timebox[np.logical_and(cond, np.logical_not(reached))] = iterc
        reached = np.logical_or(cond, reached)
        # flag = "%.2f%%" % (((cond).sum() - outside) / to_be_reached * 100.)
        # print((iterc, flag))
        iterc += 1

    timebox[np.isinf(timebox)] = np.nan

    return timebox


def calc_nuclear_distances(dist_type, mask, aspect):
    """Calculate distance from lamina and center for each voxel in a nucleus.

    Args:
        dist_type (str): any string from gp.const.LD_ARG_LABELS

    Returns:
        (np.ndarray, np.ndarray): lamina_distance, center_distance
    """

    assert_msg = "expected one of %s, got '%s'" % (
        str(gp.const.LD_ARG_LABELS),
        dist_type,
    )
    assert dist_type in range(len(gp.const.LD_ARG_LABELS)), assert_msg

    if dist_type == gp.const.LD_DIFFUSION:
        laminD = simulate_diffusion(mask, 1, aspect)
        centrD = np.absolute(laminD - np.nanmax(laminD))
    else:
        laminD = calc_lamina_distance(mask, aspect)
        centrD = calc_center_distance(
            laminD, aspect, dist_type == gp.const.LD_CENTER_PERC
        )
    return (laminD, centrD)


def normalize_nuclear_distance(dist_type, laminD, centrD):
    """Normalize lamina distnace."""

    return (
        quick_normalize(laminD)
        if dist_type == gp.const.LD_DIFFUSION
        else laminD / (laminD + centrD)
    )


# END ==========================================================================

################################################################################

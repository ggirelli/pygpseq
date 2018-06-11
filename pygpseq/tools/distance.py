# -*- coding: utf-8 -*-

'''
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
@description: distance calculation library.
'''

# DEPENDENCIES =================================================================

import os
import sys

import numpy as np

# FUNCTIONS ====================================================================

def calc_lamina_distance(mask, aspect):
    '''Calculate lamin distance. If 3D image, add empty slices to top and
    bottom of the stack.

    Args:
        mask (np.ndarray): binary image.
        aspect (tuple[float]): pixel/voxel dimension proportion.
    '''

    # Check for binary image
    assert_msg = "binary image expected, instead got [%f, %f]." % (
        mask.min(), mask.max())
    assert 1 >= mask.max() and 0 <= mask.min(), assert_msg

    # Calculate lamin distance
    if 3 == len(mask.shape):
        # Add top/bottom empty slices
        zero_slice = np.zeros((1, mask.shape[1], mask.shape[2]))
        laminD = distance_transform_edt(
            np.concatenate([zero_slice, mask, zero_slice]),
            aspect[3-len(mask.shape):])[1:-1, :, :]
    else:
        laminD = distance_transform_edt(mask,
            aspect[3-len(mask.shape):])
    return(laminD)

def calc_center_distance(laminD, aspect, asPercentile = False):
    '''Calculate center distance. Center is by default defined as the voxel(s)
    with the highest lamin distance, otherwise as the top percentile.

    Args:
        mask (np.ndarray): binary image.
        aspect (tuple[float]): pixel/voxel dimension proportion.
        asPercentile (bool): define center as percentile.
    '''

    # Center as top percentile
    if asPercentile:
        centrD = distance_transform_edt(
            laminD < np.percentile(laminD, 99.),
            aspect[3-len(laminD.shape):])

    # Center as top value
    else:
        centrD = distance_transform_edt(laminD != laminD.max(),
            aspect[3-len(laminD.shape):])

    return(centrD)

def mkGaussianKernel(size, sigma):
    """Generate a 1D Gaussian kernel of given size and sigma."""
    import numpy as np
    import scipy.stats as st
    assert 1 == size % 2
    D = np.array(range(size)) - size // 2
    mu = 0
    kernel = st.norm(mu, sigma).pdf(D)
    kernel /= kernel.sum()
    return kernel

def ndGuassianSmooth(T, sigma, aspect, normalized = False):
    """Perform iterative 1-D Gaussian convolution over each dimension of T."""
    radius = round(4 * sigma + 2)
    if 0 == radius % 2:
        radius += 1
    d = int(2 * radius + 1)
    V = T.copy()
    for di in range(len(V.shape)):
        new_shape = np.ones(len(V.shape)).astype('i')
        new_shape[di] = d
        k = np.reshape(mkGaussianKernel(d, sigma * aspect[di] / aspect[0]), new_shape)
        V = convolve(V, k, mode = 'constant')
        if normalized:
            V /= convolve(np.ones(V.shape), k, mode = 'constant')
    return V

def simulate_diffusion(mask, sigma, aspect, simthr = .7):
    """Simulates enzyme diffusion with constant external concentration and
    iterative anisotropic Gaussian blurring."""
    timebox = np.zeros(mask.shape)
    timebox[1 == mask] = -np.inf
    reached = 0 == mask
    simbox = 1 - mask

    iterc = 1
    while np.isinf(timebox.sum()):
        simbox = ndGuassianSmooth(simbox, sigma, aspect, True)
        simbox[mask == 0] = 1
        cond = simbox >= simthr
        timebox[np.logical_and(cond, np.logical_not(reached))] = iterc
        reached = np.logical_or(cond, reached)
        iterc += 1

    timebox[np.isinf(timebox)] = np.nan

    return timebox

# END ==========================================================================

################################################################################

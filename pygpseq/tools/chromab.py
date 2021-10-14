# -*- coding: utf-8 -*-

"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
@description: chromatic aberration correction library.
"""

# DEPENDENCIES =================================================================

import numpy as np
from skimage import filters
from scipy.interpolate import griddata
from skimage.morphology import dilation

# FUNCTIONS ====================================================================


def correct_stack(stack, chname, m):
    """Apply chromatic aberraction correction to a 3D image (stack),
    using the specified channel as a reference.

    Args:
      stack (ndarray)
      chname (string): name of the stack's channel
      m (dict): CA correction measurements {cnames, P:ndarray, N:int}}
    """

    # Identify coefficients of interes
    ids = filter(lambda i: m["chan"][i] == chname, range(len(m["chan"])))
    dz = m["dz"][ids]
    Cx = m["Cx"][:, ids]
    Cy = m["Cy"][:, ids]

    # Setup transformation matrices
    # ==================================

    # Coordinates
    Y, X = np.meshgrid(range(stack.shape[1]), range(stack.shape[2]))
    X = X.reshape((np.prod(X.shape), 1))
    Y = Y.reshape((np.prod(Y.shape), 1))

    # Polynomial
    PD = poly2mat(np.concatenate((X, Y), 1), 3)
    QDx = PD.dot(Cx).reshape((len(X),))
    QDy = PD.dot(Cy).reshape((len(Y),))

    # Correct every single plane
    jstack = stack.copy()
    for i in range(stack.shape[0]):
        rstack = stack[i, :, :].reshape((np.prod(stack.shape[1:3]),))
        jstack[i, :, :] = griddata(
            (Y.reshape((len(Y),)), X.reshape((len(X),))),
            rstack,
            np.transpose(np.array([QDy, QDx])),
            method="cubic",
        ).reshape(stack[i, :, :].shape)

    return jstack


def dot_candidates(I, sigma=None, sigma_diff=None, max_n_dots=None, padding=None):
    """Extract dots from the stack I as local maximas and order them
    according to the DoG (difference of Gaussians).

    Note:
      Function built ONLY for 3D images.

    Args:
      I (ndarray): 3D image (stack)
      sigma (float): sigma for DoG calculation (opt, def 1.2)
      sigma_diff (float): sigmas distance for DoG (opt, def 1e-2)
      max_n_dots (int): maximum number of dots to extract (opt, def 1e4)
      padding (tuple[int]): padding to clear borders ZYX (opt, def (1,5,5))

    Returns:
      np.array: a record array with the following fields:
        dog (f): DoG value for that maxima
        int (f): intensity for that maxima (actually, I.dtype)
        x, y, z (u4): coordinates in I
    """

    # Chek image dimensionality
    if len(I.shape) != 3:
        return I

    # Default values
    if sigma is None:
        sigma = 1.2
    if sigma_diff is None:
        sigma_diff = 1e-2
    if max_n_dots is None:
        max_n_dots = 1e4
    if padding is None:
        padding = (1, 5, 5)

    # Just look at the intensity
    J = I.copy().astype("float")
    if J.shape[0] > 1:
        A = dilation(J, np.ones((3, 3, 3)))
    else:
        A = dilation(J, np.ones((3, 3)))

    # Clear borders
    A = set_borders(J, padding, -1)

    # Don't consider saturated pixels
    A[A == np.iinfo(I.dtype).max] = -1

    # Identify maxima
    Pos = np.where(J == A)

    # Difference of Gaussians
    DoG = np.zeros(I.shape)
    for i in range(I.shape[0]):
        DoG[:, :, i] = filters.gaussian(J[:, :, i], sigma)
        DoG[:, :, i] -= filters.gaussian(J[:, :, i], sigma + sigma_diff)

    # Prepare output
    P = np.zeros(
        (len(Pos[0]),),
        dtype=[
            ("dog", "float"),
            ("int", "float"),
            ("x", "u4"),
            ("y", "u4"),
            ("z", "u4"),
        ],
    )
    P["dog"] = DoG[Pos]
    P["int"] = I[Pos]
    P["x"] = Pos[0]
    P["y"] = Pos[1]
    P["z"] = Pos[2]

    # Sort output
    P.sort(order="dog")

    # Subselect output
    P = P[range(max_n_dots)]

    # Output
    return P


def dot_fitting(I, P, cluster_min_dist, figmas, use_clustering=None):
    """Maximum-likelihood localization of points P in image I.

    Args:
      I (ndarray): 3D image (stack)
      P (ndarray): table with extracted dots
      cluster_min_dist (int): shortest distance for clustering
      sigmas ({dim[float]}): size of Gaussian profile to be fitted [x,y,z]
      use_clustering (bool): True to use clustering
    """


def poly2mat(m, i):
    """Returns the i-th polynomial matrix of m."""

    if i == 3:
        out = np.concatenate(
            (
                np.ones((m.shape[0], 1)),
                m[:, [0]],
                m[:, [1]],
                m[:, [0]] ** 2,
                m[:, [0]] * m[:, [1]],
                m[:, [1]] ** 2,
                m[:, [0]] ** 3,
                m[:, [0]] ** 2 * m[:, [1]],
                m[:, [0]] * m[:, [1]] ** 2,
                m[:, [1]] ** 3,
            ),
            1,
        )
    elif i == 2:

        out = np.concatenate(
            (
                np.ones((m.shape[0], 1)),
                m[:, [0]],
                m[:, [1]],
                m[:, [0]] ** 2,
                m[:, [0]] * m[:, [1]],
                m[:, [1]] ** 2,
            ),
            1,
        )
    elif i == 1:
        out = np.ones((m.shape[0], 3))
        out[:, 1:3] = m[:, :2]
    else:
        out = m
    return out


def set_borders(I, paddings=None, value=None):
    """Sets the edge pixels to the provided value.

    Args
      I (np.array): 3D image
      paddings (tuple[int]): padding for every dimension (opt, def 0)
      value (I.dtype): border pixels will be set to this value (opt, def 0)
    """

    # Avoid in-place change
    J = I.copy()

    # Default values
    if paddings is None:
        paddings = np.zeros(len(I.shape)).tolist()
    if value is None:
        value = 0

    # Set padding for missing dimensions
    paddings = list(paddings)
    paddings.extend(np.zeros(len(I.shape) - len(paddings)).tolist())

    # Set values
    for j in range(len(I.shape)):
        p = int(paddings[j])

        # Right
        right_pad = [
            tuple(range(I.shape[i])) if i != j else tuple(range(p))
            for i in range(len(I.shape))
        ]
        J[np.ix_(*right_pad)] = value

        # Left
        left_pad = [
            tuple(range(I.shape[i]))
            if i != j
            else tuple(range(I.shape[i] - p, I.shape[i]))
            for i in range(len(I.shape))
        ]
        J[np.ix_(*left_pad)] = value

    return J


# END ==========================================================================

################################################################################

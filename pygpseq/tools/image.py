# -*- coding: utf-8 -*-

"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
@description: image manipulation library.
"""

# DEPENDENCIES =================================================================

import os
import sys

import numpy as np
from scipy import ndimage as ndi
from scipy.ndimage.morphology import distance_transform_edt
from skimage import filters
from skimage.io import imread
from skimage.measure import label, marching_cubes_lewiner, mesh_surface_area
from skimage.morphology import closing, convex_hull_image, cube
from skimage.morphology import dilation, erosion, square
from skimage.segmentation import clear_border
import warnings

from pygpseq import const
from pygpseq.tools import vector as vt
from pygpseq.tools.io import printout

# FUNCTIONS ====================================================================


def add_top_bottom_slides(i):
    zero_slice = np.zeros((1, i.shape[1], i.shape[2]))
    return np.concatenate([zero_slice, i, zero_slice])


def apply_box(i, box):
    """Apply square/box selection to an image.

    Args:
      i (np.array): image.
      box (list): selection square/box corner coordinates.

    Returns:
      np.array: square/box selection of the provided image.
    """

    # Check box
    box = check_box(i.shape, box)
    # Apply box
    return i[np.ix_(*[tuple(range(t[0], t[1] + 1)) for t in box])]


def autoselect_time_frame(im):
    """Selects the first non-empty time frame found.

    Args:
      im (np.array): image.
    """

    if len(im.shape) != 4:
        return im
    if im.shape[0] == 1:
        return im

    selected = None

    for i in range(im.shape[3]):
        if im[:, :, :, i].max() != 0:
            selected = i
            break

    return im[:, :, :, selected]


def binarize(i, thr):
    """Binarize an image using the provided threshold.

    Args:
      i (np.array): image.
      thr (float or int): intensity threshold.

    Returns:
      np.array: thresholded image.
    """

    if len(i.shape) == 2:
        i = closing(i > thr, square(3))
    elif len(i.shape) == 3:
        i = closing(i > thr, cube(3))
    return i


def calc_surface(mask, spacing=None):
    """Calculate the surface of a binary mask.
    The provided mask is expected to have only one object.

    Args:
      mask (np.array): thresholded image.
      spacing (tuple[float]): pixel/voxel side sizes.

    Returns:
      float: mesh surface area of the provided object.
    """

    # Aspect ratio for 3d surface calculation
    if spacing is None:
        spacing = [1.0 for d in mask.shape]

    # Force binary type
    mask = mask.astype("bool")

    # Check number of objects
    if label(mask).max() != 1:
        return 0

    # Add top/bottom slices
    shape = list(mask.shape)
    shape[0] = 1
    mask = np.vstack((np.zeros(shape), mask, np.zeros(shape)))

    # Calculate sphericity
    verts, faces, ns, vs = marching_cubes_lewiner(mask, 0.0, spacing)
    return mesh_surface_area(verts, faces)


def check_box(shape, box):
    """Check if a square/box selection can be applied to an image.
    If not enough corners are specified, the whole dimension is selected.

    Args:
      shape (tuple[int]): image shape.
      box (list): selection square/box corner coordinates.
    """

    # Every dimension in the box should have two coordinates
    if len(shape) >= len(box):
        for i in range(len(box)):
            if len(box[i]) != 2:
                box[i] = (0, shape[i] - 1)
    else:
        box = box[: len(shape)]

    # Fill missing dimensions
    if len(shape) > len(box):
        [box.extend([0, shape[d]]) for d in range(int(len(box) / 2), len(shape))]

    # Remove extra dimensions
    box = box[: len(shape)]

    return box


def clear_borders(img, clean_z=None):
    """Remove objects touching the borders of the image.

    Args:
      img (np.array): binary image.
      clean_z (bool): True to remove the objects touching the Z borders.

    Returns:
      np.array: cleaned image.
    """

    if len(img.shape) == 2:
        img = clear_border(img)
    elif len(img.shape) == 3:
        for slide_id in range(img.shape[0]):
            img[slide_id, :, :] = clear_border(img[slide_id, :, :])
        if clean_z == True:
            for slide_id in range(img.shape[1]):
                img[:, slide_id, :] = clear_border(img[:, slide_id, :])
    return img


def clear_borders2(img, clean_z=None):
    """Remove objects touching the borders of the image.

    Args:
      img (np.array): labeled image.
      clean_z (bool): True to remove the objects touching the Z borders.

    Returns:
      np.array: cleaned image.
    """

    if len(img.shape) == 2:
        img = clear_border(img)
    elif len(img.shape) == 3:
        # Identify 3D objects touching X/Y borders
        l = []
        l.extend(np.unique(img[:, 0, :]).tolist())
        l.extend(np.unique(img[:, -1, :]).tolist())
        l.extend(np.unique(img[:, :, 0]).tolist())
        l.extend(np.unique(img[:, :, -1]).tolist())
        ii = set(l)

        # Remove them
        for i in ii:
            img[img == i] = 0

        # Apply also to Z borders
        if clean_z == True:
            # Identify
            l = []
            l.extend(np.unique(img[0, :, :]).tolist())
            l.extend(np.unique(img[-1, :, :]).tolist())
            ii = set(l)

            # Remove
            for i in ii:
                img[img == i] = 0

        # Re-label
        img = label(img)

    # Output
    return img


def describe_shape(mask, spacing=None):
    """Calculate sphericity (3d) or solidity (2d) of the provided mask.
    The provided mask is expected to have only one object.

    Args:
      mask (np.array): thresholded image.

    Returns:
      float: shape descriptor of the provided object.
    """

    # Aspect ratio for 3d surface calculation
    if spacing is None:
        spacing = [1.0 for d in mask.shape]

    # Force binary type
    mask = mask.astype("bool")

    # Check number of objects
    if label(mask).max() != 1:
        return 0

    # Calculate shape descriptor
    if len(mask.shape) == 2:
        # Calculate solidity
        return float(mask.sum()) / convex_hull_image(mask).sum()
    elif len(mask.shape) == 3:
        # Add top/bottom slices
        shape = list(mask.shape)
        shape[0] = 1
        mask = np.vstack((np.zeros(shape), mask, np.zeros(shape)))

        # Calculate sphericity
        verts, faces, ns, vs = marching_cubes_lewiner(mask, 0.0, spacing)
        s = mesh_surface_area(verts, faces)
        return (np.pi * (6.0 * mask.sum()) ** 2) ** (1 / 3.0) / s
    else:
        return 0


def dilate_fill_erode(mask, strel):
    """Performs dilation-fill-erosion of mask with the provided structuring
    element."""

    assert_msg = "structuring element and mask must have the same dimensions."
    assert len(mask.shape) == len(strel.shape), assert_msg

    mask = dilation(mask, strel)
    mask = fill_holes(mask)
    mask = erosion(mask, strel)

    return mask


def estimate_background(i, mask, seg_type):
    """Estimates background median.

    Args:
      i (np.array): image.
      seg_type (string): segmentation type as defined in pygpseq.const.

    Returns:
      float: estimated background.
    """

    if const.SEG_3D != seg_type and len(i.shape) != 2:
        return np.median(mk_z_projection(i, seg_type)[mask == 0])
    else:
        return np.median(i[mask == 0])


def fill_holes(mask):
    """Fill mask holes."""
    mask = ndi.binary_fill_holes(mask)
    if len(mask.shape) == 3:
        for sliceid in range(mask.shape[0]):
            slide = mask[sliceid, :, :]
            mask[sliceid, :, :] = ndi.binary_fill_holes(slide)
    return mask


def get_dtype(i):
    """
    Identify bit depth for a matrix of maximum intensity i.
    """
    depths = [8, 16]
    for depth in depths:
        if i <= 2 ** depth - 1:
            return "uint%d" % (depth,)
    return "uint"


def get_mid_section_idx(i, mask, mid_type=None):
    """Identify mid-section index.

    Note:
      Possible mid-section definitions:
        const.MID_SEC_CENTRAL : middle section
        const.MID_SEC_LARGEST : largest (size) section
        const.MID_SEC_MAXSUMI : section with maximum DNA intensity sum
        const.MID_SEC_DEFAULT : const.MID_SEC_LARGEST

    Args:
      i (np.array): image (DNA channel).
      mask (np.array): binarized version of i.
      mid_type (int): mid-section definition.

    Returns:
      int: mid-section index using provided mid-section definition.
    """

    # Back to default mid-section definition.
    if mid_type is None:
        mid_type = const.MID_SEC_DEFAULT

    # Need a stack, otherwise get the full image
    if len(i.shape) < 3:
        return 0

    # Select central slide
    if const.MID_SEC_CENTRAL == mid_type:
        return i.shape[0] / 2

    if mid_type in [const.MID_SEC_LARGEST, const.MID_SEC_DEFAULT]:
        idx = [mask[idx, :, :].sum() for idx in range(mask.shape[0])]
        return idx.index(max(idx))

    if mid_type in [const.MID_SEC_MAXSUMI]:
        idx = [i[idx, :, :].sum() for idx in range(mask.shape[0])]
        return idx.index(max(idx))


def get_objects_zsize(L):
    """Retrieve objects size (2/3D).

    Args:
      L (np.array): labelled thresholded image.

    Returns:
      list: Z size of every object in the labelled image.
    """

    return [(L == i).astype("int").sum(0).max() for i in range(1, L.max())]


def get_objects_xysize(L):
    """Retrieve objects size (2/3D).

    Args:
      L (np.array): labelled thresholded image.

    Returns:
      list: XY size of every object in the labelled image.
    """

    Larray = L.reshape([np.prod(L.shape)]).tolist()
    return [t[1] for t in vt.uniquec(Larray) if t[0] != 0]


def get_partial_nuclear_volume(mask, i, erosion):
    """Retrieve the partial volume of a nucleus.

    Note:
      Top half (from highest intensity sum slice).
      Central portion (through erosion).

    Args:
      mask (np.array): thresholded image.
      i (np.array): image.
      erosion (float): maximum distance allowed.

    Returns:
      sel_mask (np.array): partial volume mask.
      err_log (string): error log if anything bad occurred.
    """

    # Empty error log
    err_log = ""

    # Copy original mask
    sel_mask = mask.copy()

    assert 0 != mask.sum()

    # Identify middle section
    mid = [i[sliceid, :, :].sum() for sliceid in range(mask.shape[0])]
    mid = mid.index(max(mid))

    # Select top half
    sel_mask[0:mid, :, :] = False

    # Select erosion on max projection
    d2d = distance_transform_edt(sel_mask.max(0))

    # Normalize distance
    d2d = d2d - d2d.min()
    if d2d.max() != 0:
        d2d = d2d / float(d2d.max())
    else:
        # Log error
        err_log += "Found 0 maximum distance."

    # Make 2d mask
    mask2d = d2d >= erosion

    # Apply 2d erosion mask to every slice
    for i in range(sel_mask.shape[0]):
        sel_mask[i, :, :] = np.logical_and(sel_mask[i, :, :], mask2d)

    # Convert in unsigned integers
    sel_mask = sel_mask.astype("u4")

    # Output
    return (sel_mask, err_log)


def get_rescaling_factor(path, **kwargs):
    """Get rescaling factor for deconvolved.

    Args:
      path (string):
      **kwargs: basedir additional argument

    Returns:
      float: scaling factor
    """

    # Set basedir if not provided
    if "basedir" not in kwargs.keys():
        kwargs["basedir"] = os.path.dirname(path)

    # Build proper path to the deconvolution log file
    path = os.path.join(kwargs["basedir"], os.path.basename(path))
    path = list(os.path.splitext(path))
    path[0] = path[0] + "_history"
    path[-1] = ".txt"
    path = "".join(path)

    # Means that the image was not deconvolved
    if not os.path.exists(path):
        factor = 1
    else:
        # Identify line with scaling factor
        with open(path, "r") as fhistory:
            frows = fhistory.readlines()

        # Retrieve factor
        needle = "Stretched to Integer type"
        factor = [x for x in frows if needle in x]
        factor = 1 if not factor else float(factor[0].strip().split(" ")[-1])
    # Output
    return factor


def get_unit(shape):
    """Get size unity.

    Args:
      shape (tuple[int]): image shape.

    Returns:
      string: "vx" for 3d images, "px" for 2d images.
    """
    if len(shape) == 2:
        return "px"
    elif len(shape) == 3:
        return "vx"
    else:
        return ""


def in_3d_box(box, coords):
    """
    Check if point is in a box

    Args:
      box (tuple): ((x0, x1), (y0, y1), (z0, z1)).
      coords (tuple): (x, y, z).

    Returns
      bool
    """
    cx = coords[0] >= box[0][0] and coords[0] <= box[0][1]
    cy = coords[1] >= box[1][0] and coords[1] <= box[1][1]
    cz = coords[2] >= box[2][0] and coords[2] <= box[2][1]
    return cx and cy and cz


def in_mask(coords, imbin):
    """Check if a pixel in a binary mask is foreground.

    Args:
        coords (tuple): coordinates.
        imbin (np.ndarray): binary image.
    """

    # Check the pixel is inside the image boundaries
    inbound = imbin.shape[0] > coords[0]
    inbound = inbound and imbin.shape[1] > coords[1]
    inbound = inbound and imbin.shape[2] > coords[2]
    inbound = inbound and all(np.array(coords) >= 0)
    if not inbound:
        return False

    # Check the pixel is foreground
    return imbin[coords[0], coords[1], coords[2]] == 1


def mkIsoStruct(dilate_factor, aspect):
    """
    Builds isotropic structuring element for dilation.

    Args:
      dilate_factor (int): number of px for isotropic 2D dilation.
      aspect (tuple(float)): voxel side ratios.

    Returns:
      np.ndarray: structureing element for 3D anisotropic dilation.
    """

    # XY dilation factor
    df_xy = int(dilate_factor * 2 + 1)
    if aspect[0] == 0:
        se = cube(df_xy)
        se = se[0]
        new_shape = [1]
        [new_shape.append(d) for d in se.shape]
        se = se.reshape(new_shape)
        return se

    # Z dilation factor
    df_z = int(dilate_factor * aspect[1] / aspect[0] * 2 + 1)

    if df_z == df_xy:
        # Isotropic
        return cube(df_z)
    elif df_z > df_xy:
        # Larger Z side
        se = cube(df_z)
        se = se[:, 0:df_xy, 0:df_xy]
    else:
        # Larger XY side
        se = cube(df_xy)
        se = se[0:df_z]

    # Output
    return se


def mk_z_projection(i, p_type):
    """Make Z-projection.
    Allowed p_types: SUM and MAX. Otherwise return the original image.

    Args:
      i (np.array): image.
      p_type (string): Z-projection type according to pygpseq.const.

    Returns:
      np.array: Z-projection.
    """

    if const.SUM_PROJ == p_type:
        i = i.sum(0).astype(i.dtype)
    elif const.MAX_PROJ == p_type:
        i = i.max(0).astype(i.dtype)
    return i


def read_tiff(impath, k=None, noSelection=False, rescale=1):
    """Read tiff image.

    Args:
      impath (str): path to tiff image.
      k (int): number of dimensions in output image for re-slicing.
      noSelection (bool): whether to discard empty dimensions.
      rescale (float): scaling factor.

    Returns:
      np.ndarray: image.
      None: file is possibly corrupt.
    """

    assert os.path.isfile(impath), "trying to read missing file"

    # Read TIFF (capture any parsing issues)
    try:
        with warnings.catch_warnings(record=True) as wlist:
            im = imread(impath)
            if len(wlist) != 0 and "axes do not match shape" in str(wlist[0]):
                printout(
                    "image axes do not match metadata in '%s'. %s"
                    % (impath, "Using the image axes."),
                    -1,
                )
    except (ValueError, TypeError) as e:
        msg = "Something went wrong while trynig to read a file"
        printout("%s (possibly corrupt):\n%s\n" % (msg, impath), -2, canAbort=False)
        return None

    # Reshape and re-slice
    while im.shape[0] == 0 and not noSelection:
        im = im[0]
    if type(0) == type(k):
        im = slice_k_d_img(im, k)

    # Rescale
    if rescale != 1:
        im = (im / rescale).astype("float")

    return im


def rm_from_mask(L, torm):
    # Remove elements from a mask.
    #
    # Args:
    #     L (np.array[int]): labelled objects.
    #     torm (list): list of objects indexes (to remove).

    if len(torm) <= L.max() - len(torm):
        # Update list of objects to be discarded
        torm = [e + 1 for e in torm]

        # Identify which objects to discard
        rm_mask = np.vectorize(lambda x: x in torm)(L)

    else:
        # Select objects to be kept
        tokeep = [e + 1 for e in range(L.max()) if e not in torm]

        # Identify which objects to discard
        rm_mask = np.vectorize(lambda x: x not in tokeep)(L)

    # Discard and re-label
    L[rm_mask] = 0
    # Output
    return L > 0


def slice_k_d_img(img, k):
    """Select one k-d image from a n-d image.
    Note: n >= k

    Args:
      img (np.array): image.
      k (int): number of dimensions to keep.

    Returns:
      np.array: k-d image.
    """

    # Check that k is lower than or equal to n
    if k > len(img.shape):
        return img

    # Slice image
    idxs = [(0,) for i in range(len(img.shape) - k)]
    for i in range(k):
        idxs.append(tuple(range(img.shape[len(img.shape) - k + i])))
    img = img[np.ix_(*idxs)].reshape(img.shape[len(img.shape) - k :])

    # Output
    return img


def threshold_adaptive(
    i, block_size, doClosing=True, method=None, mode=None, param=None
):
    """Adaptive threshold.

    Args:
      i (np.array): image.
      block_size (int): neighbourhood, if even then it is incremented by 1.
      doClosing (bool): to trigger closing operation after local thresholding.
      method, mode, param: additional parameters for threshold_local.

    Returns:
      np.array: thresholded image.
    """

    # Increment neighbourhood size
    if 0 == block_size % 2:
        block_size += 1

    # Local threshold mask
    lmask = i.copy()

    # Apply threshold per slice
    if len(i.shape) == 3:
        for slice_id in range(i.shape[0]):
            lmask[slice_id, :, :] = filters.threshold_local(
                i[slice_id, :, :], block_size, method=method, mode=mode, param=param
            )
        lmask = i >= lmask
        if doClosing:
            lmask = closing(lmask, cube(3))
    else:
        lmask = filters.threshold_local(
            i, block_size, method=method, mode=mode, param=param
        )
        lmask = i >= lmask
        if doClosing:
            lmask = closing(lmask, square(3))

    # Output
    return lmask


# END ==========================================================================

################################################################################

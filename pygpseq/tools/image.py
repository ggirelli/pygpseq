# -*- coding: utf-8 -*-

'''
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
@description: image manipulation library.
'''

# DEPENDENCIES =================================================================

import os
import sys

import numpy as np
from scipy.ndimage.morphology import distance_transform_edt
from skimage import filters
from skimage.io import imread
from skimage.measure import label, marching_cubes_lewiner, mesh_surface_area
from skimage.morphology import closing, convex_hull_image, cube
from skimage.morphology import square
from skimage.segmentation import clear_border

from pygpseq import const
from pygpseq.tools import vector as vt
from pygpseq.tools.io import printout

# FUNCTIONS ====================================================================

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
    return(i[np.ix_(*[tuple(range(t[0], t[1]+1)) for t in box])])

def autoselect_time_frame(im):
    """Selects the first non-empty time frame found.

    Args:
      im (np.array): image.
    """

    if 4 == len(im.shape):
        if im.shape[0] == 1:
            return(im)
        
        selected = None

        if 4 == len(im.shape):
            for i in range(im.shape[3]):
                if 0 != im[:, :, :, i].max():
                    selected = i
                    break

        return(im[:, :, :, selected])
    else:
        return(im)

def binarize(i, thr):
    """Binarize an image using the provided threshold.

    Args:
      i (np.array): image.
      thr (float or int): intensity threshold.

    Returns:
      np.array: thresholded image.
    """

    if 2 == len(i.shape):
        i = closing(i > thr, square(3))
    elif 3 == len(i.shape):
        i = closing(i > thr, cube(3))
    return(i)

def calc_surface(mask, spacing = None):
    """Calculate the surface of a binary mask.
    The provided mask is expected to have only one object.

    Args:
      mask (np.array): thresholded image.
      spacing (tuple[float]): pixel/voxel side sizes.

    Returns:
      float: mesh surface area of the provided object.
    """

    # Aspect ratio for 3d surface calculation
    if None == spacing:
        spacing = [1.0 for d in mask.shape]

    # Force binary type
    mask = mask.astype('bool')

    # Check number of objects
    if 1 != label(mask).max():
        return(0)

    # Add top/bottom slices
    shape = list(mask.shape)
    shape[0] = 1
    mask = np.vstack((np.zeros(shape), mask, np.zeros(shape)))

    # Calculate sphericity
    verts, faces, ns, vs = marching_cubes_lewiner(mask, 0.0, spacing)
    return(mesh_surface_area(verts, faces))

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

    # Fill missing dimensions
    if len(shape) > len(box):
        [box.extend([0, shape[d]]) for d in range(int(len(box)/2), len(shape))]

    # Remove extra dimensions
    box = box[:len(shape)]

    return(box)

def clear_borders(img, clean_z = None):
    '''Remove objects touching the borders of the image.
    
    Args:
      img (np.array): binary image.
      clean_z (bool): True to remove the objects touching the Z borders.
    
    Returns:
      np.array: cleaned image.
    '''

    if 2 == len(img.shape):
        img = clear_border(img)
    elif 3 == len(img.shape):
        for slide_id in range(img.shape[0]):
            img[slide_id, :, :] = clear_border(img[slide_id, :, :])
        if True == clean_z:
            for slide_id in range(img.shape[1]):
                img[:, slide_id, :] = clear_border(img[:, slide_id, :])
    return(img)

def clear_borders2(img, clean_z = None):
    '''Remove objects touching the borders of the image.
    
    Args:
      img (np.array): labeled image.
      clean_z (bool): True to remove the objects touching the Z borders.
    
    Returns:
      np.array: cleaned image.
    '''

    if 2 == len(img.shape):
        img = clear_border(img)
    elif 3 == len(img.shape):
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
        if True == clean_z:
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
    return(img)

def describe_shape(mask, spacing = None):
    """Calculate sphericity (3d) or solidity (2d) of the provided mask.
    The provided mask is expected to have only one object.

    Args:
      mask (np.array): thresholded image.
    
    Returns:
      float: shape descriptor of the provided object.
    """

    # Aspect ratio for 3d surface calculation
    if None == spacing:
        spacing = [1.0 for d in mask.shape]

    # Force binary type
    mask = mask.astype('bool')

    # Check number of objects
    if 1 != label(mask).max():
        return(0)

    # Calculate shape descriptor
    if 2 == len(mask.shape):
        # Calculate solidity
        return(float(mask.sum()) / convex_hull_image(mask).sum())
    elif 3 == len(mask.shape):
        # Add top/bottom slices
        shape = list(mask.shape)
        shape[0] = 1
        mask = np.vstack((np.zeros(shape), mask, np.zeros(shape)))

        # Calculate sphericity
        verts, faces, ns, vs = marching_cubes_lewiner(mask, 0.0, spacing)
        s = mesh_surface_area(verts, faces)
        return((np.pi * (6.0 * mask.sum())**2)**(1/3.0) / s)
    else:
        return(0)

def estimate_background(i, mask, seg_type):
    """Estimates background median.

    Args:
      i (np.array): image.
      seg_type (string): segmentation type as defined in pygpseq.const.

    Returns:
      float: estimated background.
    """

    if const.SEG_3D != seg_type and 2 != len(i.shape):
        bg = np.median(mk_z_projection(i, seg_type)[mask == 0])
    else:
        bg = np.median(i[mask == 0])

    return(bg)

def get_dtype(i):
    '''
    Identify bit depth for a matrix of maximum intensity i.
    '''
    depths = [8, 16]
    for depth in depths:
        if i <= 2**depth-1:
            return("uint%d" % (depth,))
    return("uint")

def get_mid_section_idx(i, mask, mid_type = None):
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
    if None == mid_type:
        mid_type = const.MID_SEC_DEFAULT

    # Need a stack, otherwise get the full image
    if 3 > len(i.shape):
        return(0)

    # Select central slide
    if const.MID_SEC_CENTRAL == mid_type:
        return(i.shape[0] / 2)

    if mid_type in [const.MID_SEC_LARGEST, const.MID_SEC_DEFAULT]:
        idx = [mask[idx, :, :].sum() for idx in range(mask.shape[0])]
        return(idx.index(max(idx)))

    if mid_type in [const.MID_SEC_MAXSUMI]:
        idx = [i[idx, :, :].sum() for idx in range(mask.shape[0])]
        return(idx.index(max(idx)))

def get_objects_zsize(L):
    """Retrieve objects size (2/3D).

    Args:
      L (np.array): labelled thresholded image.

    Returns:
      list: Z size of every object in the labelled image.
    """

    sizes = [(L == i).astype('int').sum(0).max() for i in range(1, L.max())]
    return(sizes)

def get_objects_xysize(L):
    """Retrieve objects size (2/3D).

    Args:
      L (np.array): labelled thresholded image.

    Returns:
      list: XY size of every object in the labelled image.
    """

    Larray = L.reshape([np.prod(L.shape)]).tolist()
    sizes = [t[1] for t in vt.uniquec(Larray) if t[0] != 0]
    return(sizes)

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
    err_log = ''

    # Copy original mask
    sel_mask = mask.copy()

    # Identify middle section
    mid = [i[sliceid, :, :].sum() for sliceid in range(mask.shape[0])]
    mid = mid.index(max(mid))

    # Select top half
    sel_mask[0:mid, :, :] = False

    # Select erosion on max projection
    d2d = distance_transform_edt(sel_mask.max(0))

    # Normalize distance
    d2d = d2d - d2d.min()
    if not 0 == d2d.max():
        d2d = d2d / float(d2d.max())
    else:
        # Log error
        err_log += 'Found 0 maximum distance at '
        err_log += str(self.s) + '.' + str(self.n) + '\n\n'

    # Make 2d mask
    mask2d = d2d >= erosion

    # Apply 2d erosion mask to every slice
    for i in range(sel_mask.shape[0]):
        sel_mask[i, :, :] = np.logical_and(sel_mask[i, :, :], mask2d)

    # Convert in unsigned integers
    sel_mask = sel_mask.astype('u4')

    # Output
    return((sel_mask, err_log))

def get_rescaling_factor(path, **kwargs):
    """Get rescaling factor for deconvolved.

    Args:
      path (string):
      **kwargs: basedir additional argument
    
    Returns:
      float: scaling factor
    """

    # Set basedir if not provided
    if not 'basedir' in kwargs.keys():
        kwargs['basedir'] = './'

    # Build proper path to the deconvolution log file
    path = kwargs['basedir'] + path[0]
    path = path.split('.')
    path[len(path) - 2] = path[len(path) - 2] + '_history'
    path[len(path) - 1] = 'txt'
    path = '.'.join(path)
    
    if not os.path.exists(path):
        # Means that the image was not deconvolved
        factor = 1
    else:
        # Identify line with scaling factor
        fhistory = open(path, 'r')
        frows = fhistory.readlines()
        fhistory.close()

        # Retrieve factor
        needle = 'Stretched to Integer type'
        factor = [x for x in frows if needle in x]
        if 0 == len(factor):
            factor = 1
        else:
            factor = factor[0]
            factor = factor.strip()
            factor = factor.split(' ')
            factor = float(factor[len(factor) - 1])

    # Output
    return(factor)

def get_unit(shape):
    """Get size unity.

    Args:
      shape (tuple[int]): image shape.

    Returns:
      string: "vx" for 3d images, "px" for 2d images.
    """
    if 2 == len(shape):
        return('px')
    elif 3 == len(shape):
        return('vx')
    else:
        return('')

def in_3d_box(box, coords):
    '''
    Check if point is in a box

    Args:
      box (tuple): ((x0, x1), (y0, y1), (z0, z1)).
      coords (tuple): (x, y, z).

    Returns
      bool
    '''
    cx = coords[0] >= box[0][0] and coords[0] <= box[0][1]
    cy = coords[1] >= box[1][0] and coords[1] <= box[1][1]
    cz = coords[2] >= box[2][0] and coords[2] <= box[2][1]
    return(cx and cy and cz)

def in_mask(coords, imbin):
    '''Check if a pixel in a mask is foreground.'''
    
    # Check the pixel is inside the image boundaries
    inbound = imbin.shape[0] > coords[0]
    inbound = inbound and imbin.shape[1] > coords[1]
    inbound = inbound and imbin.shape[2] > coords[2]
    inbound = inbound and all(np.array(coords) >= 0)
    if not inbound:
        return(False)

    # Check the pixel is foreground
    return(1 == imbin[coords[0], coords[1], coords[2]])

def mkIsoStruct(dilate_factor, aspect):
    '''
    Builds isotropic structuring element for dilation.
    
    Args:
      dilate_factor (int): number of px for isotropic 2D dilation.
      aspect (tuple(float)): voxel side ratios.
    
    Returns:
      np.ndarray: structureing element for 3D anisotropic dilation.
    '''

    # XY dilation factor
    df_xy = int(dilate_factor * 2 + 1)
    if 0 == aspect[0]:
        se = cube(df_xy)
        se = se[0]
        new_shape = [1]
        [new_shape.append(d) for d in se.shape]
        se = se.reshape(new_shape)
        return(se)
    
    # Z dilation factor
    df_z = int(dilate_factor * aspect[1] / aspect[0] * 2 + 1)

    if df_z == df_xy:
        # Isotropic
        return(cube(df_z))
    elif df_z > df_xy:
        # Larger Z side
        se = cube(df_z)
        se = se[:, 0:df_xy, 0:df_xy]
    else:
        # Larger XY side
        se = cube(df_xy)
        se = se[0:df_z]

    # Output
    return(se)

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
    return(i)

def read_tiff(impath):
    '''Read tiff image.

    Args:
      impath (string): path to tiff image.
    '''

    # Check that the file exists
    if not os.path.isfile(impath):
        msg = "Trying to read missing file"
        printout("%s:\n%s\n" % (msg, impath), -2)

    # Read TIFF (capture any parsing issues)
    try: imch = imread(impath)
    except ValueError as e:
        msg = "Something went wrong while trying to read a file"
        printout("%s (possibly corrupt):\n%s\n" % (msg, impath), -2)

    # Reshape
    imch = autoselect_time_frame(imch)

    return(imch)

def rm_from_mask(L, torm):
    # Remove elements from a mask.
    # 
    # Args:
    #     L (np.array[int]): labelled objects.
    #     torm (list): list of objects indexes (to remove).

    if len(torm) <= L.max() - len(torm):
        # Update list of objects to be discarded
        torm = [e + 1  for e in torm]

        # Identify which objects to discard
        rm_mask = np.vectorize(lambda x: x in torm)(L)

        # Discard and re-label
        L[rm_mask] = 0
    else:
        # Select objects to be kept
        tokeep = [e + 1 for e in range(L.max()) if not e in torm]

        # Identify which objects to discard
        rm_mask = np.vectorize(lambda x: not x in tokeep)(L)

        # Discard and re-label
        L[rm_mask] = 0

    # Output
    return(L > 0)

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
        return(img)

    # Slice image
    idxs = [(0,) for i in range(len(img.shape) - k)]
    for i in range(k):
        idxs.append(tuple(range(img.shape[len(img.shape) - k + i])))
    img = img[np.ix_(*idxs)].reshape(img.shape[len(img.shape) - k:])

    # Output
    return(img)

def threshold_adaptive(i, block_size):
    """Adaptive threshold.

    Args:
      i (np.array): image.
      block_size (int): neighbourhood, if even then it is incremented by 1.

    Returns:
      np.array: thresholded image.
    """

    # Increment neighbourhood size
    if 0 == block_size % 2:
        block_size += 1

    # Local threshold mask
    lmask = np.zeros(i.shape)

    # Apply threshold per slice
    if 3 == len(i.shape):
        for slice_id in range(i.shape[0]):
            lmask[slice_id, :, :] = filters.threshold_local(
                i[slice_id, :, :], block_size)
        lmask = closing(lmask, cube(3))
    else:
        lmask = filters.threshold_local(i, block_size)
        lmask = closing(lmask, square(3))

    # Output
    return(lmask)

# END ==========================================================================

################################################################################

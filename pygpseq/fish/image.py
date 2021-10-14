# -*- coding: utf-8 -*-

"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
@description: methods for FISH image manipulation.
"""

# DEPENDENCIES =================================================================

import numpy as np
import os

from skimage.measure import label
from skimage.morphology import dilation

from pygpseq.fish import dot
from pygpseq.fish import nucleus

from pygpseq.tools import Binarize
from pygpseq.tools import image as imt
from pygpseq.tools.io import IOinterface, printout
from pygpseq.tools import plot

# FUNCTIONS ====================================================================


def analyze_field_of_view(
    sid,
    data,
    im2fov,
    dilate_factor,
    istruct,
    aspect,
    mask_dir,
    mask_prefix,
    plotCompartments,
    pole_fraction,
    outdir,
    noplot,
    labeled,
    compressed,
    dist_type,
    nbins,
    discard_dilation_mode,
    an_type,
    seg_type,  # Required by the Binarize class
    mask2d_dir=None,
    verbose=False,
    debug=False,
    debug_dir="",
):
    """Given a table with FISH data, add information on:
        - lamin/center absolute/normalized distance
        - angle between homogue pairs
        - compartment assignment
        - dots coordinates in normalized ellipsoid

    Args:
        sid (int): series ID.
        data (pd.DataFrame): FISH data table.
        im2fov (dict): sid-to-absPath dictionary.
        dilate_factor (int): number of pixels for dilation
        istruct (tuple): 3D isotropic structuring element for dilation.
        aspect (tuple): ZYX voxel aspect.
        mask_dir (string): path to folder for TIFF masks import/export.
        mask_prefix (string): prefix for TIFF masks.
        plotCompartments (bool): generate compartment plots.
        pole_fraction (float):.
        outdir (str): path to analysis output folder.
        noplot (bool): turn plotting off.
        labeled (bool): import/export masks as labeled.
        compressed (bool): export masks as compressed TIFFs.
        dist_type (str): nuclear distance calculation mode.
        nbins (int): number of bins for density profile.
        an_type
        seg_type
        verbose (bool): display action log.
        debug (bool): debugging mode.
    """

    # ASSERT ===================================================================

    reqcols = ["File"]
    for c in reqcols:
        assert c in data.columns, "missing '%s' column." % c

    # INPUT ====================================================================

    v = verbose
    msg = printout("Job '%s'..." % (im2fov[sid],), 1, v)
    subt = data.loc[np.where(data["File"] == sid)[0], :]

    # Get DNA scaling factor and rescale
    sf = imt.get_rescaling_factor(im2fov[sid])
    msg += printout("Re-scaling factor: %f" % sf, 2, v)

    # Read image
    msg += printout("Reading image ...", 2, v)
    im = imt.read_tiff(im2fov[sid], k=3, rescale=sf)
    if type(None) == type(im):
        return None

    # SEGMENTATION =============================================================

    Segmenter = Binarize(an_type=an_type, seg_type=seg_type, verbose=verbose)

    # Check if already segmented
    already_segmented = False
    if type(None) != type(mask_dir):
        mpath = os.path.join(mask_dir, mask_prefix + os.path.basename(im2fov[sid]))
        already_segmented = os.path.isfile(mpath)

    # Skip or binarize
    if already_segmented:
        msg += printout("Skipped binarization, using provided mask.", 3, v)
        imbin = imt.read_tiff(mpath, k=3)

        # Check input mask for shape match
        if imbin.shape != im.shape:
            errmsg = "Inconsistent shape of provided mask and image: "
            errmsg += "%s and %s\n" % (str(im.shape), str(imbin.shape))
            msg += printout(errmsg, -2, v, False)

            msg += printout("Reverting to binarization.", 2, v)
            imbin = None
            already_segmented = False
        else:
            if not labeled:
                imbin = imbin != 0  # Binarize
            thr = 0
            already_segmented = type(None) != type(imbin)

    if not already_segmented:
        msg += printout("Binarizing...", 2, v)
        (imbin, thr, log) = Segmenter.run(im)
        msg += log

        # Filter based on object size
        imbin, tmp = Segmenter.filter_obj_XY_size(imbin)
        imbin, tmp = Segmenter.filter_obj_Z_size(imbin)

        if type(None) != type(mask2d_dir):
            mask2d_path = os.path.join(mask2d_dir, os.path.basename(im2fov[sid]))
            if os.path.isfile(mask2d_path):
                mask2d = imt.read_tiff(mask2d_path, k=2)
                # If labeled, inherit nuclei labels
                imbin = Segmenter.combine_2d_mask(imbin, mask2d, labeled2d=labeled)

    # Estimate background
    dna_bg = imt.estimate_background(im, imbin, seg_type)
    msg += printout("Estimated background: %.2f a.u." % (dna_bg,), 3, v)

    # NUCLEI ===================================================================

    msg += printout("Retrieving nuclei...", 2, v)
    L = label(imbin) if np.max(imbin) == 1 else imbin
    # Save mask ----------------------------------------------------------------

    # Export binary mask as TIF
    if type(None) != type(mask_dir) and not already_segmented:
        msg += printout("Exporting mask as tif...", 4, v)
        if not os.path.isdir(mask_dir):
            os.mkdir(mask_dir)

        # Export labeled mask
        if labeled:
            plot.save_tif(mpath, L, "uint8", compressed)
        else:  # Export binary mask (min/max)
            L[np.nonzero(L)] = 255
            plot.save_tif(mpath, L, "uint8", compressed)
            L = label(imbin)

    # Export mask as PNG
    if not noplot:
        # Create png masks output directory
        maskdir = os.path.join(outdir, "masks/")
        if not os.path.isdir(maskdir):
            os.mkdir(maskdir)
        imbname = os.path.splitext(os.path.basename(im2fov[sid]))[0]

        # Save default mask
        msg += printout("Saving default binary mask...", 3, v)
        plot.export_mask_png(
            "%smask.%s.default.png" % (maskdir, imbname), imbin, "Default mask."
        )

        # Export dilated mask
        if dilate_factor != 0:
            msg += printout("Saving dilated mask...", 3, v)
            plot.export_mask_png(
                "%smask.%s.dilated%d.png" % (maskdir, imbname, dilate_factor),
                dilation(imbin, istruct),
                "Dilated mask, %d factor." % (dilate_factor,),
            )

        # Export labeled mask
        msg += printout("Saving nuclear ID mask...", 3, v)
        plot.export_mask_png(
            "%smask.%s.nuclei.png" % (maskdir, imbname),
            L,
            'Nuclei in "%s" [%d objects]' % (os.path.basename(im2fov[sid]), L.max()),
        )

    # Store nuclei -------------------------------------------------------------
    msg += printout("Building nuclei...", 3, v)
    msg, curnuclei, dp, nv = nucleus.build_nuclei(
        msg,
        L,
        dilate_factor,
        series_id=sid,
        thr=thr,
        dna_bg=dna_bg,
        sig_bg=0,
        aspect=aspect,
        offset=(1, 1, 1),
        logpath=IOinterface().logpath,
        dist_type=dist_type,
        discard_dilation_mode=discard_dilation_mode,
        i=im,
        istruct=istruct,
        nbins=nbins,
        debug=debug,
        debug_dir=debug_dir,
    )

    # ANALYSIS =================================================================

    msg += printout("Analyzing...", 2, v)

    # Assign dots to cells -----------------------------------------------------

    msg += printout("Assigning dots to cells...", 3, v)
    subt = dot.dots2cells(subt, curnuclei, dilate_factor)

    # Distances ----------------------------------------------------------------

    msg += printout("Calculating lamina distance...", 3, v)
    subt, msg = dot.calc_dot_distances(
        msg,
        subt,
        curnuclei,
        aspect,
        dist_type,
        discard_dilation_mode=discard_dilation_mode,
    )

    # Compartments -------------------------------------------------------------

    msg += printout("Annotating compartments...", 3, v)

    # Setup condition for compartment plotting
    compdir = None
    aggdir = None
    if plotCompartments and not noplot:
        compdir = os.path.join(outdir, "compartments/")
        if not os.path.isdir(compdir):
            os.mkdir(compdir)
        aggdir = os.path.join(outdir, "agg_vis/")
        if not os.path.isdir(aggdir):
            os.mkdir(aggdir)

    # Perform annotation
    subt, tvcomp, msg = nucleus.annotate_compartments(
        msg, subt, curnuclei, compdir, pole_fraction, aspect
    )

    # Plot aggregated visualization
    nucleus.plot_nuclei_aggregated(subt, tvcomp, aspect, aggdir)

    # CONCLUDE =================================================================

    # Remove masks from curnuclei to free some memory
    for k in curnuclei.keys():
        del curnuclei[k].mask

    # Output
    msg += printout("< Finished job.", 0, v)
    return (curnuclei, subt, tvcomp, dp, nv, msg)


# END ==========================================================================

################################################################################

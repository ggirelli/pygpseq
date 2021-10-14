# -*- coding: utf-8 -*-

"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
@description: methods for nuclear channels manipulation.
"""

# DEPENDENCIES =================================================================

import matplotlib
import matplotlib.pyplot as plt

import math
import numpy as np
import os
import pandas as pd
from skimage import draw
import skimage.io as io
from skimage.morphology import dilation
import warnings

from pygpseq import const
from pygpseq.anim import Nucleus
from pygpseq.tools import image as imt, plot
from pygpseq.tools import distance as dist, stat as stt

# FUNCTIONS ====================================================================


def annotate_compartments(msg, t, nuclei, outdir, pole_fraction, aspect):
    """
    Add compartment status to dots table (by DOTTER).
    For each nucleus: the major three axes are identified, the nucleus is
    centered and rotated. Then, the dots are also centered and rotated,
    and assigned to different compartments based on the fitted ellipsoid.
    Information on the goodness of ellipsoid fit is added to the main log and
    can be extracted by grepping lines starting with "   >>>> GoF_ellipse:".

    Args:
       msg (string): log message, to be continued.
       t (pd.DataFrame): DOTTER output table.
       aspect (tuple): Z,Y,X voxel sides in real units.

    Returns:

    """

    # ASSERT ===================================================================

    reqcols = ["File", "cell_ID", "x", "y", "z"]
    for c in reqcols:
        assert c in t.columns, "missing '%s' column." % c

    # RUN ======================================================================

    # Temporarily remove dots outside cells
    nan_cond = np.isnan(t.loc[:, "cell_ID"])
    vcomp_table = pd.DataFrame()

    subt = t.loc[np.logical_not(nan_cond), :].copy()
    if 0 == subt.shape[0]:
        print(
            "!WARNING! All dots in FoV#%d are outside cells." % (t["File"].values[0],)
        )
        return (t, vcomp_table, msg)

    # Extract field
    fid = subt["File"].values[0]

    # Add missing columns
    t["compartment"] = np.nan
    t["xnorm"] = np.nan
    t["ynorm"] = np.nan
    t["znorm"] = np.nan

    # Create empty table to host compartment volume data
    vcomp_table = pd.DataFrame(index=range(1, int(subt["cell_ID"].max()) + 1))
    vcomp_table["File"] = fid
    vcomp_table["cell_ID"] = range(1, int(subt["cell_ID"].max()) + 1)
    vcomp_table["center_bot"] = np.nan
    vcomp_table["center_top"] = np.nan
    vcomp_table["poles"] = np.nan
    vcomp_table["ndots_center_bot"] = np.nan
    vcomp_table["ndots_center_top"] = np.nan
    vcomp_table["ndots_poles"] = np.nan
    vcomp_table["a"] = np.nan
    vcomp_table["b"] = np.nan
    vcomp_table["c"] = np.nan
    vcomp_table["a_slice_component"] = np.nan
    vcomp_table["a_row_component"] = np.nan
    vcomp_table["a_col_component"] = np.nan
    vcomp_table["b_slice_component"] = np.nan
    vcomp_table["b_row_component"] = np.nan
    vcomp_table["b_col_component"] = np.nan
    vcomp_table["c_slice_component"] = np.nan
    vcomp_table["c_row_component"] = np.nan
    vcomp_table["c_col_component"] = np.nan

    for cid in range(int(subt["cell_ID"].max()) + 1):
        if cid in nuclei.keys():
            msg += "    >>> Working on cell #%d...\n" % (cid,)
            cell_cond = cid == subt["cell_ID"]

            # Extract dots coordinates -----------------------------------------
            dot_coords = np.vstack(
                [
                    subt.loc[cell_cond, "z"] - nuclei[cid].box_origin[0],
                    subt.loc[cell_cond, "x"] - nuclei[cid].box_origin[1],
                    subt.loc[cell_cond, "y"] - nuclei[cid].box_origin[2],
                ]
            )

            # Center coordinates -----------------------------------------------
            x, y, z, xd, yd, zd = stt.centered_coords_3d(nuclei[cid].mask, dot_coords)
            coords = np.vstack([x, y, z])
            dot_coords = np.vstack([xd, yd, zd])

            # Rotate data ------------------------------------------------------

            # First round
            xv, yv, zv = stt.extract_3ev(coords)

            # Store axes components in compartment table
            axes_labels = ["a", "b", "c"]
            axes = [xv, yv, zv]
            for i in range(len(axes)):
                cols = [
                    "%s_row_component" % axes_labels[i],
                    "%s_col_component" % axes_labels[i],
                    "%s_slice_component" % axes_labels[i],
                ]
                vcomp_table.loc[cid, cols] = axes[i]

            # Rotate nuclei once for compartment analysis
            theta1 = stt.calc_theta(xv[0], yv[0])
            xt, yt, zt = stt.rotate3d(coords, theta1, 2)
            tcoords = np.vstack([xt, yt, zt])

            # Keep rotating for semi-axes length
            # Second round
            xv, yv, zv = stt.extract_3ev(tcoords)
            theta3 = stt.calc_theta(xv[2], zv[2])
            if np.abs(theta3) > np.pi / 2.0:
                if theta3 > 0:
                    theta3 = -np.abs(theta3 - np.pi / 2.0)
                else:
                    theta3 = np.abs(theta3 + np.pi / 2.0)
            else:
                theta3 = -np.abs(theta3 + np.pi / 2.0)
            xt, yt, zt = stt.rotate3d(tcoords, theta3, 1)
            t2coords = np.vstack([xt, yt, zt])

            # Third round
            xv, yv, zv = stt.extract_3ev(t2coords)
            theta2 = stt.calc_theta(yv[1], zv[1])
            xt, yt, zt = stt.rotate3d(t2coords, theta2, 0)
            t2coords = np.vstack([xt, yt, zt])

            # Calculate semi-axes length ---------------------------------------

            # Round up rotated coordinates
            trcoords = t2coords.astype("i")

            # Convert to rotated image
            icoords = np.transpose(trcoords) + abs(trcoords.min(1))
            trbin = np.zeros((icoords.max(0) + 1).tolist()[::-1])
            trbin[icoords[:, 2], icoords[:, 1], icoords[:, 0]] = 1

            # Calculate axes size
            zax_true_size, yax_true_size, xax_true_size = trbin.shape

            # Fit ellipsoid ----------------------------------------------------

            # Round up rotated coordinates
            trcoords = tcoords.astype("i")

            # Convert to rotated image
            icoords = np.transpose(trcoords) + abs(trcoords.min(1))
            trbin = np.zeros((icoords.max(0) + 1).tolist()[::-1])
            trbin[icoords[:, 2], icoords[:, 1], icoords[:, 0]] = 1

            # Calculate axes size
            zax_size, yax_size, xax_size = trbin.shape

            el = draw.ellipsoid(zax_size / 2.0, yax_size / 2.0, xax_size / 2.0)
            el = el[2 : (zax_size + 2), 1 : (yax_size + 1), 1 : (xax_size + 1)]

            # Calculate intersection with fitting ellipsoid
            inter_size = np.logical_and(trbin, el).sum()

            # Log intersection
            comments = []
            comments.append(
                "%s%%%s [%s.%s]."
                % (
                    round(
                        inter_size / float(trbin.sum()) * 100,
                        2,
                    ),
                    " of the nucleus is in the ellipsoid",
                    fid,
                    cid,
                )
            )
            comments.append(
                "%s%%%s [%s.%s]."
                % (
                    round(
                        inter_size / float(el.sum()) * 100,
                        2,
                    ),
                    " of the ellipsoid is in the nucleus",
                    fid,
                    cid,
                )
            )
            msg += "".join(["   >>>> GoF_ellipse: %s\n" % (s,) for s in comments])

            # Rotate dots ------------------------------------------------------

            dot_coords_t = np.vstack(stt.rotate3d(dot_coords, theta1, 2))
            # dot_coords_t = np.vstack(stt.rotate3d(dot_coords_t, theta2, 0))
            # dot_coords_t = np.vstack(stt.rotate3d(dot_coords_t, theta3, 1))

            # Assign compartments ----------------------------------------------
            # Compartment code:
            # 0 = center-top
            # 1 = center-bottom
            # 2 = pole
            c = zax_size / 2.0
            b = yax_size / 2.0
            a = xax_size / 2.0
            cf = 1 - 2 * pole_fraction
            status = np.zeros(dot_coords.shape[1])
            status[dot_coords_t[2] < 0] = 1
            status[dot_coords_t[0] > cf * a] = 2
            status[dot_coords_t[0] < -(cf * a)] = 2
            subt.loc[cell_cond, "compartment"] = status

            # Rescale coords ---------------------------------------------------

            dot_coords_t2 = dot_coords_t.copy()
            dot_coords_t2[0] = dot_coords_t[0] / a
            dot_coords_t2[1] = dot_coords_t[1] / b
            dot_coords_t2[2] = dot_coords_t[2] / zax_size * 2

            # Store rescaled coords --------------------------------------------

            subt.loc[cell_cond, "znorm"] = dot_coords_t2[2]
            subt.loc[cell_cond, "xnorm"] = dot_coords_t2[0]
            subt.loc[cell_cond, "ynorm"] = dot_coords_t2[1]

            # Calculate compartment volume -------------------------------------

            # Round up coordinates
            xt = xt.astype("i")
            zt = zt.astype("i")

            # Count voxels in compartments
            vpole = sum(xt > cf * a) + sum(xt < -(cf * a))
            centr_cond = np.logical_and(xt < (cf * a), xt > -(cf * a))
            vctop = np.logical_and(centr_cond, zt >= 0).sum()
            vcbot = np.logical_and(centr_cond, zt < 0).sum()
            vcomp_table.loc[cid, "center_top"] = vctop
            vcomp_table.loc[cid, "center_bot"] = vcbot
            vcomp_table.loc[cid, "poles"] = vpole

            # Count dots
            vcomp_table.loc[cid, "ndots_center_top"] = (status == 0).sum()
            vcomp_table.loc[cid, "ndots_center_bot"] = (status == 1).sum()
            vcomp_table.loc[cid, "ndots_poles"] = (status == 2).sum()

            # Store nucleus dimensions
            vcomp_table.loc[cid, "a"] = xax_true_size / 2.0
            vcomp_table.loc[cid, "b"] = yax_true_size / 2.0
            vcomp_table.loc[cid, "c"] = zax_true_size / 2.0

            # Assign volume information
            volume = np.zeros(dot_coords.shape[1])
            volume[:] = vctop
            volume[dot_coords_t[2] < 0] = vcbot
            volume[dot_coords_t[1] > cf * a] = vpole
            volume[dot_coords_t[1] < -(cf * a)] = vpole
            subt.loc[cell_cond, "compartment_volume"] = volume

            # Generate compartment plot with dots ------------------------------

            if not type(None) == type(outdir):
                outpng = open(
                    os.path.join(
                        outdir,
                        "%s.%s.png"
                        % (
                            fid,
                            cid,
                        ),
                    ),
                    "wb",
                )
                plt.close("all")
                plot.ortho_3d(
                    tcoords,
                    dot_coords=dot_coords_t,
                    aspect=aspect,
                    c=a * cf,
                    channels=subt.loc[cell_cond, "Channel"].values,
                )
                plt.suptitle("\n".join(comments))
                plt.savefig(outpng, format="png")
                plt.close("all")
                outpng.close()

            t.loc[np.logical_not(nan_cond), :] = subt

    return (t, vcomp_table, msg)


def build_nuclei(
    msg,
    L,
    dilate_factor,
    series_id,
    thr,
    dna_bg,
    sig_bg,
    aspect,
    offset,
    logpath,
    i,
    istruct,
    discard_dilation_mode,
    dist_type=const.LD_ARG_LABELS[const.LD_DEFAULT],
    nbins=200,
    debug=False,
    debug_dir="",
):
    """
    Build nuclei objects

    Args:
      msg (string): log message, to be continued.
      L (np.ndarray): labeled mask.
      dilate_factor (int): dilation factor.
      series_id (int): series ID.
      thr (float): global threshold value.
      dna_bg (float): DNA channel background.
      sig_bg (float): signal channel background.
      aspect (tuple): Z,Y,X voxel sides in real units.
      offset (tuple): tuple with pixel offset for bounding box.
      logpath (string): path to log file.
      i (np.array): image.
      dist_type (str): nuclear distance calculation mode.
      nbins (int): number of bins for density profile.

    Returns:
      (string, list): log message and list of Nucleus objects.
    """

    # Prepare input for Nucleus class
    kwargs = {
        "series_id": series_id,
        "thr": thr,
        "dna_bg": dna_bg,
        "sig_bg": sig_bg,
        "aspect": aspect,
        "offset": offset,
        "logpath": logpath,
        "i": i,
    }

    # Default nuclear ID list and empty dictionary
    curnuclei = {}
    dp = []
    nv = []

    # Log operation
    msg += "   - Saving %d nuclei" % L.max()
    if 0 != dilate_factor:
        msg += " with dilation [%d]" % dilate_factor
    msg += "...\n"

    # Iterate through nuclei
    for n in range(1, L.max() + 1):
        if 0 == np.sum((L == n).astype(int)):
            continue

        original_mask = L == n

        # Make nucleus
        if 0 != dilate_factor:
            # With dilated mask
            mask = dilation(L == n, istruct)
        else:
            mask = L == n
        nucleus = Nucleus(n=n, mask=mask, **kwargs)

        # Apply box
        msg += "    > Applying nuclear box [%d]...\n" % (n,)
        mask = imt.apply_box(mask, nucleus.box)
        original_mask = imt.apply_box(original_mask, nucleus.box)

        # Store nucleus
        nucleus.mask = mask
        nucleus.original_mask = original_mask
        nucleus.dilate_factor = dilate_factor
        curnuclei[n] = nucleus

        # Density profile ------------------------------------------------------
        if discard_dilation_mode:
            mask = original_mask

        laminD, centrD = dist.calc_nuclear_distances(dist_type, mask, aspect)
        laminD_norm = dist.normalize_nuclear_distance(dist_type, laminD, centrD)

        if debug:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                ipath = "series%d.nucleus%d.tif" % (series_id, n)
                io.imsave(os.path.join(debug_dir, "mask.%s" % ipath), mask.astype("u4"))
                io.imsave(
                    os.path.join(debug_dir, "laminD.%s" % ipath),
                    laminD.astype(np.uint32),
                )
                io.imsave(
                    os.path.join(debug_dir, "laminD_norm.%s" % ipath),
                    (laminD_norm * 255).astype(np.uint32),
                )
                io.imsave(
                    os.path.join(debug_dir, "centrD.%s" % ipath),
                    centrD.astype(np.uint32),
                )
        laminD_norm = laminD_norm[mask].flatten()

        dna = imt.apply_box(i, nucleus.box)[mask].flatten()
        dp.append(nucleus.calc_density_profile(dna, laminD_norm, nbins))
        hist, edges = np.histogram(laminD_norm, nbins)
        out = [nucleus.c, nucleus.s, nucleus.n]
        out.extend(hist.tolist())
        nv.append(out)

    col_labs = ["c", "s", "n"]
    col_labs.extend(["nd_%f" % b for b in np.linspace(0, 1, nbins + 1)[1:]])
    if 0 != len(dp):
        dp = pd.DataFrame(np.vstack(dp))
        dp.columns = col_labs
    else:
        dp = pd.DataFrame()

    if 0 != len(nv):
        nv = pd.DataFrame(np.vstack(nv))
        nv.columns = col_labs
    else:
        nv = pd.DataFrame()

    return (msg, curnuclei, dp, nv)


def flag_G1_cells(t, nuclei, outdir, dilate_factor, dot_file_name):
    """
    Assign a binary flag identifying the predominant cell population
    based on flatten size and intensity sum

    Args:
      t (pd.DataFrame): DOTTER output table.
      nuclei (list(gp.Nucleus)): identified nuclei.
      outdir (string): path to output folder.
      dilate_factor (int): number of dilation operations.
      dot_file_name (string): output file name.

    Returns:
      pd.DataFrame:.
    """

    print("> Flagging G1 cells...")

    # Retrieve nuclei summaries ------------------------------------------------
    print("   > Retrieving nuclear summary...")
    summary = np.zeros(len(nuclei), dtype=const.DTYPE_NUCLEAR_SUMMARY_3D)
    for i in range(len(nuclei)):
        summary[i] = nuclei[i].get_summary()

    # Filter nuclei ------------------------------------------------------------
    print("   > Filtering nuclei based on flatten size and intensity...")
    cond_name = "none"
    sigma = 0.1
    nsf = (const.NSEL_FLAT_SIZE, const.NSEL_SUMI)
    out_dir = "."

    # Filter features
    sel_data = {}
    ranges = {}
    plot_counter = 1
    for nsfi in nsf:
        # Identify Nuclear Selection Feature
        nsf_field = const.NSEL_FIELDS[nsfi]
        nsf_name = const.NSEL_NAMES[nsfi]
        print("   >> Filtering %s..." % (nsf_name,))

        # Start building output
        d = {"data": summary[nsf_field]}

        # Calculate density
        d["density"] = stt.calc_density(d["data"], sigma=sigma)

        # Identify range
        args = [d["density"]["x"], d["density"]["y"]]
        d["fwhm_range"] = stt.get_fwhm(*args)
        ranges[nsf_name] = d["fwhm_range"]

        # Plot
        sel_data[nsf_field] = d

    # Select based on range
    f = lambda x, r: x >= r[0] and x <= r[1]
    for nsfi in nsf:
        nsf_field = const.NSEL_FIELDS[nsfi]
        nsf_name = const.NSEL_NAMES[nsfi]
        print("   > Selecting range for %s ..." % (nsf_name,))

        # Identify nuclei in the FWHM range
        nsf_data = sel_data[nsf_field]
        nsf_data["sel"] = [f(i, nsf_data["fwhm_range"]) for i in nsf_data["data"]]
        sel_data[nsf_field] = nsf_data

    # Select those in every FWHM range
    print("   > Applying selection criteria")
    nsfields = [const.NSEL_FIELDS[nsfi] for nsfi in nsf]
    selected = [sel_data[f]["sel"] for f in nsfields]
    g = lambda i: all([sel[i] for sel in selected])
    selected = [i for i in range(len(selected[0])) if g(i)]
    sub_data = np.array(summary[selected])

    # Identify selected nuclei objects
    sel_nuclei_labs = ["_%d.%d_" % (n, s) for (n, s) in sub_data[["s", "n"]]]
    sel_nucl = [n for n in nuclei if "_%d.%d_" % (n.s, n.n) in sel_nuclei_labs]

    # Check which dots are in which nucleus and update flag --------------------
    print("   > Matching DOTTER cells with GPSeq cells...")
    t["G1"] = 0
    t.loc[np.where(np.isnan(t["cell_ID"]))[0], "G1"] = np.nan
    t["universalID"] = [
        "_%s.%s_" % x for x in zip(t["File"].values, t["cell_ID"].values)
    ]
    g1ids = [i for i in range(t.shape[0]) if t.loc[i, "universalID"] in sel_nuclei_labs]
    t.loc[g1ids, "G1"] = 1
    t = t.drop("universalID", 1)

    # Add G1 status to summary -------------------------------------------------
    summary = pd.DataFrame(summary)
    summary["G1"] = np.zeros((summary.shape[0],))
    summary["universalID"] = [
        "_%s.%s_" % x for x in zip(summary["s"].values, summary["n"].astype("f").values)
    ]
    g1ids = [
        i
        for i in range(summary.shape[0])
        if summary.loc[i, "universalID"] in sel_nuclei_labs
    ]
    summary.loc[g1ids, "G1"] = 1
    summary = summary.drop("universalID", 1)

    # Estimate radius ----------------------------------------------------------
    summary["sphere_radius"] = summary["size"].values * 3 / (4 * math.pi)
    summary["sphere_radius"] = (summary["sphere_radius"]) ** (1 / 3.0)

    # Export -------------------------------------------------------------------

    # Export feature ranges
    s = "".join(["%s\t%f\t%f\n" % (k, v[0], v[1]) for (k, v) in ranges.items()])
    with open("%s/feature_ranges.txt" % (outdir,), "w+") as f:
        f.write(s)

    # Export summary
    outname = "%s/nuclei.out.dilate%d.%s" % (outdir, dilate_factor, dot_file_name)
    summary.to_csv(outname, sep="\t", index=False)

    # Output -------------------------------------------------------------------
    print("> Flagged G1 cells...")
    return t


def plot_nuclei_aggregated(t, nt, aspect, outdir=None):
    """Generate aggregated visualization for the provided data.

    Args:
            t (pd.DataFrame): FISH data frame.
            nt (pd.DataFrame): nuclear compartment data frame.
            outdir (str): output directory.
    """

    # Don't print if no output folder is provided
    if type(None) == type(outdir):
        return
    if 0 == t.shape[0] or 0 == nt.shape[0]:
        return

    # Asserts ------------------------------------------------------------------

    assert os.path.isdir(outdir), "output folder does not exist."
    for c in ["File", "Channel", "xnorm", "ynorm", "znorm"]:
        assert c in t.columns, "missing '%s' column." % c
    for c in ["a", "b", "c"]:
        assert c in nt.columns, "missing '%s' column." % c

    if not np.logical_not(np.isnan(t["xnorm"].values)).any():
        return

    for colname in ["xnorm", "ynorm", "znorm"]:
        t = t.loc[np.isfinite(t[colname].values), :]
        if 0 == t.shape[0]:
            return

    for colname in ["a", "b", "c"]:
        nt = nt.loc[np.isfinite(nt[colname].values), :]
        if 0 == nt.shape[0]:
            return

    # Make aggregated visualization --------------------------------------------

    # Calculate median a/b/c
    a = np.nanmedian(nt["a"].values)
    b = np.nanmedian(nt["b"].values)
    c = np.nanmedian(nt["c"].values)

    coords = np.vstack(
        [t["xnorm"].values * a, t["ynorm"].values * b, t["znorm"].values * c]
    )

    if 1 == len(set(t["File"].values)):
        fid = "f_%d" % t["File"].values[0]
    else:
        fid = "f_all"

    # Plot ---------------------------------------------------------------------

    # All-channels visualization
    plot.dots_in_ellipsoid(
        a,
        b,
        c,
        coords,
        aspect=aspect,
        channels=t["Channel"].values,
        title="[%s] Aggregated FISH visualization" % fid,
        outpng=os.path.join(outdir, "%s.aggregated.all.png" % fid),
    )

    # All-channels folded visualization
    plot.dots_in_ellipsoid(
        a,
        b,
        c,
        coords,
        aspect=aspect,
        channels=t["Channel"].values,
        fold=True,
        title="[%s] Aggregated-folded FISH visualization" % fid,
        outpng=os.path.join(outdir, "%s.aggregated.all.folded.png" % fid),
    )

    # Single channel visualization
    for channel in set(t["Channel"].values.tolist()):
        title = "[%s] Aggregated FISH visualization" % fid
        title += " for channel '%s'" % channel
        plot.dots_in_ellipsoid(
            a,
            b,
            c,
            coords[:, t["Channel"].values == channel],
            aspect=aspect,
            channels=t["Channel"].values[t["Channel"].values == channel],
            title=title,
            outpng=os.path.join(outdir, "%s.aggregated.%s.png" % (fid, channel)),
        )

        title = "[%s] Aggregated-folded FISH visualization" % fid
        title += " for channel '%s'" % channel
        plot.dots_in_ellipsoid(
            a,
            b,
            c,
            coords[:, t["Channel"].values == channel],
            channels=t["Channel"].values[t["Channel"].values == channel],
            aspect=aspect,
            fold=True,
            title=title,
            outpng=os.path.join(outdir, "%s.aggregated.%s.folded.png" % (fid, channel)),
        )


# END ==========================================================================

################################################################################

# -*- coding: utf-8 -*-

"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
@description: methods for FISH dots manipulation.
"""

# DEPENDENCIES =================================================================

import numpy as np
import pandas as pd

from scipy.interpolate import RegularGridInterpolator
from scipy.ndimage.morphology import distance_transform_edt

from pygpseq.tools import distance as dist, image as imt, stat as stt

# FUNCTIONS ====================================================================


def add_allele(data):
    """
    Add allele labels to DOTTER-based table with GPSeq-like centrality.

    Labels:
      NaN : dot outside of cells.
      -1  : more than 2 dots per cell.
      0   : less than 2 dots per cell.
      1   : central dot.
      2   : peripheral dot.

    Args:
      data (pd.DataFrame): DOTTER-based table with GPSeq-like centrality.
                                               Required columns:
                                                      cell_ID, lamin_dist_norm, File, Channel
            Returns:
      pd.DataFrame: input data table with added Allele column (label).
    """

    # Initial checks -----------------------------------------------------------

    # Check that the format corresponds
    assert_msg = "input should be a DataFrame from the pandas library."
    assert type(data) == type(pd.DataFrame()), assert_msg

    # Check that required columns are present
    req_cols = ["cell_ID", "lamin_dist_norm", "File", "Channel"]
    check_cols = [True for c in req_cols if c in data.columns.tolist()]
    miss_cols = [req_cols[i] for i in range(len(req_cols)) if not check_cols[i]]
    assert 0 == len(miss_cols), "Some required columns are missing: %s" % (
        ", ".join(miss_cols),
    )

    # Universal index and dots in cells ----------------------------------------

    # Default value of np.nan for dots outside of nuclei
    data["Allele"] = np.nan

    # Identify dots within cells
    validIdx = np.where(np.logical_not(np.isnan(data["cell_ID"])))[0]
    subt = data.loc[validIdx, :]

    # Assemble universal index
    subt["universalID"] = [
        "%s_%s_%s" % t
        for t in zip(
            subt["File"].values, subt["Channel"].values, subt["cell_ID"].values
        )
    ]

    # Count dots per universalID
    uID, uCount = np.unique(
        subt.loc[validIdx, "universalID"], return_index=False, return_counts=True
    )
    IDmap = zip(
        subt.loc[validIdx, "universalID"],
        [dict(zip(uID, uCount))[ID] for ID in subt.loc[validIdx, "universalID"]],
    )
    IDmap = np.array(list(IDmap))

    # Stop if no dots are inside a cell
    if sum(IDmap.shape) == 0:
        return data

    # Fill Allele column -------------------------------------------------------

    # -1 if more than 2 dots
    cond = IDmap[:, 1].astype("i") > 2
    if 0 != sum(cond):
        subt.loc[validIdx[cond], "Allele"] = -1

    #  0 if less than 2 dots
    cond = IDmap[:, 1].astype("i") == 1
    if sum(cond) != 0:
        subt.loc[validIdx[cond], "Allele"] = 0

    # Iterate over 2-dots cases
    cond = IDmap[:, 1].astype("i") == 2
    if sum(cond) != 0:
        uID = np.unique(IDmap[cond, 0]).tolist()
        for ID in uID:
            dotPair = subt.loc[subt["universalID"] == ID, :]
            ldn = dotPair["lamin_dist_norm"].tolist()
            if ldn[0] == ldn[1]:  # Same centrality
                # Central
                subt.loc[dotPair.index[0], "Allele"] = 1
                # Peripheral
                subt.loc[dotPair.index[1], "Allele"] = 2
            else:  # Different centrality
                # Peripheral
                subt.loc[dotPair["lamin_dist_norm"].idxmin(), "Allele"] = 2
                # Central
                subt.loc[dotPair["lamin_dist_norm"].idxmax(), "Allele"] = 1

    # Output -------------------------------------------------------------------
    data.loc[validIdx, "Allele"] = subt["Allele"]
    return data


def add_allele_polarity(t, nuclei, aspect):
    """Add inter-homologous angle in respect to the nucleus center of mass
    for homologous couples.

    Args:
            t (pd.DataFrame): FISH data frame.
            nuclei (list): list of pygpseq.anim.Nucleus instances.
            aspect (tuple): Z,Y,X voxel sides in real units.

    Returns:
            pd.DataFrame: updated FISH data frame.
    """

    # Check input
    reqcols = ["File", "Channel", "cell_ID", "Allele"]
    for c in reqcols:
        assert c in t.columns, "missing '%s' column." % c
    az, ay, ax = aspect

    # Assemble universal index
    t.loc[:, "universalID"] = [
        "%s_%s_%s" % x
        for x in zip(t["File"].values, t["Channel"].values, t["cell_ID"].values)
    ]

    # Subset data to Allele columns
    subt = t.loc[t["Allele"] > 0, :]

    # Set default value for angle and com columns
    t["angle"] = np.nan
    # t['com'] = np.nan

    # Go through cells ---------------------------------------------------------
    for uid in subt["universalID"]:
        idx = subt[subt["universalID"] == uid].index

        # Retrieve allele coordinates
        focus = subt.loc[subt["universalID"] == uid, ("x", "y", "z")]
        if sum(focus.shape) == 0:
            continue

        # Identify nucleus
        cid = subt.loc[subt["universalID"] == uid, "cell_ID"].values[0]
        sid = subt.loc[subt["universalID"] == uid, "File"].values[0]
        if np.isnan(cid) or np.isnan(sid):
            continue

        # Retrieve nucleus
        ncond = [n.s == sid and n.n == cid for n in nuclei]
        if not any(ncond):
            print(
                "Nucleus not found for %s.%s"
                % (
                    sid,
                    cid,
                )
            )
            continue
        else:
            nucleus = [nuclei[i] for i in range(len(ncond)) if ncond[i]][0]

        # Calculate nucleus center of mass coordinates
        C = (nucleus.box_mass_center + nucleus.box_origin).astype("i")
        C = C[[1, 2, 0]]
        # t.loc[idx, 'com'] = "_".join([str(x) for x in C.tolist()])

        # Calculate angle
        P1 = focus.loc[focus.index[0], :]
        P2 = focus.loc[focus.index[1], :]

        if all(P2 == C) or all(P1 == C):
            t.loc[idx, "angle"] = 0
        else:
            # Calculate angle
            xyz_aspect = np.array((ax, ay, az))
            t.loc[idx, "angle"] = stt.angle_between_points(
                P1 * xyz_aspect, C * xyz_aspect, P2 * xyz_aspect
            )

    # Remove universal ID
    t = t.drop("universalID", 1)

    return t


def calc_dot_distances(msg, t, nuclei, aspect, dist_type, discard_dilation_mode=False):
    """
    Calculate distance of dots from lamina and central area

    Args:
      msg (string): log message, to be continued.
      t (pd.DataFrame): DOTTER output table.
      nuclei (list(gp.Nucleus)): identified nuclei.
      aspect (tuple): Z,Y,X voxel sides in real units.
      dist_type (str): nuclear distance calculation mode.

    Returns:
      pd.DataFrame: updated dotter table.
      str: message log..
    """

    t["lamin_dist"] = np.nan
    t["lamin_dist_norm"] = np.nan
    t["centr_dist"] = np.nan
    t["centr_dist_norm"] = np.nan

    # Skip if no cells are present
    if np.all(np.isnan(t["cell_ID"].values)):
        return (t, msg)

    # Calculate distances ------------------------------------------------------
    max_cell_ID = int(np.nanmax(t["cell_ID"].values))
    for cid in (i for i in range(max_cell_ID + 1) if i in nuclei.keys()):
        msg += "    >>> Working on cell #%d...\n" % (cid,)

        cell_cond = cid == t["cell_ID"]
        if sum(cell_cond) == 0:
            continue

        mask = nuclei[cid].original_mask if discard_dilation_mode else nuclei[cid].mask
        # Perform EDT and normalize
        laminD, centrD = dist.calc_nuclear_distances(dist_type, mask, aspect)
        laminD_norm = dist.normalize_nuclear_distance(dist_type, laminD, centrD)

        # Interpolate EDT maps
        gridShape = np.array(laminD.shape)
        regularGrid = [
            np.linspace(0, gridShape[i] - 1, gridShape[i])
            for i in range(len(gridShape))
        ]
        laminD_gradient = RegularGridInterpolator(
            regularGrid, laminD, method="linear", bounds_error=False, fill_value=0
        )
        centrD_gradient = RegularGridInterpolator(
            regularGrid, centrD, method="linear", bounds_error=False, fill_value=0
        )
        laminD_norm_gradient = RegularGridInterpolator(
            regularGrid, laminD_norm, method="linear", bounds_error=False, fill_value=0
        )

        # Store distances
        st = t.loc[cell_cond, :].copy()
        for i in t.loc[cell_cond, :].index:
            absolute_coords = t.loc[i, ["z", "x", "y"]].tolist()
            nuclear_box_origin = nuclei[cid].box_origin
            relative_coords = absolute_coords - nuclear_box_origin
            st.loc[i, "lamin_dist"] = laminD_gradient(relative_coords)
            st.loc[i, "centr_dist"] = centrD_gradient(relative_coords)
            st.loc[i, "lamin_dist_norm"] = laminD_norm_gradient(relative_coords)
        cols = ["lamin_dist", "centr_dist", "lamin_dist_norm"]
        t.loc[cell_cond, cols] = st.loc[:, cols].values

        lamin_dist_norm_values = t.loc[cell_cond, "lamin_dist_norm"].values
        t.loc[cell_cond, "centr_dist_norm"] = np.absolute(1 - lamin_dist_norm_values)

    # Output
    return (t, msg)


def dots2cells(t, nuclei, dilate_factor):
    """
    Assign dots to cells

    Args:
      t (pd.DataFrame): DOTTER output subset.
      nuclei (list(gp.Nucleus)): identified nuclei.
      dilate_factor (int): number of dilation operations.

    Returns:
      pd.DataFrame: updated DOTTER output.
    """

    t["cell_ID"] = np.nan
    for idx in t.index:
        coords = (t.loc[idx, "zi"], t.loc[idx, "xi"], t.loc[idx, "yi"])
        for (nid, n) in nuclei.items():
            if imt.in_mask(coords - n.box_origin, n.mask):
                t.loc[idx, "cell_ID"] = nid
                continue

    # Output
    return t


# END ==========================================================================

################################################################################

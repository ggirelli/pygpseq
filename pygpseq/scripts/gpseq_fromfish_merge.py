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
# Date: 2018
# Project: GPSeq/FISH
# Description: script to merge multiple gpseq_fromfish outputs.
#
# ------------------------------------------------------------------------------


# DEPENDENCIES =================================================================

import argparse
import datetime
from joblib import Parallel, delayed
import multiprocessing
import numpy as np
import os
import pandas as pd
import pygpseq as gp
import sys
import time

from pygpseq.tools.io import printout

# FUNCTION =====================================================================


def merge_DataFrameDict_byKey(dfd, k):
    out = [d[k] for d in dfd if type(d[k]) == type(pd.DataFrame())]
    out = pd.DataFrame() if not out else pd.concat(out)
    out.index = range(out.shape[0])
    return out


def look_for_data(d, flist, k, flag, dataset, ipath):
    """Look for tables in input folder.
    If the table was not found, store np.nan.

    Args:
            d (dict): dictionary to update.
            flist (list): possible file list.
            k (str): dictionary key for storage.
            flag (str): dataset_series flag.
            dataset (str): dataset label.
            ipath (str): input folder path.

    Return
            dict: updated dictionary.
    """
    if len(flist) == 0:
        d["error"] = "Warning! Cannot find %s information in %s. %s" % (
            k,
            flag,
            "Dataset %s. Folder: %s" % (dataset, ipath),
        )
        d["good"] = False
        d["partial"] = True
        d[k] = np.nan
    else:
        d[k] = pd.read_csv("%s/%s" % (ipath, flist[0]), args.delim)
    return d


def add_dataset_info(data, did, date, sid, ct):
    data["dataset"] = did
    data["date"] = date
    data["session"] = sid
    data["cell_type"] = ct
    return data


def extract_data(did, date, sid):
    # Subset metadata
    time.sleep(1)  # To avoid automatic batch resizing
    log = []

    sub_cond = np.logical_and.reduce(
        (
            md["dataset"].values == did,
            md["date"].values == date,
            md["session"].values == sid,
        )
    )
    subt = md.loc[
        sub_cond,
    ]
    subt.index = range(subt.shape[0])

    # Extract dataset and series information -----------------------------------
    dataset = did
    date = date
    session = "%03d" % sid
    flag = "%s_%s_%s" % (dataset, date, session)
    cell_type = subt.loc[0, "cell_line"]

    # Log current status
    log.append("Working on %s..." % flag)

    # Identify input folder and skip if missing --------------------------------
    outl = []
    for indir in args.indir:
        ipath = "%s/%s" % (indir, flag)

        # Prepare flag output
        out = {"good": True, "partial": False}
        if not os.path.isdir(ipath):
            msg = "Warning! Cannot find folder for %s." % flag
            msg += "\nSkipped. Folder: %s" % ipath
            log.append(msg)
            out["good"] = False
            continue

        # Identify input files -------------------------------------------------
        flist = os.listdir(ipath)
        nuclei = [x for x in flist if "nuclei.out" in x]
        dots = [x for x in flist if "wCentr.out" in x and "noAllele" not in x]
        comps = [x for x in flist if "nuclear_compartment." in x]
        dens = [x for x in flist if "density_profile.tsv" in x]
        vols = [x for x in flist if "volume_profile.tsv" in x]

        out = look_for_data(out, nuclei, "nuclei", flag, dataset, ipath)
        out = look_for_data(out, dots, "dots", flag, dataset, ipath)
        out = look_for_data(out, comps, "compartments", flag, dataset, ipath)
        out = look_for_data(out, dens, "density_profile", flag, dataset, ipath)
        out = look_for_data(out, vols, "volume_profile", flag, dataset, ipath)

        outl.append(out)

    # Check for completeness
    if any(d["good"] for d in outl):
        d = [d for d in outl if d["good"]][0]

    elif any(d["partial"] for d in outl):
        # Print error due to partial information present
        [log.append(d["error"]) for d in outl if d["partial"]]
        d = [d for d in outl if d["partial"]][0]
    else:
        # Print error due to not-found information
        msg = "Warning! Cannot find information on dataset %s %s" % (
            dataset,
            "in any of the input directories.",
        )
        msg += "\nSkipped %s.\n" % flag
        print(msg)
        return np.nan
    # Extract dot data ---------------------------------------------------------
    if type(pd.DataFrame()) == type(d["dots"]):
        dots = d["dots"]
        dots = add_dataset_info(dots, did, date, sid, cell_type)
        dch = {
            subt.loc[i, "channel"].lower(): subt.loc[i, "probe_label"]
            for i in subt.index
        }

        dots["probe_label"] = np.nan
        for i in dots.index:
            c = dots.loc[i, "Channel"].lower()
            if c in dch:
                dots.loc[i, "probe_label"] = dch[c]
    else:
        dots = np.nan

    # Extract nuclei data ------------------------------------------------------
    if type(pd.DataFrame()) == type(d["nuclei"]):
        nuclei = d["nuclei"]
        nuclei = add_dataset_info(nuclei, did, date, sid, cell_type)
    else:
        nuclei = np.nan

    # Extract compartments data ------------------------------------------------
    if type(pd.DataFrame()) == type(d["compartments"]):
        comps = d["compartments"]
        comps = add_dataset_info(comps, did, date, sid, cell_type)
    else:
        comps = np.nan

    # Extract density profile --------------------------------------------------
    if type(pd.DataFrame()) == type(d["density_profile"]):
        dens = d["density_profile"]
        dens = add_dataset_info(dens, did, date, sid, cell_type)
    else:
        dens = np.nan

    # Extract volume profile ---------------------------------------------------
    if type(pd.DataFrame()) == type(d["volume_profile"]):
        vols = d["volume_profile"]
        vols = add_dataset_info(vols, did, date, sid, cell_type)
    else:
        vols = np.nan

    # Prepare allele by channel table ------------------------------------------
    if type(pd.DataFrame()) == type(dots):
        aldata = dots.loc[
            np.logical_not(np.isnan(dots["Allele"].values)),
        ]
        aldata = aldata.loc[
            aldata["Allele"].values > 0,
        ]

        if aldata.shape[0] != 0:
            al_uniID = set(zip(aldata["File"], aldata["Channel"], aldata["cell_ID"]))

            dl = []
            cols = []

            for (fid, chid, cid) in al_uniID:
                al_cond = np.logical_and(
                    np.logical_and(
                        aldata["File"].values == fid, aldata["Channel"].values == chid
                    ),
                    aldata["cell_ID"].values == cid,
                )
                alt = aldata.loc[
                    al_cond,
                ]
                alt.index = range(alt.shape[0])

                d = alt.loc[:, ["File", "Channel", "cell_ID", "G1"]]

                d_3d = alt.loc[0, ["x", "y", "z"]].values
                d_3d -= alt.loc[1, ["x", "y", "z"]].values

                d["d_3d"] = np.sqrt(np.sum(np.power(d_3d * args.aspect, 2)))
                d["d_lamin"] = np.abs(np.diff(alt["lamin_dist"].values))[0]
                d["d_lamin_norm"] = np.abs(np.diff(alt["lamin_dist_norm"].values))[0]
                d["d_centr"] = np.abs(np.diff(alt["centr_dist"].values))[0]
                d["d_centr_norm"] = np.abs(np.diff(alt["centr_dist_norm"].values))[0]
                d["angle"] = alt["angle"].values[0]
                d = add_dataset_info(d, did, date, sid, cell_type)
                d["probe_label"] = alt["probe_label"].values[0]

                dl.append(d.loc[0, :])
                cols = d.columns

            alleles = pd.concat(dl, 1).transpose()
            alleles.index = range(alleles.shape[0])
            alleles.columns = cols
        else:
            msg = "Warning! No homologue copy pairs found in %s" % flag
            log.append("%s, homologue copy calculation skipped.\n" % msg)
            alleles = np.nan
    else:
        log.append("Warning! No dots found in %s.\n" % flag)
        dots = np.nan
        alleles = np.nan

    log.append("Finished %s" % flag)
    return {
        "dots": dots,
        "nuclei": nuclei,
        "comps": comps,
        "dens": dens,
        "vols": vols,
        "alleles": alleles,
    }


def run():

    # PARAMETERS ===================================================================

    # Add script description
    parser = argparse.ArgumentParser(
        description="""
    Merge gpseq_fromfish outputs, add dataset and cell type information. The script
    looks for the output in subfolders of the specified input directories. These
    subfolders must be named as dataset_date_session, with session being in XXX
    format with leading zeros (e.g., iJC001_YYYYMMDD_001, iJC001_YYYYMMDD_002,...).
    Please, note that the channel is enforced as lower-case by the merge operation.

    The script will not only merge the input, but also produce a table of
    informative homologue copy pairs features, if any pair is found in the input.

    Example 1: output in current directory.
    gpseq_fromfish_merge -m meta.tsv -i HAP1/dots_auto IMR90/dots_auto

    Example 2: output to "/home/user/out" directory.
    gpseq_fromfish_merge -m meta.tsv -i HAP1/dots_auto IMR90/dots_auto
        -o /home/user/out
    """,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # Positional parameters
    parser.add_argument(
        "-m",
        "--meta",
        metavar="metadata",
        required=True,
        type=str,
        help="""Path to metadata table. Required columns:
        dataset, date, session, cell_line, probe_label, channel.""",
    )
    parser.add_argument(
        "-i",
        "--indir",
        metavar="indir",
        required=True,
        type=str,
        nargs="+",
        help="""List of space-separated input folder
        paths with gpseq_fromfish output as subfolders.""",
    )

    # Optional parameters
    parser.add_argument(
        "-o",
        "--outdir",
        metavar="outdir",
        type=str,
        help="""Path to output folder, created if missing. Defaults to current
        position.""",
        default=".",
    )
    parser.add_argument(
        "-a",
        "--aspect",
        type=float,
        nargs=3,
        help="""Physical size of Z, Y and X voxel sides.
        Default: 300.0 130.0 130.0""",
        metavar=("Z", "Y", "X"),
        default=[300.0, 130.0, 130.0],
    )
    parser.add_argument(
        "-d",
        "--delim",
        metavar="sep",
        type=str,
        help="""Column separator in input tables. Default: TAB""",
        default="\t",
    )
    parser.add_argument(
        "-t",
        "--threads",
        metavar="threads",
        type=int,
        default=1,
        help="""Number of threads to be used for parallelization.
        Increasing the number of threads might increase the required amount of RAM.
        """,
    )
    parser.add_argument(
        "--no-date",
        action="store_const",
        dest="addDate",
        const=False,
        default=True,
        help="Do not add date as prefix to output.",
    )

    # Version flag
    version = "4.0.1"
    parser.add_argument(
        "--version",
        action="version",
        version="%s v%s"
        % (
            sys.argv[0],
            version,
        ),
    )

    # Parse arguments
    args = parser.parse_args()

    # Additional checks
    if type([]) != type(args.indir):
        args.indir = [args.indir]
    assert os.path.isfile(args.meta), "metadata file not found: %s" % args.meta
    for idp in args.indir:
        assert os.path.isdir(idp), "input folder not found: %s" % idp
    maxncores = multiprocessing.cpu_count()
    if maxncores < args.threads:
        printout(
            "Lowered number of threads to maximum available: %d\n" % (maxncores), -1
        )
        args.threads = maxncores
    if 0 >= args.threads:
        args.threads = 1

    # RUN ==========================================================================

    # Read metadata table
    md = pd.read_csv(args.meta, args.delim)

    # Check that all mandatory metadata columns
    req_col = ["dataset", "date", "session", "cell_line", "probe_label", "channel"]
    for c in req_col:
        assert_msg = "%s '%s' in metadata table.\nMeta: %s\n" % (
            "Missing mandatory column",
            c,
            args.meta,
        )
        assert c in md.columns, assert_msg

    # Iterate through Fields of View (FoVs)
    uniID = set(zip(md["dataset"].values, md["date"].values, md["session"].values))

    print("Start parallel")
    l2 = Parallel(n_jobs=args.threads, verbose=11)(
        delayed(extract_data)(did, date, sid) for (did, date, sid) in uniID
    )
    print("End parallel")

    # Remove skipped
    l2 = [x for x in l2 if type(x) == type({})]

    # Merge dataset outputs --------------------------------------------------------

    print("Merging copy information...")
    alleles = merge_DataFrameDict_byKey(l2, "alleles")

    print("Merging dot information...")
    dots = merge_DataFrameDict_byKey(l2, "dots")

    print("Merging compartment information...")
    comps = merge_DataFrameDict_byKey(l2, "comps")

    print("Merging nuclear information...")
    nuclei = merge_DataFrameDict_byKey(l2, "nuclei")

    print("Merging density profile information...")
    dens = merge_DataFrameDict_byKey(l2, "dens")

    print("Merging volume profile information...")
    vols = merge_DataFrameDict_byKey(l2, "vols")

    # Write output
    if not os.path.isdir(args.outdir):
        os.mkdir(args.outdir)
    curdate = datetime.datetime.now().isoformat().split("T")[0] + "_"
    if not args.addDate:
        curdate = ""
    dots.to_csv(
        "%s/%s%s" % (args.outdir, curdate, "dots.merged.tsv"),
        args.delim,
        index=False,
        na_rep="NA",
    )
    alleles.to_csv(
        "%s/%s%s" % (args.outdir, curdate, "copies.merged.tsv"),
        args.delim,
        index=False,
        na_rep="NA",
    )
    nuclei.to_csv(
        "%s/%s%s" % (args.outdir, curdate, "nuclei.merged.tsv"),
        args.delim,
        index=False,
        na_rep="NA",
    )
    dens.to_csv(
        "%s/%s%s" % (args.outdir, curdate, "density_profile.merged.tsv"),
        args.delim,
        index=False,
        na_rep="NA",
    )
    vols.to_csv(
        "%s/%s%s" % (args.outdir, curdate, "volume_profile.merged.tsv"),
        args.delim,
        index=False,
        na_rep="NA",
    )
    comps.to_csv(
        "%s/%s%s" % (args.outdir, curdate, "ncomps.merged.tsv"),
        args.delim,
        index=False,
        na_rep="NA",
    )

    # End --------------------------------------------------------------------------

    ################################################################################

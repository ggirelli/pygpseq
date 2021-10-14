# -*- coding: utf-8 -*-

"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
@description: contains Main wrapper, which should be used to run a standard
              GPSeq image dataset analysis.
"""

# DEPENDENCIES =================================================================

import pickle as cp
import datetime
import multiprocessing
import os
import pkg_resources
import time
import warnings

from jinja2 import Environment, FileSystemLoader
import logging
import matplotlib.pyplot as plt
import numpy as np
from numpy.lib.recfunctions import append_fields
import pandas as pd
from weasyprint import HTML, CSS

from pygpseq import const

from pygpseq.tools import io as iot
from pygpseq.tools import path as pt
from pygpseq.tools import plot
from pygpseq.tools import stat as stt
from pygpseq.tools import string as st
from pygpseq.tools import vector as vt

from pygpseq.anim.analysis import Analyzer
from pygpseq.anim.condition import Condition

# CLASSES ======================================================================


class Main(Analyzer):
    """GPSeq data manager (extends Analyzer).

    Args:
      __version__ (string): package version.
      basedir (string): input data directory path.
      outdir (string): output directory path.
      skip (list): steps to be skipped.
        1  : skip instantiation (need gpi.inst.cpickle, otherwise unskips).
        2  : skip segmentation (need gpi.seg.cpickle, otherwise unskips).
        3  : skip analysis (need gpi.an.cpickle, otherwise unskips).
        3.5: skip single-nuclei boxplot (end of step 3).
        4  : final plots.
        5  : final report.
      conds (list[pygpseq.wraps.Condition]): conditions.
      ext (string): image file extension.
      font_size (int): plot fontsize.
      ncores (int): number of cores for parallel computation.
      reg (string): regular expression for image file identification.
      debugging (bool): True for debugging mode. Debugging mode, effects:
                        Save intermediate nuclear tif images
                        of DNA, Sig, D and mask.
      plotting (bool): True to generate plot.
      suffix (string): suffix to every output file.
      notes (string): user-provided notes.
    """

    __version__ = const.VERSION
    basedir = "./"
    outdir = "./output/"
    skip = []
    conds = []
    ext = ".tif"
    reg = "^(?P<" + const.REG_CHANNEL_NAME + ">[^/]*)"
    reg += "\.(?P<" + const.REG_CHANNEL_ID + ">channel[0-9]+)"
    reg += "\.(?P<" + const.REG_SERIES_ID + ">series[0-9]+)"
    reg += "(?P<" + const.REG_EXT + ">\.tif)$"
    suffix = ""
    font_size = 8
    debugging = False
    plotting = True
    ncores = 1
    notes = "..."

    def __init__(self, ncores=None, font_size=None, **kwargs):
        """Run Analyzer __init__ method.

        Args:
          ncores (int): number of cores (opt, def: 1).
          font_size (int): plotting fontsize (opt, def: 8).
        """
        super(Main, self).__init__()

        # Read additional parameters
        if type(None) != type(ncores):
            self.ncores = ncores
        if type(None) != type(font_size):
            self.font_size = font_size

        self.printout("", 0)

    def __getitem__(self, key):
        """Allow get item."""
        if key in dir(self):
            return getattr(self, key)
        else:
            return None

    def __setitem__(self, key, value):
        """Allow set item."""
        if key in dir(self):
            self.__setattr__(key, value)

    def __setattr__(self, name, value):
        """Check the attribute and set it."""

        # Check the attribute
        check = self.check_attr(name, value)

        if check == True:
            # Set the attribute
            return super(Main, self).__setattr__(name, value)
        else:
            # Don't set the attribute
            return None

    def check_attr(self, name, value):
        """Check attribute format and value.

        Args:
          name (string): attribute name
          value: attribute value

        Returns:
          bool: whether the provided (name,value) couple passed its
                name-specific test
        """

        # Default answer
        super(Main, self).check_attr(name, value)
        checked = True

        if name in ["basedir", "outdir"]:
            # Require a string
            if type("") != type(value):
                checked = False
            # # Require an existing folder
            # elif not os.path.isdir(value):
            #   checked = False

            if not checked:
                msg = '"' + name + "\" must be an existing folder's path.\n"
                msg += "Keeping previous value."
                self.printout(msg, -1)
        elif name == "ncores":
            # Require an integer
            if type(0) != type(value):
                checked = False
            elif value <= 0:
                checekd = False

            if not checked:
                msg = '"' + name + '" must be a positive non-zero integer.\n'
                msg += "Keeping previous value."
                self.printout(msg, -1)
        elif name == "ext":
            # Require a string
            if type("") != type(value):
                checked = False
            elif not value.startswith(".") or len(value) <= 1:
                checked = False

            if not checked:
                msg = '"' + name + '" must be a file extension'
                msg += ' (start with ".").\n'
                msg += "Keeping previous value."
                self.printout(msg, -1)
        elif name == "skip":
            # Require a list
            if type([]) != type(value):
                checked = False
            elif len(value) != 0:
                # Require integer list
                types = [-1 for i in value if type(i) not in [type(0), type(0.0)]]
                types = len(types) == 0
                if not types:
                    checked = False

            if not checked:
                msg = '"' + name + '" must be an integer list (can be empty).\n'
                msg += "Keeping previous value."
                self.printout(msg, -1)
        elif name == "font_size":
            # Require positive number
            if type(value) not in [type(0), type(0.0)]:
                checked = False
            elif value <= 0:
                checked = False

            if not checked:
                msg = '"' + name + '" must be a positive non-zero number.\n'
                msg += "Keeping previous value."
                self.printout(msg, -1)
        elif name == "reg":
            # Require a string
            if type("") != type(value):
                checked = False

                msg = '"' + name + '" must be a regular expression (string).\n'
                msg += "Keeping previous value."
                self.printout(msg, -1)

        # Output
        return checked

    def is_skipped(self, step):
        """Check if a run step should be skipped."""
        return step in self.skip

    def load(self, f, more_keys=None):
        """Re-load a dumped pyGPSeq instance.
        The function keeps some parameters from the current instance.

        Args:
          f (file): file pointer to dumped instance.
          more_keys (list[string]): list of attribute names to be kept
                                    from the current instance.

        Returns:
          pygpseq.Main: dumped instance with some attributes from the current.
        """

        # Parameters to keep from current instance
        keys = list(const.PARAM_STATIC)
        if more_keys is not None:
            keys.extend(list(more_keys))
        vals = [self[k] for k in keys]

        # Load
        loaded = cp.load(f)

        # If the loaded cpickle bundle contained a Main instance
        if type(loaded) == type(self):
            # Set parameters back
            for i in range(len(keys)):
                loaded[keys[i]] = vals[i]

            # Propagate parameters
            for param in const.PARAM_PROPAGATE:
                loaded.propagate_attr(param)
        elif len(loaded) != 0:
            # Identify the Main instance
            i = [i for i in range(len(loaded)) if type(loaded[i]) == type(self)]
            tmp = loaded[i[0]]

            # Set parameters back
            for j in range(len(keys)):
                tmp[keys[j]] = vals[j]

            # Propagate parameters
            for param in const.PARAM_PROPAGATE:
                tmp.propagate_attr(param)

            # Prepare output
            loaded = [tmp if j == i else loaded[j] for j in range(len(loaded))]
            loaded = tuple(loaded)

        return loaded

    def mk_general_boxplots(self, profiles, sumd, md, **kwargs):
        """Generate general boxplots."""

        # Common destinations
        out_pdf = self.outdir + const.OUTDIR_PDF
        out_png = self.outdir + const.OUTDIR_PNG

        # Multi-condition single-nucleus boxplots
        self.printout("Preparing general single-nuclei boxplot...", 0)
        fig = plot.multi_condition_boxplot(profiles, sumd, md, **kwargs)

        # Export PDF
        fname = out_pdf + "boxplots" + kwargs["suffix"] + ".pdf"
        if self.plotting:
            plot.export(fname, "pdf")

        # Export PNG
        fname = out_png + "boxplots" + kwargs["suffix"] + ".png"
        if self.plotting:
            plot.export(fname, "png")

        # Close figure
        plt.close(fig)

        # Calculate boxplot relative widths
        bp_widths = [p["n"] for p in profiles]
        bp_widths = np.asarray(bp_widths, dtype="float")
        bp_widths = bp_widths / bp_widths.max() * 0.9

        # Per nucleus boxplots
        bpitems = [
            # Per nucleus size boxplots
            (
                plot.get_nsf_label(const.NSEL_SIZE, kwargs["seg_type"])
                + " [per nucleus]",
                "size",
            ),
            # Per nucleus shape boxplots
            (
                plot.get_nsf_label(const.NSEL_SHAPE, kwargs["seg_type"])
                + " [per nucleus]",
                "shape",
            ),
            # Per nucleus intensity average boxplots
            ("mean(DNA [a.u.]) [per nucleus]", "meanI"),
            # Per nucleus intensity sum boxplots
            ("sum(DNA [a.u.]) [per nucleus]", "sumI"),
        ]
        for (ylab, field) in bpitems:
            plot.multi_condition_single_boxplot(
                profiles,
                sumd,
                field,
                bp_widths,
                ylab,
                self.outdir + const.OUTDIR_PNG_REPORT,
                **kwargs
            )

        # Calculate boxplot relative widths
        bp_widths = np.asarray([m.shape[0] for m in md], dtype="float")
        bp_widths = bp_widths / bp_widths.max() * 0.9

        # Per pixel boxplots
        bpitems = [
            # Per pixel DNA intensity
            ("DNA [a.u.] [per pixel]", "dna"),
            # Per pixel Signal intensity
            ("Signal [a.u.] [per pixel]", "sig"),
        ]
        for (ylab, field) in bpitems:
            plot.multi_condition_single_boxplot(
                profiles,
                md,
                field,
                bp_widths,
                ylab,
                self.outdir + const.OUTDIR_PNG_REPORT,
                **kwargs
            )

        # Compare per pixel distributions --------------------------------------

        # # Add condition column to nuclear data
        # for i in range(len(md)):
        #   md[i] = append_fields(md[i],
        #       names = 'condition',
        #       data = np.tile(profiles[i]['condition'], md[i].shape[0]),
        #       dtypes = 'S100')
        # md = pd.DataFrame(vt.merge_nparrays(md))

        # # Run WMW U test
        # pvals = vt.merge_nparrays([
        #   stt.wilcox_sets(md, 'condition', 'dna'),
        #   stt.wilcox_sets(md, 'condition', 'sig')
        # ])

        # # Export as CSV
        # fname = self.outdir + const.OUTDIR_CSV
        # fname += 'single_pixel_wmw' + kwargs['suffix'] + '.csv'
        # if self.plotting: pd.DataFrame(pvals).to_csv(fname)

    def mk_general_plots(self, profiles, sumd, **kwargs):
        """Generate final plots.

        Args:
          self (pygpseq.main): current class instance.
          profiles (list): list of profile data dictionaries.

        Returns:
          Produces multiple plots and CSV files.
          profeat (list): profile features.
        """

        # Common destinations
        out_pdf = self.outdir + const.OUTDIR_PDF
        out_png = self.outdir + const.OUTDIR_PNG

        # Multi condition profile plot
        for yfield in ["mean", "median", "mode", "max"]:
            msg = "Preparing multi-condition profiles plot [" + yfield + "]..."
            self.printout(msg, 0)

            # Plot
            fig = plot.multi_condition_profiles(profiles, yfield=yfield, **kwargs)

            # Export PDF
            common_name = "profiles." + yfield + kwargs["suffix"]
            fname = out_pdf + common_name + ".pdf"
            if self.plotting:
                plot.export(fname, "pdf")

            # Export PNG
            fname = out_png + common_name + ".png"
            if self.plotting:
                plot.export(fname, "png")

            # Close figure
            plt.close(fig)

        # Export profiles to CSV
        self.printout("Exporting profiles to CSV...", 0)
        mprof = iot.merge_profiles(profiles)
        fname = self.outdir + const.OUTDIR_CSV
        fname += "profiles" + kwargs["suffix"] + ".csv"
        if self.plotting:
            pd.DataFrame(mprof).to_csv(fname)

        # Export nuclear summaries to CSV
        self.printout("Exporting summaries to CSV...", 0)
        msum = iot.merge_summaries(sumd)
        fname = self.outdir + const.OUTDIR_CSV
        fname += "summaries" + kwargs["suffix"] + ".csv"
        if self.plotting:
            pd.DataFrame(msum).to_csv(fname)

        # Background plot
        self.printout("Plotting background levels...", 0)
        plot.bgplot(self.conds, self.outdir + const.OUTDIR_PNG_REPORT, **kwargs)

    def mk_report(self, start_time, end_time, profeat, path=None):
        """Produce PDF report.

        Args:
          start_time (int): start time UNIX timestamp.
          end_time (int): end time UNIX timestamp.
          profeat (dict): profile features.
          path (string): report file path (opt).

        Returns:
          None: writes PDF report.
        """

        # From time to timestamp
        start_time = datetime.datetime.fromtimestamp(start_time)
        start_time = start_time.strftime("%Y-%m-%d %H:%M:%S")
        end_time = datetime.datetime.fromtimestamp(end_time)
        end_time = end_time.strftime("%Y-%m-%d %H:%M:%S")

        # Load template
        gpdir = pkg_resources.resource_filename(const.PACK_NAME, "")
        gpdir += "/static/"
        env = Environment(loader=FileSystemLoader(gpdir))
        env.filters["type"] = type
        env.filters["has_key"] = lambda a, b: b in a.keys()
        template = env.get_template("report_template.html")

        # Prepare variables to fill the template
        tempv = {
            "profeat": profeat,
            "starttime": start_time,
            "endtime": end_time,
            "basedir": self.basedir,
            "outdir": self.outdir,
            "logpath": self.logpath,
            "reg": self.reg,
            "ext": self.ext,
            "skip": [
                const.STEP_DESCR[i - 1]
                for i in self.skip
                if i in range(len(const.STEP_DESCR))
            ],
            "ncores": self.ncores,
            "verbose": self.verbose,
            "debugging": self.debugging,
            "suffix": self.suffix,
            "dna_names": self.dna_names,
            "sig_names": self.sig_names,
            "plotting": self.plotting,
            "fontsize": self.font_size,
            "nbins": self.nbins,
            "sigma_smooth": self.sigma_smooth,
            "sigma_density": self.sigma_density,
            "seg_type": const.SEG_LABELS[self.seg_type],
            "adp_thr": self.adaptive_neighbourhood,
            "radius_interval": self.radius_interval,
            "min_z_size": self.min_z_size,
            "offset": self.offset,
            "rm_z_tips": self.do_clear_Z_borders,
            "an_type": const.AN_LABELS[self.an_type],
            "mid_type": const.MID_SEC_LABELS[self.mid_type],
            "dist_type": self.dist_type,
            "aspect": self.aspect,
            "nsf": [
                const.NSEL_NAMES[i]
                for i in self.nsf
                if i in range(len(const.NSEL_NAMES))
            ],
            "part_n_erosion": self.part_n_erosion,
            "norm_d": self.normalize_distance,
            "rescale_deconvolved": self.rescale_deconvolved,
            "notes": self.notes,
            "conds": self.conds,
            "cnuclei": [sum(len(s.nuclei) for s in c.series) for c in self.conds],
            "cdescr": self.cdescr,
        }

        # Escape characters
        for (k, v) in tempv.items():
            if type("") == type(v):
                tempv[k] = v.replace("<", "&lt;").replace(">", "&gt;")

        # Fill template
        html_out = template.render(tempv)

        # Hide CSS warnings
        logger = logging.getLogger("weasyprint")
        logger.handlers = []
        logger.addHandler(logging.FileHandler("/tmp/weasyprint.log"))

        # Output
        suffix = datetime.datetime.fromtimestamp(time.time())
        suffix = suffix.strftime("%Y-%m-%d %H:%M:%S")
        fname = self.outdir + "report." + suffix + ".pdf"
        HTML(string=html_out).write_pdf(fname)

        # f = open(self.outdir + 'report.' + suffix + '.htm', 'w')
        # f.write(html_out)
        # f.close()

    def propagate_attr(self, key):
        """Propagate attribute current value to every condition."""
        for i in range(len(self.conds)):
            self.conds[i][key] = self[key]

    def run(self, **kwargs):
        """Run the GPSeq manager."""

        # INIT =================================================================

        # Suppress RuntimeWarning(s)
        warnings.simplefilter("ignore", category=FutureWarning)
        warnings.simplefilter("ignore", category=RuntimeWarning)
        warnings.simplefilter("ignore", category=UserWarning)

        # Check number of cores
        if self.ncores > multiprocessing.cpu_count():
            self.ncores = multiprocessing.cpu_count()
            msg = "Decreased core number to maximum allowed: %i" % self.ncores
            msg += "\nPlease, don't ask for the impossible... ಠ_ಠ"
            self.printout(msg, -1)
            self.printout("", 0)

        # Warn for log freezing if parallelization is on
        if self.ncores > 1:
            msg = "pyGPSeq log might freeze due to parallelization.\n"
            msg += "But do NOT dispair. Everything will work out, eventually..."
            msg += " Just be patient.\n┬─┬﻿ ノ( ゜-゜ノ)"
            self.printout(msg, -1)
            self.printout("", 0)

        # Check parameters
        start_time = time.time()
        kwargs["suffix"] = st.add_leading_dot(self.suffix)
        self.basedir = pt.add_trailing_slash(self.basedir)
        self.outdir = pt.add_trailing_slash(self.outdir)

        # Create basedir and outdir if missing
        for d in [self.basedir, self.outdir]:
            if not os.path.isdir(d):
                os.makedirs(d)

        # Make profile output folder
        for d in [
            const.OUTDIR_PDF,
            const.OUTDIR_PNG,
            const.OUTDIR_CSV,
            const.OUTDIR_PNG_REPORT,
            const.OUTDIR_MASK,
        ]:
            if not os.path.exists(self.outdir + d):
                os.makedirs(self.outdir + d)
        debug_dir = self.outdir + const.OUTDIR_DEBUG
        if self.debugging and not os.path.exists(debug_dir):
            os.makedirs(debug_dir)

        # Plot font-size
        if not "font_size" in kwargs.keys():
            kwargs["font_size"] = self.font_size

        # Number of cores for parallelization
        if not "ncores" in kwargs.keys():
            kwargs["ncores"] = self.ncores

        # Output directory
        kwargs["out_dir"] = self.outdir

        # Distance field for profiles
        if self.normalize_distance:
            kwargs["dfield"] = const.DLAMIN_NORM_LABEL
            kwargs["dlabel"] = "Relative distance from nuclear lamina"
            kwargs["dlabel"] += " [a.u.]"
        else:
            kwargs["dfield"] = const.DLAMIN_LABEL
            msg = "Absolute distance from nuclear lamina"
            kwargs["dlabel"] = "%s [%s]" % (msg, self.umes)

        # Update kwargs with self attributes
        kwargs.update(
            [
                (n, getattr(self, n))
                for n in dir(self)
                if type(getattr(self, n)) in const.KWARGS_TYPELIST
                and not n.startswith("__")
                and not n in const.KWARGS_AVOIDLIST
            ]
        )

        # SINGLE STEPS =========================================================

        # INSTANTIATION --------------------------------------------------------
        # Check whether to skip instantiation
        if self.is_skipped(1):
            fname = self.outdir + "gpi.inst" + kwargs["suffix"] + ".cpickle"
            if os.path.exists(fname):
                self.printout("Skipping instantiation...", 0)
                self.printout("Loading dumped instance...\n", 0)

                try:
                    with open(fname, "rb") as f:
                        keys = list(const.PARAM_SEG)
                        keys.extend(list(const.PARAM_AN))
                        self = self.load(f, keys)
                    # Dump
                    f = open(fname, "wb")
                    cp.dump(self, f)
                    f.close()
                except:
                    self.printout("Unable to load dumped instance...", 0)
                    self.printout("Unskipping instantiation...", 0)
                    self.unskip(1)
            else:
                self.printout("Unskipping instantiation...", 0)
                self.unskip(1)

        if not self.is_skipped(1):
            # Run instantiation if not skipped
            self.run_initialization(**kwargs)

            # Dump
            fname = self.outdir + "gpi.inst" + kwargs["suffix"] + ".cpickle"
            with open(fname, "wb") as f:
                cp.dump(self, f)
        # SEGMENTATION ---------------------------------------------------------
        # Check whether to skip segmentation
        if self.is_skipped(2):
            fname = self.outdir + "gpi.seg" + kwargs["suffix"] + ".cpickle"
            if os.path.exists(fname):
                self.printout("Skipping segmentation...", 0)
                self.printout("Loading dumped instance...\n", 0)

                try:
                    # Load
                    f = open(fname, "rb")
                    self = self.load(f, const.PARAM_AN)
                    f.close()

                    with open(fname, "wb") as f:
                        cp.dump(self, f)
                except:
                    self.printout("Unable to load dumped instance...", 0)
                    self.printout("Unskipping segmentation...", 0)
                    self.unskip(2)
            else:
                self.printout("Unskipping segmentation...", 0)
                self.unskip(2)

        if not self.is_skipped(2):
            # Run segmentation if not skipped
            self.run_segmentation(**kwargs)

            # Dump
            fname = self.outdir + "gpi.seg" + kwargs["suffix"] + ".cpickle"
            with open(fname, "wb") as f:
                cp.dump(self, f)
        # ANALYSIS -------------------------------------------------------------
        # Check whether to skip analysis
        if self.is_skipped(3):
            fname = self.outdir + "gpi.an" + kwargs["suffix"] + ".cpickle"
            if os.path.exists(fname):
                self.printout("Skipping analysis...", 0)
                self.printout("Loading dumped instance...\n", 0)

                try:
                    with open(fname, "rb") as f:
                        self, profiles, profeat, sumd = self.load(f)
                    with open(fname, "wb") as f:
                        cp.dump((self, profiles, profeat, sumd), f)
                except:
                    self.printout("Unable to load dumped instance...", 0)
                    self.printout("Unskipping analysis...", 0)
                    self.unskip(3)
            else:
                self.printout("Unskipping analysis...", 0)
                self.unskip(3)

        if not self.is_skipped(3):
            # Run analysis if not skipped
            profiles, profeat, sumd, md = self.run_analysis(**kwargs)

            # Dump
            fname = self.outdir + "gpi.an" + kwargs["suffix"] + ".cpickle"
            with open(fname, "wb") as f:
                cp.dump((self, profiles, profeat, sumd), f)
            # Generate general boxplots if not skipped
            if not self.is_skipped(3.5):
                self.mk_general_boxplots(profiles, sumd, md, **kwargs)

        # FINAL PLOTS ----------------------------------------------------------
        # Produce final plots if not skipped
        if not self.is_skipped(4):
            self.mk_general_plots(profiles, sumd, **kwargs)
        else:
            self.printout("Skipping final plots...", 0)

        # FINAL REPORT ---------------------------------------------------------
        # Execution end
        end_time = time.time()

        # Generate final report if not skipped
        if not self.is_skipped(5):
            self.printout("Generating final report...", 0)
            self.mk_report(start_time, end_time, profeat)
        else:
            self.printout("Skipping final report...", 0)

        # CONCLUSION ===========================================================

        # Final remarks
        self.printout("", 0)
        self.printout("~ took %s s ~" % (end_time - start_time), 0)
        self.printout("\n ~ FIN ~", 0)
        self.printout("\n└[∵┌] └[ ∵ ]┘ [┐∵]┘\n\n", 0)

        return self

    def run_analysis(self, **kwargs):
        """Run analysis."""

        # Retrieve and save nuclear data
        self.printout("* Retrieving nuclear data *", 0)
        self.printout("", 0)

        # [{dtype:{x:float, y:float}, n:int, condition:string}]
        data = [c.analyze_nuclei(**kwargs) for c in self.conds]
        profiles = [d[0] for d in data if not type(None) == type(d)]
        sumd = [d[1] for d in data if not type(None) == type(d)]
        md = [d[2] for d in data if not type(None) == type(d)]
        dp = [d[3] for d in data if not type(None) == type(d)]
        vp = [d[4] for d in data if not type(None) == type(d)]

        # Assemble and export density profile
        dp = pd.concat(dp)
        dp.to_csv(
            "%s%s/density_profiles%s.csv"
            % (kwargs["outdir"], const.OUTDIR_CSV, kwargs["suffix"])
        )

        # Assemble and export volume profile
        vp = pd.concat(vp)
        vp.to_csv(
            "%s%s/volume_profiles%s.csv"
            % (kwargs["outdir"], const.OUTDIR_CSV, kwargs["suffix"])
        )

        # Calculate profile-specific features
        self.printout("* Calculating profile features *", 0)
        profeat = []
        for ip in range(len(profiles)):
            self.printout('Working on "' + self.conds[ip].name + '"', 1)
            profile = profiles[ip]

            # Will contain current profile features
            cpf = {}

            # For profile type
            for k1 in ["dna", "sig", "ratio"]:
                self.printout('"' + k1 + '" profile...', 2)
                # Intercepts
                cis = {}

                # Area
                cia = {}

                # For statistic measure
                for k2 in ["mean", "median", "mode", "max"]:
                    cprof = profile[k1]

                    # Get the intercepts ---------------------------------------
                    cis[k2] = []

                    cis[k2].append(stt.get_intercepts(cprof[k2], cprof["x"]))
                    if len(cis[k2][0]) == len(cprof[k2]):
                        cis[k2][0] = []

                    cis[k2].append(
                        stt.get_intercepts(
                            np.diff(cprof[k2]), cprof["x"][0 : (len(cprof["x"]) - 1)]
                        )
                    )
                    if len(cis[k2][1]) == (len(cprof[k2]) - 1):
                        cis[k2][1] = []

                    cis[k2].append(
                        stt.get_intercepts(
                            np.diff(np.diff(cprof[k2])),
                            cprof["x"][0 : (len(cprof["x"]) - 2)],
                        )
                    )
                    if len(cis[k2][2]) == (len(cprof[k2]) - 2):
                        cis[k2][2] = []

                    # Get the area
                    cia[k2] = sum(cprof[k2])

                # Save current stat features
                cpf[k1] = [cis, cia]

            # Save current profile features
            profeat.append(cpf)

        return (profiles, profeat, sumd, md)

    def run_initialization(self, **kwargs):
        """Initialize run."""

        self.printout("Starting GPSeq manager...", 0)

        self.printout('Looking into folder: "' + self.basedir + '"', 0)
        self.printout("", 0)

        # Instantiate OOP architecture
        self.printout("* Building OOP Architecture *", 0)

        # Select condition folders
        self.conds = pt.select_folders(self.basedir, self.ext)

        # If no conditions, stop and trigger error
        if 0 == len(self.conds):
            msg = "No condition folders found, i.e., no tif files found in any "
            msg += "subfolders of: %s\n" % (self.basedir,)
            self.printout(msg, -2)
        else:
            msg = "Found " + str(len(self.conds)) + " condition folder(s)..."
            self.printout(msg, 0)
            self.printout("", 0)

        # Instantiate conditions
        self.conds = [
            Condition(c, self.dna_names, self.sig_names, main=self) for c in self.conds
        ]
        self.printout("", 0)

    def run_segmentation(self, **kwargs):
        """Run segmentation."""

        # Check analysis/segmentation types
        self.check_anseg_types()
        kwargs["seg_type"] = self.seg_type
        kwargs["an_type"] = self.an_type

        # Identify nuclei
        self.printout("* Looking for nuclei *", 0)
        self.printout("", 0)
        [c.find_nuclei(**kwargs) for c in self.conds]
        self.printout("", 0)

    def unskip(self, step):
        """Unskips a run step that was supposed to be skipped."""
        self.skip = [i for i in self.skip if step != i]


# END ==========================================================================

################################################################################

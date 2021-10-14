# -*- coding: utf-8 -*-

"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
@description: contains Condition wrapper, which in turn contains Series.
"""

# DEPENDENCIES =================================================================

from joblib import Parallel, delayed
import multiprocessing
import os
import time

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from pygpseq import const
from pygpseq.tools import path as pt, io as iot, plot, stat as stt, string as st

from pygpseq.anim.series import Series

# CLASSES ======================================================================


class Condition(iot.IOinterface):
    """GPSeq condition, i.e., ensemble of series with nuclei.

    Args:
      __version__ (string): package string.
      path (string): condition folder path.
      name (string): condition name.
      ext (string): condition series extension.
      reg (string): condition series regexp.
      series (list[series]): condition series.
    """

    __version__ = const.VERSION
    path = "."
    name = ""
    ext = ".tif"
    reg = "^(?P<channel_name>[^/]*)"
    reg += "\.(?P<channel_str>channel[0-9]+)"
    reg += "\.(?P<series_str>series[0-9]+)"
    reg += "(?P<ext>\.tif)$"
    series = []

    def __init__(self, path, dna_channels, sig_channels, main=None):
        """Run IOinterface __init__ method.

        Args:
          path (string): path to the condition folder.
          main (pyGPSeq.main): main wrapper (opt).
        """

        # If required, inherit from `main` wrap
        if main is None:
            super(Condition, self).__init__()

        else:
            logpath = main.logpath
            super(Condition, self).__init__(path=logpath, append=True)
            self.ext = main.ext
            self.verbose = main.verbose
            self.reg = main.reg
        # Save input parameters
        self.path = pt.add_trailing_slash(os.path.abspath(path))
        self.name = self.path[: len(self.path) - 1].split("/")
        self.name = self.name[len(self.name) - 1]
        self.printout('Initializing condition: "' + self.name + '"', 0)

        # Select condition's series
        self.series = pt.select_files(self.path, self.ext)
        self.series = pt.select_series(self.series, self.reg).items()

        # Check that each series has at least one dna_channel and sig_channel
        for s in self.series:
            if all(
                cdata["channel_name"] not in dna_channels
                for (cname, cdata) in s[1].items()
            ):
                self.printout(
                    "No DNA channel found in '%s' of '%s'" % (s[0], self.name), -2
                )
            if all(
                cdata["channel_name"] not in sig_channels
                for (cname, cdata) in s[1].items()
            ):
                self.printout(
                    "No Signal channel found in '%s' of '%s'" % (s[0], self.name), -2
                )

        # If no series, stop and trigger error
        if len(self.series) == 0:
            msg = "No series found in condition %s." % (self.name,)
            self.printout(msg, -2)
        else:
            msg = "Found %d series..." % (len(self.series),)
            self.printout(msg, 1)

        # Instantiate series
        self.series = [list(ds) for ds in self.series]
        [self.series[i].append(i + 1) for i in range(len(self.series))]
        self.series = [Series(s, condition=self) for s in self.series]

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

    def adjust_options(self, **kwargs):
        """Adjust options to be passed to the Series class.

        Args:
          **kwargs

        Returns:
          dict: adds the following kawrgs:
            cond_name (string): condition wrapper name.
        """

        kwargs["cond_name"] = self.name

        return kwargs

    def analyze_nuclei(self, **kwargs):
        """Export current condition nuclei.

        Args:
          sigma (float): sigma for smoothing and covariance calculation (opt).
          **kwargs

        Returns:
          tuple: profiles, summaries and single-pixel tables.
        """

        # CHECK PARAMS =========================================================

        # Get default number of cores
        ncores = 1 if not "ncores" in kwargs.keys() else kwargs["ncores"]
        # Set output suffix
        if not "suffix" in kwargs.keys():
            suffix = ""
        else:
            suffix = st.add_leading_dot(kwargs["suffix"])

        # Check plotting
        if not "plotting" in kwargs.keys():
            kwargs["plotting"] = True

        # Check number of cores
        if ncores > multiprocessing.cpu_count():
            ncores = multiprocessing.cpu_count()
            msg = "Decreased core number to maximum allowed: %i" % ncores
            msg += "\nPlease, don't ask for the impossible... ಠ_ಠ"
            self.printout(msg, -1)

        # Add necessary options
        self.printout('Current condition: "' + self.name + '"...', 0)
        kwargs = self.adjust_options(**kwargs)

        # Create condition nuclear data directory if necessary
        if not os.path.isdir(kwargs["out_dir"]):
            os.mkdir(kwargs["out_dir"])

        # GET NUCLEAR DATA =====================================================

        # Retrieve nuclei
        nuclei = self.get_nuclei()

        # Check that the condition contains nuclei
        if len(nuclei) == 0:
            return (None, None, None, None)

        if kwargs["seg_type"] == const.SEG_3D:
            DTYPE_NUCLEAR_SUMMARY = const.DTYPE_NUCLEAR_SUMMARY_3D
        else:
            DTYPE_NUCLEAR_SUMMARY = const.DTYPE_NUCLEAR_SUMMARY_2D

        # Retrieve nuclei summaries
        self.printout("Retrieving nuclear summary...", 1)
        summary = np.zeros(len(nuclei), dtype=DTYPE_NUCLEAR_SUMMARY)
        for i in range(len(nuclei)):
            summary[i] = nuclei[i].get_summary()

        # Filter nuclei
        msg = "Filtering nuclei based on size, intensity and shape..."
        self.printout(msg, 1)
        selected = self.multi_threshold_nuclei(data=summary, **kwargs)

        # Check that nuclei are selected
        if len(selected) == 0:
            return (None, None, None, None)

        # Apply selection
        summary = np.asarray(
            [summary[i] for i in selected], dtype=DTYPE_NUCLEAR_SUMMARY
        )

        # Retrieve selected nuclei single-pixel data
        data_nested = Parallel(n_jobs=ncores)(
            delayed(get_series_nuclear_data)(self, summary, sidx, **kwargs)
            for sidx in list(set(summary["s"]))
        )

        # Un-nest nuclear data
        data = []
        [data.extend(nested["spx_data"]) for nested in data_nested]

        # Assemble into a single array
        self.printout("Merging into a single table...", 1)
        merged = np.zeros(
            sum([d.shape[0] for d in data]), dtype=const.DTYPE_NUCLEAR_DATA
        )
        currpos = 0
        for d in data:
            merged[currpos : (currpos + d.shape[0])] = d
            currpos += d.shape[0]

        # Remove rows with no DNA signal
        self.printout("Removing pixels without DNA signal...", 1)
        self.printout("Identifying pixels...", 2)
        toKeep = np.where(merged["dna"] != 0)[0]
        nToRemove = len(merged) - len(toKeep)
        if not nToRemove == 0:
            merged = merged[toKeep]
            msg = "Removed %i pixels without DNA signal..." % nToRemove
            self.printout(msg, 2)

        # Density profile ------------------------------------------------------

        dp = np.vstack([nested["density"] for nested in data_nested])
        dp = pd.DataFrame(dp)
        col_labs = ["c", "s", "n"]
        col_labs.extend(
            ["nd_%f" % b for b in np.linspace(0, 1, kwargs["nbins"] + 1)[1:]]
        )
        dp.columns = col_labs

        # Volume profile -------------------------------------------------------

        vp = np.vstack([nested["volume"] for nested in data_nested])
        vp = pd.DataFrame(vp)
        col_labs = ["c", "s", "n"]
        col_labs.extend(
            ["nd_%f" % b for b in np.linspace(0, 1, kwargs["nbins"] + 1)[1:]]
        )
        vp.columns = col_labs

        # PLOT =================================================================

        # EVERY PIXEL ----------------------------------------------------------

        # Produce profile plot
        self.printout("Generating profiles...", 1)
        profiles = self.make_profiles(merged, len(data), **kwargs)

        # Export single profile study
        self.printout("Studying single-pixel behaviour...", 1)
        self.check_single_pixels(merged, profiles, **kwargs)

        # Export single-condition plot
        self.printout("Exporting profiles...", 1)

        # Mean/median/mode profile plot
        fig = plot.single_condition_profiles(profiles, n_nuclei=len(data), **kwargs)
        plot.single_condition_profiles(
            profiles, n_nuclei=len(data), yfield="median", new_figure=False, **kwargs
        )
        plot.single_condition_profiles(
            profiles, n_nuclei=len(data), yfield="mode", new_figure=False, **kwargs
        )
        plot.single_condition_profiles(
            profiles, n_nuclei=len(data), yfield="max", new_figure=False, **kwargs
        )

        # Add legend
        plt.subplot(3, 2, 1)
        plot.set_font_size(12)
        plt.legend(
            labels=["mean", "median", "mode", "max"],
            bbox_to_anchor=(0.0, 1.12, 1.0, 0.102),
            loc=3,
            ncol=2,
            mode="expand",
            borderaxespad=0.0,
        )

        # Export PDF
        fname = kwargs["out_dir"] + const.OUTDIR_PDF + self.name
        fname += ".profiles" + suffix + ".pdf"
        if kwargs["plotting"]:
            plot.export(fname, "pdf")

        # Export PNG
        fname = kwargs["out_dir"] + const.OUTDIR_PNG + self.name
        fname += ".profiles" + suffix + ".png"
        if kwargs["plotting"]:
            plot.export(fname, "png")

        # Close figure
        plt.close(fig)

        # Output
        self.printout("", 0)
        return (profiles, summary, merged, dp, vp)

    def check_single_pixels(
        self, indata, profiles, partial=None, supcomm=None, **kwargs
    ):
        """Produce single pixel behaviour study plot.

        Args:
          indata (np.array): single-pixel table, const.DTYPE_NUCLEAR_DATA.
          profiles (dict): smoothened and raw profiles (I ~ d).
          partial (bool): True if working on partial volume.
          supcomm (string): a comment to be add to the plot main title.
          **kwargs
        """

        # CHECK PARAMS =========================================================

        # Set output suffix
        if not "suffix" in kwargs.keys():
            suffix = ""
        else:
            suffix = st.add_leading_dot(kwargs["suffix"])

        # Partial volume data
        if None == partial:
            partial = False
        partial = False

        # Check plotting
        if not "plotting" in kwargs.keys():
            kwargs["plotting"] = True

        # Output file pointers
        fname = kwargs["out_dir"] + const.OUTDIR_PDF
        out_png = kwargs["out_dir"] + const.OUTDIR_PNG
        out_png += self.name + ".pixel_study."
        if partial:
            fname += self.name + ".pixel_study.part" + suffix + ".pdf"
        else:
            fname += self.name + ".pixel_study" + suffix + ".pdf"
        if kwargs["plotting"]:
            pp = PdfPages(fname)

        # PREPARE DATA =========================================================

        # Setup data for plotting
        dna = indata["dna"].astype("float")
        rat = indata["sig"] / dna
        pltitems = [
            ("DNA channel...", indata["dna"], "DNA [a.u.]", "dna"),
            ("Signal channel...", indata["sig"], "Signal [a.u.]", "sig"),
            ("Signal/DNA ratio...", rat[rat != np.inf], "Signal/DNA", "ratio"),
        ]

        # PLOT =================================================================

        # Set plot super title
        if self.name in kwargs["cdescr"].keys():
            suptitle = 'Analysis "' + kwargs["cdescr"][self.name] + '"'
        else:
            suptitle = 'Analysis "' + self.name + '"'
        if None != supcomm:
            suptitle += supcomm
        suptitle += " [" + str(kwargs["an_type"]) + "]"
        suptitle += " [sigma = " + str(kwargs["sigma_smooth"]) + "]"
        suptitle += " [nbins = " + str(kwargs["nbins"]) + "]"

        # Plot
        for (msg, y, ylab, lab) in pltitems:
            self.printout(msg, 2)

            # Actual plot
            fig = plot.single_pixel_study(
                indata[kwargs["dfield"]],
                y,
                ylab,
                profiles[lab],
                partial=partial,
                **kwargs
            )
            fig.tight_layout()
            plt.subplots_adjust(top=0.95)
            plt.suptitle(suptitle)

            # Export PDF
            if kwargs["plotting"]:
                plt.savefig(pp, format="pdf")

            # Export PNG
            if partial:
                fname = out_png + lab + ".part" + suffix + ".png"
            else:
                fname = out_png + lab + suffix + ".png"
            if kwargs["plotting"]:
                plt.savefig(fname, format="png")

            # Close figure
            plt.close(fig)

        # Close file pointer
        if kwargs["plotting"]:
            pp.close()

    def export_nuclei(self, **kwargs):
        """Export current condition nuclei."""

        # Set output suffix
        if not "suffix" in kwargs.keys():
            suffix = ""
        else:
            suffix = st.add_leading_dot(kwargs["suffix"])

        # Add necessary options
        self.printout('Current condition: "' + self.name + '"...', 0)
        kwargs = self.adjust_options(**kwargs)

        # Create condition nuclear data directory if necessary
        if not os.path.isdir(kwargs["out_dir"]):
            os.mkdir(kwargs["out_dir"])

        # Segment every series in the condition
        logs = [s.export_nuclei(**kwargs) for s in self.series]

        if kwargs["seg_type"] == const.SEG_3D:
            DTYPE_NUCLEAR_SUMMARY = const.DTYPE_NUCLEAR_SUMMARY_3D
        else:
            DTYPE_NUCLEAR_SUMMARY = const.DTYPE_NUCLEAR_SUMMARY_2D

        # Will contain the summaries
        summary = np.zeros(sum([i.shape[0] for i in logs]), dtype=DTYPE_NUCLEAR_SUMMARY)

        # Log counter
        c = 0
        for log in logs:
            # Merge logs
            summary[c : (c + log.shape[0]), :] = log

            # Increase log counter
            c += log.shape[0]

        # Export condition summary
        np.savetxt(
            kwargs["out_dir"] + "summary" + suffix + ".csv",
            summary,
            delimiter=",",
            comments="",
            header=",".join([h for h in summary.dtype.names]),
        )

    def find_nuclei(self, **kwargs):
        """Segment current condition.

        Args:
          **kwargs: all Main attributes.
        """

        # Get default number of cores
        if not "ncores" in kwargs.keys():
            ncores = 1
        else:
            ncores = kwargs["ncores"]

        # Suffix for output
        if not "suffix" in kwargs.keys():
            suffix = ""
        else:
            suffix = st.add_leading_dot(kwargs["suffix"])

        # Check plotting
        if not "plotting" in kwargs.keys():
            kwargs["plotting"] = True

        # Check number of cores
        if ncores > multiprocessing.cpu_count():
            ncores = multiprocessing.cpu_count()
            msg = "Decreased core number to maximum allowed: %i" % ncores
            msg += "\nPlease, don't ask for the impossible... ಠ_ಠ"
            self.printout(msg, -1)

        # Add necessary options
        kwargs = self.adjust_options(**kwargs)

        # Segment every series in the condition
        self.printout('Current condition: "' + self.name + '"...', 0)
        self.series = Parallel(n_jobs=ncores)(
            delayed(find_series_nuclei)(self, i, **kwargs)
            for i in range(len(self.series))
        )

    def get_nuclei(self):
        """Return a list of the nuclei in the condition."""
        nuclei = []
        for s in self.series:
            nuclei.extend(s.nuclei)
        return nuclei

    def make_profiles(self, pdata, n_nuclei, **kwargs):
        """Prepare profiles for plotting.

        Args:
          pdata (np.array): single-pixel data table, const.DTYPE_NUCLEAR_DATA.
          n_nuclei (int): number of nuclei.
          **kwargs

        Returns:
          list: profiles.
        """

        # Will contain the profiles
        profiles = {}

        # DNA profile
        self.printout("DNA profile...", 2)
        profiles["dna"] = stt.smooth_sparse_gaussian(
            pdata[kwargs["dfield"]].tolist(), pdata["dna"].tolist(), **kwargs
        )

        # Signal profile
        self.printout("Signal profile...", 2)
        profiles["sig"] = stt.smooth_sparse_gaussian(
            pdata[kwargs["dfield"]].tolist(), pdata["sig"].tolist(), **kwargs
        )

        # Ratio profile
        self.printout("Signal/DNA profile...", 2)
        rat = pdata["sig"] / pdata["dna"].astype("float")
        profiles["ratio"] = stt.smooth_sparse_gaussian(
            pdata[kwargs["dfield"]][rat != np.inf].tolist(),
            rat[rat != np.inf].tolist(),
            **kwargs
        )

        # Save number of nuclei and condition name
        profiles["n"] = n_nuclei
        profiles["condition"] = self.name

        # Output
        return profiles

    def multi_threshold_nuclei(
        self,
        cond_name,
        seg_type,
        data,
        sigma_density,
        nsf,
        font_size,
        out_dir,
        wspace=None,
        hspace=None,
        **kwargs
    ):
        """Plot density with FWHM range and select rows in that range.
        Features used for the selection based on nsf.

        Args:
          seg_type (pygpseq.const): segmentation type according to pygpseq.const
          data (numpy.array): nuclear data.
          sigma_density (float): sigma for density calculation.
          xsteps (int): density distribution curve precision (opt).
          nsf (list[int]): list of features to use for nuclear selection
                           according to pygpseq.const.

        Returns:
          list: list of selected nuclei indexes.
        """

        # Set output suffix
        if not "suffix" in kwargs.keys():
            suffix = ""
        else:
            suffix = st.add_leading_dot(kwargs["suffix"])

        # Check plotting
        if not "plotting" in kwargs.keys():
            kwargs["plotting"] = True

        fig = plt.figure(figsize=[8, 8])
        suptitle = "Automatic nuclei threshold for condition "
        if self.name in kwargs["cdescr"].keys():
            suptitle += '"%s"' % (kwargs["cdescr"][cond_name],)
        else:
            suptitle += '"%s"' % (cond_name,)
            suptitle += "\n [sigma: %.2f]" % (sigma_density,)
        plt.suptitle(suptitle)

        # Setup subplots spacing
        if None == wspace:
            wspace = 0.4
        if None == hspace:
            hspace = 0.4
        plt.subplots_adjust(wspace=wspace, hspace=hspace)

        if not 0 == len(nsf):
            # Filter features
            sel_data = {}
            plot_counter = 1
            for nsfi in nsf:
                # Identify Nuclear Selection Feature
                nsf_field = const.NSEL_FIELDS[nsfi]
                nsf_name = const.NSEL_NAMES[nsfi]
                self.printout("Filtering " + nsf_name + "...", 2)

                # Select subplot
                if 1 == len(nsf):
                    plt.subplot(1, 1, 1)
                else:
                    plt.subplot(2, 2, plot_counter)

                # Plot
                sel_data[nsf_field] = self.single_threshold_nuclei(
                    data=data[nsf_field],
                    sigma_density=sigma_density,
                    xlab=plot.get_nsf_label(nsfi, seg_type),
                )

                # Edit plot
                plot.set_font_size(font_size)
                plot.density_with_range(new_figure=False, **sel_data[nsf_field])
                plot_counter += 1

            # Select based on range
            self.printout("Selecting nuclei...", 2)
            f = lambda x, r: x >= r[0] and x <= r[1]
            for nsfi in nsf:
                nsf_field = const.NSEL_FIELDS[nsfi]

                # Identify nuclei in the FWHM range
                nsf_data = sel_data[nsf_field]
                nsf_data["sel"] = [
                    f(i, nsf_data["fwhm_range"]) for i in nsf_data["data"]
                ]
                sel_data[nsf_field] = nsf_data

            # Select those in every FWHM range
            nsfields = [const.NSEL_FIELDS[nsfi] for nsfi in nsf]
            selected = [sel_data[f]["sel"] for f in nsfields]
            g = lambda i: all([sel[i] for sel in selected])
            selected = [i for i in range(len(selected[0])) if g(i)]
            sub_data = data[selected]

            # Set title
            title = "Selected " + str(len(selected)) + "/"
            title += str(len(sel_data[const.NSEL_FIELDS[nsf[0]]]["data"]))
            title += " nuclei."

            if not 1 == len(nsf):
                # General scatterplot
                plt.subplot(2, 2, 4)
                plot.set_font_size(font_size)
                plt.plot(
                    data[const.NSEL_FIELDS[nsf[0]]],
                    data[const.NSEL_FIELDS[nsf[1]]],
                    ",k",
                )
                plt.hold(True)

                # Selected scatterplot
                plt.plot(
                    sub_data[const.NSEL_FIELDS[nsf[0]]],
                    sub_data[const.NSEL_FIELDS[nsf[1]]],
                    ",r",
                )
                plt.xlabel(plot.get_nsf_label(nsf[0], seg_type))
                plt.ylabel(plot.get_nsf_label(nsf[1], seg_type))
                plt.title(title)
                plt.ticklabel_format(style="sci", axis="x", scilimits=(0, 0))
                plt.ticklabel_format(style="sci", axis="y", scilimits=(0, 0))
        else:
            # Set title
            selected = range(data.shape[0])
            title = "Selected " + str(len(selected)) + " nuclei."
            plt.title(title)

        # Export plot
        self.printout("Exporting threshold plot...", 2)
        fname = out_dir + const.OUTDIR_PDF + self.name + ".threshold_summary"
        fname += suffix + ".pdf"
        if kwargs["plotting"]:
            plot.export(fname, "pdf")
        fname = out_dir + const.OUTDIR_PNG + self.name + ".threshold_summary"
        fname += suffix + ".png"
        if kwargs["plotting"]:
            plot.export(fname, "png")
        self.printout(title, 2)
        plt.close(fig)

        # Output
        return selected

    def propagate_attr(self, key):
        """Propagate attribute current value to every series."""
        for i in range(len(self.series)):
            self.series[i][key] = self[key]

    def single_threshold_nuclei(self, data, sigma_density, xlab=None, **kwargs):
        """Select a single-feature nuclear threshold.

        Args:
          data (numpy.array): single column nuclear data.
          sigma_density (float): sigma for density profile calculation.
          xlab (string): x-label (opt).

        Returns:
          dict: data to generate the density threshold plot.
        """

        # Start building output
        t = {"data": data}

        # Calculate density
        t["density"] = stt.calc_density(t["data"], sigma=sigma_density)

        # Identify range
        args = [t["density"]["x"], t["density"]["y"]]
        t["fwhm_range"] = stt.get_fwhm(*args)

        # Add plot features
        if None != xlab:
            t["xlab"] = xlab
        else:
            t["xlab"] = "x"
        t["ylab"] = "Density"
        np.set_printoptions(precision=3)
        t["title"] = str(np.array([float(x) for x in t["fwhm_range"]]))

        # Output
        return t


# FUNCTIONS ====================================================================


def find_series_nuclei(self, i, **kwargs):
    """Function for parallelized nuclear segmentation.

    Args:
      i (int): series index.
      **kwargs: ncores.

    Returns:
      pygpseq.wraps.Series: updated i-th series instance.
    """

    # Set series verbosity
    if kwargs["ncores"] != 1:
        self.series[i].verbose = False

    # Get starting time
    start_time = time.time()

    # Find nuclei
    self.series[i], log = self.series[i].find_nuclei(**kwargs)

    # Print log all at once
    time_msg = "Took %s s.\n" % (round(time.time() - start_time, 3))
    if not 1 == kwargs["ncores"]:
        log += iot.printout(time_msg, 1, False)
        self.printout(log, 0)
    else:
        self.printout(time_msg, 1)

    # Output
    return self.series[i]


def get_series_nuclear_data(self, summary, sidx, **kwargs):
    """Function for parallelized single-pixel nuclear data retrieval.

    Args:
      summary (list): list of summaries.
      sidx (int): series index.
      **kwargs: ncores.

    Returns:
      np.array: series nuclear data.
    """

    # Set series verbosity
    if kwargs["ncores"] != 1:
        self.series[sidx - 1].verbose = False
    else:
        self.series[sidx - 1].verbose = True

    # Get starting time
    start_time = time.time()

    # Setup starting message
    msg = "Retrieving nuclear data from series #" + str(sidx) + "..."
    msg = iot.printout(msg, 1, verbose=False)

    # Get nuclei ids for current series
    ns = [i for i in range(summary.shape[0]) if summary["s"][i] == sidx]
    ns = summary["n"][ns]

    # Retrieve nuclear data
    data, dp, vp, log = self.series[sidx - 1].get_nuclei_data(ns, **kwargs)

    # Print log all at once
    time_msg = "Took %s s." % (round(time.time() - start_time, 3))
    if not 1 == kwargs["ncores"]:
        log = msg + log
        log += iot.printout(time_msg, 2, False)
        self.printout(log, 0)
    else:
        self.printout(time_msg, 2)

    # Output
    return {"spx_data": data, "density": dp, "volume": vp}


# END ==========================================================================

################################################################################

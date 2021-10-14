# -*- coding: utf-8 -*-

"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
@description: plot generation library.
@notes: If the font cache is re-built every time that matplotlib is imported,
        try removing the ~/.cache/matplotlib/ folder.
"""

# DEPENDENCIES =================================================================

import matplotlib
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import rc
from matplotlib import ticker
from matplotlib.ticker import ScalarFormatter
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import numpy as np
import os
import scipy
import seaborn as sns
import tifffile

from pygpseq import const

from pygpseq.tools import path as pt, stat as stt

# FUNCTIONS ====================================================================


def bgplot(conds, out_png, **kwargs):
    """Generate background intensity boxplot.
    The produced plot is directly exported in PNG

    Args:
      conds (list): list of pygpseq.Condition objects.
      out_png (string): path to png output folder.
      **kwargs: optional keyword arguments (suffix).
    """

    # Cycle through the FoVs and retrieve the background level -----------------

    # Empty list that will contain the backgrounds
    bgs = []

    # For every condition
    for ci in range(len(conds)):
        # Current condition
        cc = conds[ci]

        # For every series in the current condition
        for si in range(len(cc.series)):
            # Current series
            cs = cc.series[si]

            # Append
            bgs.append((ci, si, cs.dna_bg, cs.sig_bg))

    # Convert to np.array
    bgs = np.array(
        bgs, dtype=[("c", "u4"), ("s", "u4"), ("dnabg", "u8"), ("sigbg", "u8")]
    )

    # Plot ---------------------------------------------------------------------

    # Store the median-median background  for the different conditions
    dm = []

    # Iterate on the conditions
    for cidx in range(len(conds)):

        # Count FoVs in the current condition
        N = sum(bgs["c"] == cidx)
        ind = np.arange(N)
        width = 0.35

        # Prepare plot
        fig, ax = plt.subplots()
        dna_bg = bgs["dnabg"][bgs["c"] == cidx]
        rects1 = ax.bar(ind, dna_bg, width, color="b")
        sig_bg = bgs["sigbg"][bgs["c"] == cidx]
        rects2 = ax.bar(ind + width, sig_bg, width, color="g")

        # Add medians
        plt.axhline(y=np.median(dna_bg), xmin=0, xmax=N, color="b")
        plt.axhline(y=np.median(sig_bg), xmin=0, xmax=N, color="g")

        # Store medians
        dm.append((cidx, np.median(dna_bg), np.median(sig_bg)))

        # Labels and legend
        ax.set_ylabel("Intensity [a.u.]")
        title = 'Background levels per FoV in condition "'
        if conds[cidx].name in kwargs["cdescr"].keys():
            title += kwargs["cdescr"][conds[cidx].name] + '"'
        else:
            title += conds[cidx].name + '"'
        ax.set_title(title)
        ax.set_xlabel("Field of view")
        ax.set_xticks(ind + width)
        ax.set_xticklabels([str(i + 1) for i in range(N)])
        ax.legend((rects1[0], rects2[0]), ("DNA", "Signal"), loc="best")

        # Save plot
        fname = out_png + conds[cidx].name + ".backgrounds"
        if "suffix" in kwargs.keys():
            fname += kwargs["suffix"]
        fname += ".png"
        export(fname, "png")
        plt.close(fig)

    # Plot the summary
    dm = np.array(dm, dtype=[("c", "u4"), ("dnabg", "u8"), ("sigbg", "u8")])

    # Count FoVs in the current condition
    N = len(conds)
    ind = np.arange(N)
    width = 0.35

    # Prepare plot
    fig, ax = plt.subplots()
    dna_bg = dm["dnabg"]
    rects1 = ax.bar(ind, dna_bg, width, color="b")
    sig_bg = dm["sigbg"]
    rects2 = ax.bar(ind + width, sig_bg, width, color="g")

    # Add medians
    plt.axhline(y=np.median(dna_bg), xmin=0, xmax=N, color="b")
    plt.axhline(y=np.median(sig_bg), xmin=0, xmax=N, color="g")

    # Labels and legend
    ax.set_ylabel("Intensity [a.u.]")
    title = "Background levels per condition"
    ax.set_title(title)
    ax.set_xlabel("Condition")
    ax.set_xticks(ind + width)
    xticks = []
    for cidx in range(len(conds)):
        cn = conds[cidx].name
        if cn in kwargs["cdescr"].keys():
            xticks.append(kwargs["cdescr"][cn])
        else:
            xticks.append(cn)
    ax.set_xticklabels(xticks, rotation=90)
    ax.legend((rects1[0], rects2[0]), ("DNA", "Signal"), loc="best")

    # Save plot
    plt.tight_layout()
    fname = out_png + "backgrounds"
    if "suffix" in kwargs.keys():
        fname += kwargs["suffix"]
    fname += ".png"
    export(fname, "png")
    plt.close(fig)


def boxplot(
    fig,
    subplot,
    data,
    labels,
    widths,
    ylabel,
    nnuclei,
    xaxis1=None,
    xaxis2=None,
    ylab_pos=None,
    forceSciNotation=None,
):
    """Generate single condition boxplot.
    Plot in-place in the provided figure canvas object.

    Args:
      fig (plt.figure): figure canvas.
      subplot (int): subplot signature (i.e., 111, 122,...).
      data (list): list of datasets to boxplot.
      labels (list[string]): x-axis #1 tick values.
      widths (list[float]): boxplot widths.
      ylabel (string): y-axis label.
      nnuclei (list(int)): x-axis #2 tick values.
      xaxis1 (bool): True to plot x-axis #1 (opt, def True).
      xaxis2 (bool): True to plot x-axis #2 (opt, def True).
      ylab_pos (string): y-axis position (opt, def 'left').
      forceSciNotation (bool): True to force sci. notation (opt, def False).
    """

    # SET PARAMS ===============================================================

    if xaxis1 is None:
        xaxis1 = True

    if xaxis2 is None:
        xaxis2 = True

    if ylab_pos is None:
        ylab_pos = "left"

    if forceSciNotation is None:
        forceSciNotation = False

    # PLOT =====================================================================

    # Create double xaxis
    ax1 = fig.add_subplot(subplot)
    ax2 = ax1.twiny()

    # Produce boxplot
    bp = ax1.boxplot(data, labels=labels, sym=".", widths=widths, patch_artist=True)

    # Change boxplot appearance
    [b.set(color="#000000") for b in bp["boxes"]]
    [b.set(facecolor="#ffffff") for b in bp["boxes"]]
    [w.set(color="#000000") for w in bp["whiskers"]]
    [c.set(color="#000000") for c in bp["caps"]]
    [m.set(color="red") for m in bp["medians"]]
    [m.set(linewidth=2) for m in bp["medians"]]
    [f.set(color="#000000") for f in bp["fliers"]]

    # Change axes labels
    ax1.set_xlabel("Condition")
    ax1.set_ylabel(ylabel)
    ax2.set_xlim(ax1.get_xlim())
    ax2.set_xticks(np.asarray(range(len(widths))) + 1)
    ax2.set_xticklabels(nnuclei)
    ax2.set_xlabel("Number of nuclei")

    # Force scientific notation
    if forceSciNotation:
        ax1.ticklabel_format(style="sci", axis="y", scilimits=(0, 0))

    # Move y axis
    ax1.yaxis.set_label_position(ylab_pos)

    # Hide unnecessary axes
    if not xaxis1:
        ax1.get_xaxis().set_visible(False)
    if not xaxis2:
        ax2.get_xaxis().set_visible(False)


def density_with_range(
    density, fwhm_range, new_figure=None, show=None, close=None, **kwargs
):
    """Plot density curve and FWHM range.
    Density should be the output of .stat.calc_density

    Args:
      density (dict): dictionary w/ density function and coords for its plot.
      fwhm_range (list): FWHM range of major density peak.
      new_figure (plt.figure): figure object to plot to.
      show (bool): True to display the figure.
      close (bool): True to close the open canvas.
      **kwargs: optional keyword arguments (xlab, ylab, title)
    """

    # Default values
    if new_figure is None:
        new_figure = True
    if show is None:
        show = True
    if close is None:
        clos = False

    # New figure
    if new_figure:
        fig = plt.figure(figsize=[8, 8])

    # Plot
    plt.plot(density["x"], density["y"], "k")
    plt.hold(True)
    plt.plot(fwhm_range, density["f"](fwhm_range), "b.")
    plt.axvline(x=fwhm_range[0], color="red", linestyle="-.")
    plt.axvline(x=fwhm_range[1], color="red", linestyle="-.")
    plt.ticklabel_format(style="sci", axis="x", scilimits=(0, 0))
    plt.ticklabel_format(style="sci", axis="y", scilimits=(0, 0))

    # Adjust labels
    if "xlab" in kwargs.keys():
        plt.xlabel(kwargs["xlab"])
    if "ylab" in kwargs.keys():
        plt.ylabel(kwargs["ylab"])
    if "title" in kwargs.keys():
        plt.title(kwargs["title"])

    # Show figure
    if new_figure and show:
        plt.show()

    # Close
    if new_figure and close:
        plt.close(fig)


def dots_in_ellipsoid(
    a, b, c, coords, aspect=None, fold=None, channels=None, title=None, outpng=None
):
    """Plot coords in an ellipsoid, 3 orthogonal views (XT, YZ, XZ).

    Args:
        a, b, c (float): ellipsoid semi-axes.
        coords (np.ndarray): X/Y/Z coordinates array.
        fold (bool): whether to fold X/Y dimensions.
        aspect (tuple): Z/Y/X aspect ratio.
        channels (list): list of channel names for each dot in coords.
        title (str): title of the full image.
        outpng (str): path to output png file (optional).
    """

    # Set defaults -------------------------------------------------------------
    if type(None) == type(title):
        title = ""
    if type(None) == type(fold):
        fold = False
    if type(None) == type(channels):
        channels = np.repeat("b", coords.shape[1])
    if type(None) == type(aspect):
        aspect = (1, 1, 1)
    az, ay, ax = aspect

    if type(None) != type(outpng):
        outpng = open(outpng, "wb")
        plt.close("all")

    # Plot ---------------------------------------------------------------------
    plt.figure(figsize=(10, 10))

    # Extract coordinates
    x, y, z = coords

    # Calculate general axes limits
    dlim = np.array([a * 1.1, b * 1.1, c * 1.1]).max()

    # Prepare colors
    cmap = plt.get_cmap("rainbow")
    channel_names = list(set(channels))
    channel_names.sort()
    channel_colors = cmap(np.linspace(0, 1, len(channel_names)))

    # XY plot
    st = plt.subplot(2, 2, 1)
    e = Ellipse((0, 0), 2 * a, 2 * b, facecolor="none", edgecolor="k")
    st.add_artist(e)
    if fold:
        for i in range(len(channel_names)):
            plt.plot(
                abs(x)[channels == channel_names[i]],
                abs(y)[channels == channel_names[i]],
                ".",
                color=channel_colors[i],
            )
        plt.xlim(0, dlim)
        plt.ylim(0, dlim * ax / ay)
    else:
        for i in range(len(channel_names)):
            plt.plot(
                x[channels == channel_names[i]],
                y[channels == channel_names[i]],
                ".",
                color=channel_colors[i],
            )
        plt.xlim(-dlim, dlim)
        plt.ylim(-dlim, dlim)
    plt.xlabel("X")
    plt.ylabel("Y")

    # ZY plot
    st = plt.subplot(2, 2, 2)
    e = Ellipse((0, 0), 2 * c, 2 * b, facecolor="none", edgecolor="k")
    st.add_artist(e)
    if fold:
        for i in range(len(channel_names)):
            plt.plot(
                z[channels == channel_names[i]],
                abs(y)[channels == channel_names[i]],
                ".",
                color=channel_colors[i],
            )
        plt.xlim(-dlim * ay / az, dlim * ay / az)
        plt.ylim(0, dlim)
    else:
        for i in range(len(channel_names)):
            plt.plot(
                z[channels == channel_names[i]],
                y[channels == channel_names[i]],
                ".",
                color=channel_colors[i],
            )
        plt.xlim(-dlim * ay / az, dlim * ay / az)
        plt.ylim(-dlim, dlim)
    plt.xlabel("Z")
    plt.ylabel("Y")

    # XZ plot
    st = plt.subplot(2, 2, 3)
    e = Ellipse((0, 0), 2 * a, 2 * c, facecolor="none", edgecolor="k")
    st.add_artist(e)
    if fold:
        for i in range(len(channel_names)):
            plt.plot(
                abs(x)[channels == channel_names[i]],
                z[channels == channel_names[i]],
                ".",
                color=channel_colors[i],
            )
        plt.xlim(0, dlim)
    else:
        for i in range(len(channel_names)):
            plt.plot(
                x[channels == channel_names[i]],
                z[channels == channel_names[i]],
                ".",
                color=channel_colors[i],
            )
        plt.xlim(-dlim, dlim)
    plt.ylim(-dlim * ax / az, dlim * ax / az)
    plt.xlabel("X")
    plt.ylabel("Z")

    # Legend
    plt.subplot(2, 2, 4)
    plt.plot()
    if fold:
        plt.xlim(0, dlim)
        plt.ylim(0, dlim * ax / ay)
    else:
        plt.xlim(-dlim, dlim)
        plt.ylim(-dlim, dlim)
    plt.xlabel("")
    plt.ylabel("")
    channel_handles = [
        mlines.Line2D(
            [],
            [],
            color="white",
            marker="o",
            markerfacecolor=channel_colors[i],
            label=channel_names[i],
        )
        for i in range(len(channel_names))
    ]
    plt.legend(handles=channel_handles)

    # Figure title
    plt.suptitle(title)

    # Save output
    if type(None) != type(outpng):
        plt.savefig(outpng, format="png")
        plt.close("all")
        outpng.close()


def export(path, exp_format=None):
    """Export the current plot to the provided path.

    Note:
      The formats should be specified.
      Also, a plot object can be exported instead of the current one.

    Args:
      path (string): output file path.
      exp_format (string): output format (supported: pdf, png, jpg).
    """

    # SET PARAMS ===============================================================

    if exp_format is None:
        exp_format = "pdf"

    # EXPORT ===================================================================

    if exp_format in ["pdf", "png", "jpg"]:

        # Check presence of the extension
        path = pt.add_extension(path, "." + exp_format)

        # Save plot
        if exp_format == "pdf":
            pp = PdfPages(path)
            plt.savefig(pp, format=exp_format)
            pp.close()
        else:
            plt.savefig(path, format=exp_format)


def export_mask_png(outpath, im, title):
    """
    Export png of an image, max projected.

    Args:
    Returns:
    """

    # Create new plot figure
    fig = plt.figure()

    # Plot projection
    if len(im.shape) == 3:
        plt.imshow(im.max(0).astype("u4"))
    else:
        plt.imshow(im.astype("u4"))

    # Hide axes
    plt.gca().get_xaxis().set_visible(False)
    plt.gca().get_yaxis().set_visible(False)

    # Add title
    plt.title(title)
    set_font_size(8)

    # Export as png
    export(outpath, "png")

    # Close plot figure
    plt.close(fig)


def get_nsf_label(nsfi, seg_type=None, an_type=None):
    """Get the proper plot label for the selected nuclear selection feature,
    based on the current segmentation type.
    The function needs one of the two: an_type or seg_type.

    Args:
      nsfi (int): nuclear selection feature index according to pygpseq.const.
      seg_type (int): segmentation type according to pygpseq.const.
      an_type (int): analysis type according to pygpseq.const.

    Returns:
      string: the provided nuclear selection feature label.
    """

    if type(None) == type(an_type) and type(None) == type(seg_type):
        return ""

    # Check that the provided nuclear selection feature index exists
    if nsfi not in range(len(const.NSEL_FIELDS)):
        return None

    # Retrieve label
    label = const.NSEL_LABELS[nsfi]

    # Change label if required
    if label == "auto":
        if nsfi == const.NSEL_SIZE:
            if seg_type == const.SEG_3D or an_type == const.AN_3D:
                label = "Volume [vx]"
            else:
                label = "Area [px]"
        elif nsfi == const.NSEL_SHAPE:
            if seg_type == const.SEG_3D or an_type == const.AN_3D:
                label = "Sphericity"
            else:
                label = "Circularity"

    # Output
    return label


def multi_condition_boxplot(profiles, summaries, merged, an_type=None, **kwargs):
    """Generate multi condition boxplots.

    Args:
      profiles (list): list of profile dictionaries.
      summaries (list): list of condition nuclear summaries.
      merged (list): list of condition nuclear single-pixel tables.
      an_type (int): analysis type as defined in pygpseq.const.
      **kwargs: additional keyword arguments.

    Returns:
      plt.figure: figure canvas with plot.
    """

    # Create new plot figure
    fig = plt.figure(figsize=(5 + 3 * len(profiles), 12))

    # Calculate boxplot relative widths
    bp_widths = np.asarray([p["n"] for p in profiles], dtype="float")
    bp_widths = bp_widths / bp_widths.max() * 0.9

    # Set labels.
    size_ylab = get_nsf_label(const.NSEL_SIZE, an_type=an_type)
    shape_ylab = get_nsf_label(const.NSEL_SHAPE, an_type=an_type)

    # Plot size
    boxplot(
        fig,
        321,
        [s["size"] for s in summaries],
        [p["condition"] for p in profiles],
        bp_widths.tolist(),
        size_ylab + "  [per nucleus]",
        [p["n"] for p in profiles],
        xaxis1=False,
        forceSciNotation=True,
    )

    # Plot shape
    boxplot(
        fig,
        322,
        [s["shape"] for s in summaries],
        [p["condition"] for p in profiles],
        bp_widths.tolist(),
        shape_ylab + " [per nucleus]",
        [p["n"] for p in profiles],
        xaxis1=False,
        ylab_pos="right",
        forceSciNotation=True,
    )

    # Plot intensity average
    boxplot(
        fig,
        323,
        [s["meanI"] for s in summaries],
        [p["condition"] for p in profiles],
        bp_widths.tolist(),
        "mean(DNA [a.u.]) [per nucleus]",
        [p["n"] for p in profiles],
        xaxis1=False,
        xaxis2=False,
        forceSciNotation=True,
    )

    # Plot intensity sum
    boxplot(
        fig,
        324,
        [s["sumI"] for s in summaries],
        [p["condition"] for p in profiles],
        bp_widths.tolist(),
        "sum(DNA [a.u.]) [per nucleus]",
        [p["n"] for p in profiles],
        xaxis1=False,
        xaxis2=False,
        ylab_pos="right",
        forceSciNotation=True,
    )

    # Calculate boxplot relative widths
    bp_widths = np.asarray([m.shape[0] for m in merged], dtype="float")
    bp_widths = bp_widths / bp_widths.max() * 0.9

    # Plot single-pixel DNA intensity
    boxplot(
        fig,
        325,
        [m["dna"] for m in merged],
        [p["condition"] for p in profiles],
        bp_widths.tolist(),
        "DNA [a.u.] [per px]",
        [m.shape[0] for m in merged],
        xaxis2=False,
        forceSciNotation=True,
    )

    # Plot single-pixel Signal intensity
    boxplot(
        fig,
        326,
        [m["sig"] for m in merged],
        [p["condition"] for p in profiles],
        bp_widths.tolist(),
        "Signal [a.u.] [per px]",
        [m.shape[0] for m in merged],
        xaxis2=False,
        ylab_pos="right",
        forceSciNotation=True,
    )

    return fig


def multi_condition_profiles(
    profiles, font_size=None, hspace=None, wspace=None, yfield=None, **kwargs
):
    """Save three single condition profiles.

    Args:
      profiles (list): [{dtype:{x:float, y:float}, n:int, condition:string}].
      font_size (int)
      hspace (float): horizontal intra-plot spacing.
      wspace (float): vertical intra-plot spacing.
      yfield (string): y-axis label.
      **kwargs:
    """

    # SET PARAMS ===============================================================

    if type(None) == type(yfield):
        yfield = "mean"

    # PLOT =====================================================================

    fig = plt.figure(figsize=(12, 10))
    fig.subplots_adjust(top=0.9 - 0.01 * (len(profiles) / 2), bottom=0.1)

    # Prepare palette
    if len(profiles) <= 12:
        palette = sns.color_palette("Paired")
    else:
        palette = sns.color_palette("hls", len(profiles))

    # Setup subplots spacing
    for i in range(len(profiles)):
        profile = profiles[i]
        n_nuclei = profile["n"]

        single_condition_profiles(
            profile,
            color=palette[i],
            n_nuclei=n_nuclei,
            new_figure=False,
            yfield=yfield,
            **kwargs
        )

    # Add legend
    plt.subplot(3, 2, 1)
    plt.legend(
        bbox_to_anchor=(0.0, 1.12, 1.0, 0.102),
        loc=3,
        ncol=2,
        mode="expand",
        borderaxespad=0.0,
        fontsize="xx-small",
    )

    # Add title
    title = "GPSeq profiles [" + yfield + "]"
    if "title_comment" in kwargs.keys():
        title += " [%s]" % (kwargs["title_comment"],)
    if "sigma_smooth" in kwargs.keys():
        title += " [sigma = %.2f]" % (kwargs["sigma_smooth"],)
    if "nbins" in kwargs.keys():
        title += " [nbins = %d]" % (kwargs["nbins"],)
    plt.suptitle(title)

    # Output
    return fig


def multi_condition_single_boxplot(
    profiles, data, data_field, bp_widths, ylab, out_png, **kwargs
):
    """Generate multi condition boxplot.
    The produced plot is directly exported in PNG

    Args:
      profiles (list): list of profile dictionaries.
      data (list): list of datasets to boxplot.
      data_field (string): column (in the datasets) to boxplot.
      bp_widths (list(float)): boxplot widths.
      ylab (string): y-axis label
      out_png (string): path to png output folder.
      **kwargs: optional keyword arguments (suffix).
    """

    # Plot
    fig = plt.figure(figsize=(4 + len(profiles), 4))
    boxplot(
        fig,
        111,
        [d[data_field] for d in data],
        [i + 1 for i in range(len(profiles))],
        bp_widths.tolist(),
        ylab,
        [p["n"] for p in profiles],
        forceSciNotation=True,
    )

    # Save plot
    fname = out_png + "boxplot." + data_field
    if "suffix" in kwargs.keys():
        fname += kwargs["suffix"]
    fname += ".png"
    export(fname, "png")
    plt.close(fig)


def ortho_3d(coords, scale=None, dot_coords=None, c=None, aspect=None, channels=None):
    """
    Plot orthogonal projections with 3 major eigenvectors.

    Args:
        coords (nd.array): coordinate table with x, y, z columns.
        scale (int): size of plotted eigenvectors.
    """

    # Set defaults -------------------------------------------------------------
    if type(None) == type(scale):
        scale = 100
    if type(None) == type(aspect):
        aspect = (1, 1, 1)
    az, ay, ax = aspect
    if type(None) == type(channels):
        channels = np.repeat("b", coords.shape[1])

    # Plot ---------------------------------------------------------------------
    plt.figure(figsize=(10, 10))

    # Extract coordinates
    x, y, z = coords
    xv, yv, zv = stt.extract_3ev(coords)

    # Prepare colors
    cmap = plt.get_cmap("rainbow")
    channel_names = list(set(channels))
    channel_names.sort()
    channel_colors = cmap(np.linspace(0, 1, len(channel_names)))

    # XY plot
    plt.subplot(2, 2, 1)
    plt.plot(
        [xv[0] * -scale * 3, xv[0] * scale * 3],
        [yv[0] * -scale * 3, yv[0] * scale * 3],
        color="r",
    )
    plt.plot(
        [xv[1] * -scale * 2, xv[1] * scale * 2],
        [yv[1] * -scale * 2, yv[1] * scale * 2],
        color="b",
    )
    plt.plot(
        [xv[2] * -scale, xv[2] * scale], [yv[2] * -scale, yv[2] * scale], color="g"
    )
    plt.plot(x, y, "k.", alpha=0.05)
    plt.xlim(coords.min() - 10, coords.max() + 10)
    plt.ylim((coords.min() - 10) * (ax / ay), (coords.max() + 10) * (ax / ay))
    plt.xlabel("X")
    plt.ylabel("Y")

    if type(None) != type(dot_coords):
        for i in range(len(channel_names)):
            plt.plot(
                dot_coords[0][channels == channel_names[i]],
                dot_coords[1][channels == channel_names[i]],
                ".",
                color=channel_colors[i],
            )
    if type(None) != type(c):
        plt.axvline(c, color="g", linestyle=":")
        plt.axvline(-c, color="g", linestyle=":")

    # ZY plot
    plt.subplot(2, 2, 2)
    plt.plot(
        [zv[0] * -scale * 3, zv[0] * scale * 3],
        [yv[0] * -scale * 3, yv[0] * scale * 3],
        color="r",
    )
    plt.plot(
        [zv[1] * -scale * 2, zv[1] * scale * 2],
        [yv[1] * -scale * 2, yv[1] * scale * 2],
        color="b",
    )
    plt.plot(
        [zv[2] * -scale, zv[2] * scale], [yv[2] * -scale, yv[2] * scale], color="g"
    )
    plt.plot(z, y, "k.", alpha=0.05)
    plt.xlim((coords.min() - 10) * (ay / az), (coords.max() + 10) * (ay / az))
    plt.ylim(coords.min() - 10, coords.max() + 10)
    plt.xlabel("Z")
    plt.ylabel("Y")

    if type(None) != type(dot_coords):
        for i in range(len(channel_names)):
            plt.plot(
                dot_coords[2][channels == channel_names[i]],
                dot_coords[1][channels == channel_names[i]],
                ".",
                color=channel_colors[i],
            )
    if type(None) != type(c):
        plt.axvline(0, color="b", linestyle=":")

    # XZ plot
    plt.subplot(2, 2, 3)
    plt.plot(
        [xv[0] * -scale * 3, xv[0] * scale * 3],
        [zv[0] * -scale * 3, zv[0] * scale * 3],
        color="r",
    )
    plt.plot(
        [xv[1] * -scale * 2, xv[1] * scale * 2],
        [zv[1] * -scale * 2, zv[1] * scale * 2],
        color="b",
    )
    plt.plot(
        [xv[2] * -scale, xv[2] * scale], [zv[2] * -scale, zv[2] * scale], color="g"
    )
    plt.plot(x, z, "k.", alpha=0.05)
    plt.xlim(coords.min() - 10, coords.max() + 10)
    plt.ylim((coords.min() - 10) * (ax / az), (coords.max() + 10) * (ax / az))
    plt.xlabel("X")
    plt.ylabel("Z")

    if type(None) != type(dot_coords):
        for i in range(len(channel_names)):
            plt.plot(
                dot_coords[0][channels == channel_names[i]],
                dot_coords[2][channels == channel_names[i]],
                ".",
                color=channel_colors[i],
            )
    if type(None) != type(c):
        plt.axhline(0, color="b", linestyle=":")
        plt.axvline(c, color="g", linestyle=":")
        plt.axvline(-c, color="g", linestyle=":")

    # Legend
    plt.subplot(2, 2, 4)
    plt.plot()
    plt.xlim(coords.min() - 10, coords.max() + 10)
    plt.ylim((coords.min() - 10) * (ax / az), (coords.max() + 10) * (ax / az))
    plt.xlabel("")
    plt.ylabel("")
    channel_handles = [
        mlines.Line2D(
            [],
            [],
            color="white",
            marker="o",
            markerfacecolor=channel_colors[i],
            label=channel_names[i],
        )
        for i in range(len(channel_names))
    ]
    plt.legend(handles=channel_handles)


def profile_density_scatterplot(
    pp, field_label, nbins, xyrange, profile, new_figure=None, dlabel=None
):
    """Plot profile pixel density scatterplot.

    Args:
      pp (np.array): bin-normalized density.
      field_label (string): y-axis label.
      nbins (int): number of bins.
      xyrange (list): x/y-axis range [min, max].
      profile (list(np.array)): mean, median, mode, std (smoothed and raw).
      new_figure (boolean): whether to plot in a new figure (and return it).
      dlabel (string): distance label.

    Returns:
      plt.figure or None: if new_figure is True, return canvas object.
    """

    # CHECK PARAMS =============================================================

    if new_figure is None:
        new_figure = True

    if new_figure:
        fig = plt.figure(figsize=[10, 10])

    # PLOT =====================================================================

    plt.imshow(pp, origin="lower", interpolation="none", cmap="gray")

    cbar = plt.colorbar()
    cbar.set_label("Bin-relative local density", rotation=270, labelpad=20)

    plt.xlabel(dlabel)
    plt.ylabel(field_label)

    plt.hold(True)

    plt.plot(
        range(nbins),
        nbins * (profile["median"] - xyrange[1][0]) / (xyrange[1][1] - xyrange[1][0]),
        "g",
        linewidth=2,
    )
    plt.plot(
        range(nbins),
        nbins * (profile["mean"] - xyrange[1][0]) / (xyrange[1][1] - xyrange[1][0]),
        "y",
        linewidth=2,
    )
    plt.plot(
        range(nbins),
        nbins * (profile["mode"] - xyrange[1][0]) / (xyrange[1][1] - xyrange[1][0]),
        "m",
        linewidth=2,
    )
    plt.plot(
        range(nbins),
        nbins * (profile["max"] - xyrange[1][0]) / (xyrange[1][1] - xyrange[1][0]),
        "b",
        linewidth=2,
    )

    plt.xlim([0, nbins - 1])
    plt.ylim([0, nbins - 1])

    # Adjust ticks
    ax = plt.gca()

    ticks = ax.get_xticks() / float(nbins) * xyrange[0][1]
    ticks = [const.SCI_FORMAT % n for n in ticks]
    ax.set_xticklabels(ticks)

    ticks = ax.get_yticks() / float(nbins) * xyrange[1][1]
    ticks = [const.SCI_FORMAT % n for n in ticks]
    ax.set_yticklabels(ticks)

    if new_figure:
        return fig


def save_tif(
    path,
    img,
    dtype,
    compressed,
    bundled_axes="CZYX",
    inMicrons=False,
    ResolutionZ=None,
    forImageJ=False,
    **kwargs
):
    # Add channel axis for ImageJ compatibility
    if "C" not in bundled_axes:
        bundled_axes = "C" + bundled_axes

    # Add missing axes
    while len(bundled_axes) > len(img.shape):
        new_shape = [1]
        [new_shape.append(n) for n in img.shape]
        img.shape = new_shape

    assert_msg = "shape mismatch between bundled axes and image."
    assert len(bundled_axes) == len(img.shape), assert_msg

    metadata = {"axes": bundled_axes}
    if inMicrons:
        metadata["unit"] = "um"
    if type(None) != ResolutionZ:
        metadata["spacing"] = ResolutionZ

    if compressed:
        tifffile.imsave(
            path,
            img.astype(dtype),
            shape=img.shape,
            compress=9,
            dtype=dtype,
            imagej=forImageJ,
            metadata=metadata,
            **kwargs
        )
    else:
        tifffile.imsave(
            path,
            img.astype(dtype),
            shape=img.shape,
            compress=0,
            dtype=dtype,
            imagej=forImageJ,
            metadata=metadata,
            **kwargs
        )


def set_font_size(size):
    """Set plot font size."""
    rc("font", **{"size": size})


def single_condition_profiles(
    profiles,
    font_size=None,
    color=None,
    hspace=None,
    wspace=None,
    n_nuclei=None,
    new_figure=None,
    cdescr=None,
    yfield=None,
    dlabel=None,
    **kwargs
):
    """Save three single condition profiles.

    Args:
      profiles (dict): {dtype:{x:float, y:float}, n: int, condition:string}.
      font_size (int): (opt).
      color (string): (opt).
      hspace (float): horizontal intra-plot spacing (opt).
      wspace (float): vertical intra-plot spacing (opt).
      n_nuclei (int): number of nuclei (opt).
      new_figure (bool): plot line color (opt).
      cdescr (dict): better condition description (opt).
      yfield (string): which measure to plot (opt).
      dlabel (string): distance label (opt).
      **kwargs: title_comment, sigma, nbins.

    Returns:
      plt.figure or None: figure canvase if new_figure is True.
    """

    # SET PARAMS ===============================================================

    if type(None) != type(cdescr) and profiles["condition"] in cdescr.keys():
        profiles["condition"] = cdescr[profiles["condition"]]

    if type(None) == type(dlabel):
        dlabel = "x"

    if type(None) == type(font_size):
        font_size = 12

    if type(None) == type(hspace):
        hspace = 0.4

    if type(None) == type(new_figure):
        new_figure = True

    if type(None) == type(n_nuclei):
        n_nuclei = "~unknown~"

    if type(None) == type(wspace):
        wspace = 0.4

    if type(None) == type(yfield):
        yfield = "mean"

    # Plot options
    opt = {}
    if type(None) != type(color):
        opt["c"] = color

    # PLOT =====================================================================

    # Init canvas --------------------------------------------------------------
    if new_figure:
        fig = plt.figure(figsize=(12, 8))

        # Main title
        title = 'GPSeq profiles for condition "%s" (n.nuclei = %d)' % (
            profiles["condition"],
            n_nuclei,
        )

        # Add user-defined comment
        if "title_comment" in kwargs.keys():
            title += " [%s]" % (kwargs["title_comment"],)

        # Additional information
        if "sigma_smooth" in kwargs.keys():
            title += " [sigma = %.3f]" % (kwargs["sigma_smooth"],)
        if "nbins" in kwargs.keys():
            title += " [nbins = %d]" % (kwargs["nbins"],)

        # Add suptitle
        plt.suptitle(title, fontsize=font_size)

    # Setup subplots spacing
    plt.subplots_adjust(wspace=wspace, hspace=hspace)

    # Number of steps (precision)
    nsteps = len(profiles["ratio"]["x"])

    # DNA profile --------------------------------------------------------------
    plt.subplot(3, 2, 1)
    plt.hold(True)
    plt.plot(
        profiles["dna"]["x"],
        profiles["dna"][yfield],
        label=profiles["condition"],
        **opt
    )
    plt.xlabel(dlabel)
    plt.ylabel("DNA [a.u.]")
    plt.ticklabel_format(style="sci", axis="y", scilimits=(0, 0))
    set_font_size(font_size)

    # DNA first derivative -----------------------------------------------------
    plt.subplot(3, 2, 2)
    plt.hold(True)
    plt.plot(
        profiles["dna"]["x"][1:nsteps],
        np.diff(profiles["dna"][yfield]),
        label=profiles["condition"],
        **opt
    )
    plt.axhline(0, color="k", linestyle="-.")
    plt.xlabel(dlabel)
    plt.ylabel("diff(DNA [a.u.])")
    plt.ticklabel_format(style="sci", axis="y", scilimits=(0, 0))
    set_font_size(font_size)

    # Signal profile -----------------------------------------------------------
    plt.subplot(3, 2, 3)
    plt.hold(True)
    plt.plot(
        profiles["sig"]["x"],
        profiles["sig"][yfield],
        label=profiles["condition"],
        **opt
    )
    plt.xlabel(dlabel)
    plt.ylabel("Signal [a.u.]")
    plt.ticklabel_format(style="sci", axis="y", scilimits=(0, 0))
    set_font_size(font_size)

    # Signal first derivative --------------------------------------------------
    plt.subplot(3, 2, 4)
    plt.hold(True)
    plt.plot(
        profiles["sig"]["x"][1:nsteps],
        np.diff(profiles["sig"][yfield]),
        label=profiles["condition"],
        **opt
    )
    plt.axhline(0, color="k", linestyle="-.")
    plt.xlabel(dlabel)
    plt.ylabel("diff(Signal [a.u.])")
    plt.ticklabel_format(style="sci", axis="y", scilimits=(0, 0))
    set_font_size(font_size)

    # Ratio profile ------------------------------------------------------------
    plt.subplot(3, 2, 5)
    plt.hold(True)
    plt.plot(
        profiles["ratio"]["x"],
        profiles["ratio"][yfield],
        label=profiles["condition"],
        **opt
    )
    plt.xlabel(dlabel)
    plt.ylabel("Signal/DNA")
    plt.ticklabel_format(style="sci", axis="y", scilimits=(0, 0))
    set_font_size(font_size)

    # Ratio first derivative ---------------------------------------------------
    plt.subplot(3, 2, 6)
    plt.hold(True)
    plt.plot(
        profiles["ratio"]["x"][1:nsteps],
        np.diff(profiles["ratio"][yfield]),
        label=profiles["condition"],
        **opt
    )
    plt.axhline(0, color="k", linestyle="-.")
    plt.xlabel(dlabel)
    plt.ylabel("dif(Signal/DNA)")
    plt.ticklabel_format(style="sci", axis="y", scilimits=(0, 0))
    set_font_size(font_size)

    # Close figure -------------------------------------------------------------
    if new_figure:
        return fig


def single_pixel_study(xs, ys, field_label, profile, nbins=None, **kwargs):
    """Plot a study of single-pixel and general behaviour of profile data.

    Args:
      xs (np.array): pixel (relative) distance from nuclear lamina.
      ys (np.array): pixel channel intensity.
      field_label (string): y-axis label.
      profile (dict): profile to study.
      nbins (int): study precision (opt, def 200).

    Returns:
      plt.figure: current figure canvas.
    """

    # CHECK PARAMS =============================================================

    # Default values
    if nbins is None:
        nbins = 200

    # PREPARE DATA =============================================================

    # XY range
    xyrange = [[xs.min(), xs.max()], [ys.min(), ys.max()]]

    # Bin data
    assigned_bins = np.digitize(xs, np.linspace(xs.min(), xs.max(), nbins))
    binned_ydat = [[] for i in range(nbins)]
    for bin_id in range(nbins):
        binned_ydat[bin_id] = ys[np.where(assigned_bins == bin_id)]

    # Calculate local density
    hh, locx, locy = scipy.histogram2d(xs, ys, range=xyrange, bins=[nbins, nbins])
    posx = np.digitize(xs, locx)
    posy = np.digitize(ys, locy)

    # Hide empty bins
    pp = hh.T / hh.T.max(0)
    pp[:, hh.T.sum(0) == 0] = 0

    # MULTIPLOT ================================================================

    # Init figure --------------------------------------------------------------
    fig = plt.figure(figsize=[40, 15])
    rc("font", **{"size": 15})

    # Number of pixels plot ----------------------------------------------------
    plt.subplot2grid((4, 3), (0, 0))
    plt.plot(range(nbins), hh.sum(1))
    plt.ylabel("Number of pixels")
    plt.xlabel(kwargs["dlabel"])

    # Adjust ticks
    ax = plt.gca()
    ticks = ax.get_xticks() / float(nbins) * xyrange[0][1]
    ticks = [const.SCI_FORMAT % n for n in ticks]
    ax.get_xaxis().set_visible(False)
    plt.ticklabel_format(style="sci", axis="y", scilimits=(0, 0))

    # Boxplot ------------------------------------------------------------------
    ax = plt.subplot2grid((4, 3), (1, 0))
    plt.boxplot(binned_ydat, sym="k.")
    plt.ylabel(field_label)
    ax.get_xaxis().set_visible(False)
    plt.ticklabel_format(style="sci", axis="y", scilimits=(0, 0))

    # STD plot -----------------------------------------------------------------
    plt.subplot2grid((4, 3), (2, 0))
    plt.plot(profile["std_raw"], "r")
    plt.plot(profile["std_raw"], "r.")
    plt.ylabel("std(" + field_label + ")")
    plt.hold(True)
    plt.plot(profile["std"], "k")
    plt.xlabel(kwargs["dlabel"])

    # Adjust ticks
    ax = plt.gca()
    ax.set_xticklabels(ticks)
    plt.ticklabel_format(style="sci", axis="y", scilimits=(0, 0))

    # Mean plot ----------------------------------------------------------------
    plt.subplot2grid((4, 3), (0, 1))
    plt.plot(profile["mean_raw"], "b")
    plt.plot(profile["mean_raw"], "b.")
    plt.ylabel("mean(" + field_label + ")")
    plt.hold(True)
    plt.plot(profile["mean"], "y")

    # Adjust ticks
    ax = plt.gca()
    ax.get_xaxis().set_visible(False)
    plt.ticklabel_format(style="sci", axis="y", scilimits=(0, 0))

    # Median plot --------------------------------------------------------------
    plt.subplot2grid((4, 3), (1, 1))
    plt.plot(profile["median_raw"], "r")
    plt.plot(profile["median_raw"], "r.")
    plt.ylabel("median(" + field_label + ")")
    plt.hold(True)
    plt.plot(profile["median"], "g")

    # Adjust ticks
    ax = plt.gca()
    ax.get_xaxis().set_visible(False)
    plt.ticklabel_format(style="sci", axis="y", scilimits=(0, 0))

    # Mode plot ----------------------------------------------------------------
    plt.subplot2grid((4, 3), (2, 1))
    plt.plot(profile["mode_raw"], "c")
    plt.plot(profile["mode_raw"], "c.")
    plt.ylabel("mode(" + field_label + ")")
    plt.hold(True)
    plt.plot(profile["mode"], "m")
    plt.xlabel(kwargs["dlabel"])

    # Adjust ticks
    ax = plt.gca()
    ax.get_xaxis().set_visible(False)
    plt.ticklabel_format(style="sci", axis="y", scilimits=(0, 0))

    # Max plot -----------------------------------------------------------------
    plt.subplot2grid((4, 3), (3, 1))
    plt.plot(profile["max_raw"], "c")
    plt.plot(profile["max_raw"], "c.")
    plt.ylabel("max(" + field_label + ")")
    plt.hold(True)
    plt.plot(profile["max"], "b")
    plt.xlabel(kwargs["dlabel"])

    # Adjust ticks
    ax = plt.gca()
    ax.set_xticklabels(ticks)
    plt.ticklabel_format(style="sci", axis="y", scilimits=(0, 0))

    # Local density scatterplot ------------------------------------------------
    plt.subplot2grid((4, 3), (0, 2), rowspan=4, colspan=2)
    profile_density_scatterplot(
        pp,
        field_label,
        nbins,
        xyrange,
        profile,
        dlabel=kwargs["dlabel"],
        new_figure=False,
    )

    # SINGLE PLOTS =============================================================

    if kwargs["plotting"]:

        # Density scatterplot --------------------------------------------------

        # Retrieve tmp figure
        tmp_fig = profile_density_scatterplot(
            pp, field_label, nbins, xyrange, profile, dlabel=kwargs["dlabel"]
        )

        # Get filename suffix
        suffix = field_label.lower().split(" ")[0].replace("/", "_")

        # Get condition description
        if kwargs["cond_name"] in kwargs["cdescr"].keys():
            suptitle = kwargs["cdescr"][kwargs["cond_name"]]
        else:
            suptitle = kwargs["cond_name"]
        plt.suptitle(suptitle)

        # Set output name
        fname = kwargs["outdir"] + const.OUTDIR_PNG_REPORT + kwargs["cond_name"]
        fname += "." + suffix + ".density_profile"

        # Add partial suffix
        if "partial" in kwargs.keys() and kwargs["partial"]:
            fname += ".part"

        # Add general suffix
        fname += kwargs["suffix"] + ".png"

        # Export
        export(fname, "png")

        # Close tmp figure
        plt.close(tmp_fig)

        # Number of pixels plot  -----------------------------------------------

        # Init tmp figure
        tmp_fig = plt.figure()

        # Plot
        plt.plot(range(nbins), hh.sum(1))
        plt.ylabel("Number of pixels")
        plt.xlabel(kwargs["dlabel"])
        ticks = ax.get_xticks() / float(nbins) * xyrange[0][1]
        ticks = [const.SCI_FORMAT % n for n in ticks]
        plt.gca().set_xticklabels(ticks)
        plt.ticklabel_format(style="sci", axis="y", scilimits=(0, 0))

        # Set suptitle
        if "partial" in kwargs.keys() and kwargs["partial"]:
            suptitle += " [partial volume "
            suptitle += str(kwargs["part_n_erosion"]) + "]"
        plt.suptitle(suptitle)

        # Set output name
        fname = kwargs["outdir"] + const.OUTDIR_PNG_REPORT + kwargs["cond_name"]
        fname += "." + suffix + ".npx"

        # Add partial suffix
        if "partial" in kwargs.keys() and kwargs["partial"]:
            fname += ".part"

        # Add general suffix
        fname += kwargs["suffix"] + ".png"

        # Export
        export(fname, "png")

        # Close tmp figure
        plt.close(tmp_fig)

    # Return MULTIPLOT canvas
    return fig


# END ==========================================================================

################################################################################

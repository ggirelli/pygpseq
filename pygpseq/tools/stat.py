# -*- coding: utf-8 -*-

"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
@description: statistic operations library.
"""

# DEPENDENCIES =================================================================

import math
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import convolve
from scipy import stats

from pygpseq import const

from pygpseq.tools import vector as vt

# FUNCTIONS ====================================================================


def angle_between_points(p0, c, p1):
    """
    c is the center point; result is in degrees
    From http://phrogz.net/angle-between-three-points
    """
    p0 = np.array(p0)
    c = np.array(c)
    p1 = np.array(p1)

    p0c = np.sqrt(np.sum((p0 - c) ** 2))
    p1c = np.sqrt(np.sum((p1 - c) ** 2))
    p01 = np.sqrt(np.sum((p0 - p1) ** 2))

    d = round((p0c ** 2 + p1c ** 2 - p01 ** 2), 6)
    n = round((2 * p0c * p1c), 6)

    try:
        tetha = math.acos(d / n)
    except ValueError as e:
        print("Something went wrong when calculating an angle...")
        raise

    return tetha / math.pi * 180


def binned_mode(x, nbins):
    """Identify binned mode.

    Args:
      x (np.array): dataset.
      nbins (int): number of bins.

    Returns:
      int: the most occupied bin in the provided dataset.
    """

    if len(x) == 0:
        return np.nan

    # Bin breaks
    breaks = np.linspace(0, max(x), nbins)

    # Assign data to the bins
    assigned_bins = np.digitize(x, breaks)

    # Count bins occurrences
    occ = vt.uniquec(assigned_bins)
    counts = np.array([e[1] for e in occ])

    # Order counts
    ordered = np.argsort(counts).tolist()
    ordered.reverse()
    occ = [occ[i] for i in ordered]

    # Return mode
    return breaks[occ[0][0] - 1]


def binned_profile(x, y, nbins=None):
    """Produce an approximation of sparse data by binning it.

    Args:
      x (numeric): x coordinates.
      y (numeric): y coordinates.
      nbins (int): curve precision (opt, def: 200).

    Returns:
      np.array: profiles.
    """

    if nbins is None:
        nbins = 200

    # Check format
    y = np.array(y)

    # Bin breaks
    breaks = np.linspace(0, max(x), nbins)

    # Assign data to the bins
    assigned_bins = np.digitize(x, breaks)

    # Get mean and median for every bin
    data = np.zeros(
        (len(breaks),),
        dtype=[
            ("breaks", "f"),
            ("mean", "f"),
            ("median", "f"),
            ("std", "f"),
            ("mode", "f"),
            ("max", "f"),
            ("mean_raw", "f"),
            ("median_raw", "f"),
            ("std_raw", "f"),
            ("mode_raw", "f"),
            ("max_raw", "f"),
            ("n", "f"),
        ],
    )
    for bin_id in range(assigned_bins.max()):
        where = np.where(assigned_bins == bin_id)
        data["breaks"][bin_id] = breaks[bin_id]
        if where[0].shape[0] != 0:
            data["mean"][bin_id] = np.mean(y[where])
            data["median"][bin_id] = np.median(y[where])
            data["mode"][bin_id] = binned_mode(y[where], nbins)
            data["std"][bin_id] = np.std(y[where])
            data["max"][bin_id] = np.max(y[where])
            data["n"][bin_id] = len(y[where])
        else:
            data["mean"][bin_id] = np.nan
            data["median"][bin_id] = np.nan
            data["mode"][bin_id] = np.nan
            data["std"][bin_id] = np.nan
            data["max"][bin_id] = np.nan
            data["n"][bin_id] = 0

    # Output
    return data


def calc_density(data, **kwargs):
    """
    Calculate the Gaussian KDE of the provided data series.

    Args:
      sigma_density (float): standard deviation used for covariance calculation.
      nbins (int): #steps for the density curve calculation (opt, def: 1000).

    Returns:
      dict: density curve (x, y) and function (f).
    """

    # Default values
    if "sigma_density" not in kwargs.keys():
        sigma_density = 0.1
    else:
        sigma_density = kwargs["sigma_density"]

    nbins = 1000 if "nbins" not in kwargs.keys() else kwargs["nbins"]
    # If only one nucleus was found
    if len(data) == 1:
        f = eval("lambda x: 1 if x == " + str(data[0]) + " else 0")
        f = np.vectorize(f)
        return {"x": np.array([data[0]]), "y": np.array([1]), "f": f}

    # Prepare density function
    density = stats.gaussian_kde(data)

    # Re-compute covariance
    density.covariance_factor = lambda: sigma_density
    density._compute_covariance()

    # Output
    out = {"x": np.linspace(min(data), max(data), nbins)}
    out["f"] = density
    out["y"] = density(out["x"])
    return out


def calc_theta(a, b):
    """
    Calculate rotation angle based on a (opposite) and b (adjacent) sides.

    Return:
        float: theta in rad.
    """
    c = np.sqrt(a ** 2 + b ** 2)
    if a > 0 and b < 0 or (a >= 0 or b <= 0) and a < 0 and b < 0:
        return np.arccos(a / c)
    else:
        return -np.arccos(a / c)


def centered_coords_3d(img, dot_coords=None):
    """
    Extract coordinates from binary image and center them on the origin.

    Args:
        img (nd.array): binary image.
    """
    z, x, y = np.nonzero(img)

    if type(None) != type(dot_coords):
        zd, xd, yd = dot_coords
        xd = xd - np.mean(x)
        yd = yd - np.mean(y)
        zd = zd - np.mean(z)

    x = x - np.mean(x)
    y = y - np.mean(y)
    z = z - np.mean(z)

    if type(None) != type(dot_coords):
        return (x, y, z, xd, yd, zd)
    else:
        return (x, y, z, None, None, None)


def extract_3ev(coords):
    """
    Extract 3 major eigen vectors.

    Args:
        coords (nd.array): coordinates table with one point per row.

    Returns:
        tuple: major 3 eigen vectors. Coordinate order matches input columns.
    """

    cov = np.cov(coords)
    evals, evecs = np.linalg.eig(cov)

    sort_indices = np.argsort(evals)[::-1]
    a_v1, b_v1, c_v1 = evecs[:, sort_indices[0]]
    a_v2, b_v2, c_v2 = evecs[:, sort_indices[1]]
    a_v3, b_v3, c_v3 = evecs[:, sort_indices[2]]

    av = [a_v1, a_v2, a_v3]
    bv = [b_v1, b_v2, b_v3]
    cv = [c_v1, c_v2, c_v3]

    return (av, bv, cv)


def get_fwhm(xs, ys):
    """Calculate FWHM of highest peak in a curve.

    Args:
      xs (np.array): x-coordinates of the curve.
      ys (np.array): y-coordinates of the curve.

    Returns:
      list: FWHM interval x-coords of highest peak.
    """

    # CHECK PARAMS =============================================================

    # Must have the same length
    if len(xs) != len(ys):
        return None

    if len(ys) == 1:
        return [xs[0] - 1, xs[0] + 1]

    # GET FWHM =================================================================

    # Identify highest peak
    xmaxi = ys.tolist().index(max(ys))

    # Get FWHM range [left] ----------------------------------------------------

    if xmaxi != 0:
        # Get absolute difference to HM value
        x1 = abs(ys[range(xmaxi)] - max(ys) / 2)

        # Get threshold based on average distance of consecutive points
        thr = x1[0] if len(x1) == 1 else np.max(abs(np.diff(x1)))
        # Select values close to the HM (based on threshold and abs difference)
        selected = [i for i in range(len(x1)) if x1[i] <= thr]

        # Select left boundary
        if not selected:
            x1 = xs[range(xmaxi)][x1.tolist().index(min(x1))]
        else:
            x1 = xs[range(xmaxi)][max(selected)]
    else:
        x1 = min(xs)

    # Get FWHM range [right] ---------------------------------------------------

    if len(xs) != (xmaxi + 1):
        # Get absolute difference to HM value
        x2 = abs(ys[range(xmaxi + 1, len(ys))] - max(ys) / 2)

        # Get threshold based on average distance of consecutive points
        thr = x2[0] if len(x2) == 1 else np.max(abs(np.diff(x2)))
        # Select values close to the HM (based on threshold and abs difference)
        selected = [i for i in range(len(x2)) if x2[i] <= thr]

        # Select right boundary
        if not selected:
            x2 = xs[range(xmaxi + 1, len(xs))][x2.tolist().index(min(x2))]
        else:
            x2 = xs[range(xmaxi + 1, len(xs))][min(selected)]
    else:
        x2 = max(xs)

    # Output
    return [x1, x2]


def get_intercepts(ys, xs=None):
    """Given a series of y coordinates, identify the x-axis intercept.

    Args:
      ys (numeric): y coordinates series.
      xs (numeric): x coordinates series (optional).

    Returns:
      int: index of x-axis intercepting y (if xs == None).
      numeric: x coordinates interpolated point of intersection with x-axis.
    """

    # Return no intercept if no data was provided
    if len(ys) == 0:
        return []

    # Will contain the intercepts (pairs of) indexes
    out = []

    # Reformat coordinate series
    ys = np.array(ys)

    # Count ys
    nys = ys.shape[0]

    # Checking matching ys/xs size
    if type(None) != type(xs):
        xs = np.array(xs)
        if ys.shape[0] != xs.shape[0]:
            print(ys.shape[0])
            print(xs.shape[0])
            print("Discarded provided X coordinates.")
            xs = None

    # Perfect intersections ----------------------------------------------------

    # Identify zero-holding cells
    bys = ys == 0
    if bys.sum() != 0:
        # Save zero-holding cells index
        out.extend([i for i in range(nys) if bys[i]])

    # Interpolate intersections ------------------------------------------------

    # Identify positive/negative ys (convert to boolean)
    bys = np.zeros(ys.shape)
    bys[ys == 0] = np.nan
    bys[ys > 0] = 1
    bys[ys < 0] = -1

    # Identify intercepts by summing previous boolean cell value
    bys = np.array([bys[i] + bys[i + 1] for i in range(bys.shape[0] - 1)])
    bysi = [i for i in range(len(bys)) if (bys == 0)[i]]

    if type(None) == type(xs):
        # Interpolate index of intercepts
        out.extend([i + 0.5 for i in bysi])
    else:
        # Interpolate x of intercepts
        out = [xs[i] for i in out]
        out.extend([(xs[i] + xs[i + 1]) / 2.0 for i in bysi])

    # Output
    return out


def get_norm_pdf(mu, sigma, x):
    """Normal distribution N(mu, sigma) probability value at x.

    Args:
      mu (float): normal distribution mean.
      sigma (float): normal distribution standard deviation.
      x (float): x-coordinate.

    Returns:
      float: N(mu, sigma) probability value at x.
    """

    return (
        1
        / np.sqrt(2 * sigma ** 2 * np.pi)
        * np.exp(-((x - mu) ** 2) / (2 * sigma ** 2))
    )


def get_outliers(x, non=None, fig=None, close=None):
    """Identifies the outliers in a data set.

    Args:
      x (np.array): data set.
      non (bool): whether to return outliers or non-outliers.
      fig (plt.figure): for boxplot purposes.
      close (bool): whether to close the figure before return.

    Returns:
      list: (non-)outlier indexes.
    """

    if non is None:
        non = False
    if fig is None:
        fig = plt.figure()
        ax = fig.gca()
    if close is None:
        close = False

    # If no dataset is provided
    if len(x) == 0:
        return []

    # Identify outliers through boxplot
    bp = ax.boxplot(x)

    # Close figure
    if close:
        plt.close(tmp_fig)

    # Retrieve outliers
    outliers = []
    outliers.extend(bp["fliers"][0].get_data()[0].tolist())
    outliers.extend(bp["fliers"][0].get_data()[1].tolist())

    # Unique outliers
    outliers = set(outliers)

    # Output
    if not non:
        return [i for i in range(len(x)) if x[i] in outliers]
    else:
        return [i for i in range(len(x)) if x[i] not in outliers]


def gpartial(V, d, sigma):
    """Calculate the partial derivative of V along dimension d using a filter
    of size sigma.

    Based on code by Erik Wernersson, PhD.
    """

    w = round(8 * sigma + 2)
    if w % 2 == 0:
        w = w + 1
    w = 2 * w + 1
    if sigma == 1:
        w = 11

    if sigma == 0:
        dg = [0, -1, 1]
        g = [0, 0.5, 0.5]
    else:
        g = stats.norm.pdf(np.linspace(-w / 2.0, w / 2.0, w + 1), scale=sigma)
        x = np.linspace(-(w - 1) / 2, (w - 1) / 2, w + 1)
        k0 = 1 / np.sqrt(2 * np.pi * sigma ** 2.0)
        k1 = 1 / (2 * sigma ** 2)
        dg = -2 * k0 * k1 * x * np.exp(-k1 * x ** 2.0)

    # Iterate through slices
    intlist2 = []

    if len(V.shape) == 3:
        if d == 1:
            V = convolve(V, dg.reshape([1, 1, w + 1]), "same")
        else:
            V = convolve(V, g.reshape([1, 1, w + 1]), "same")

        if d == 2:
            V = convolve(V, dg.reshape([1, w + 1, 1]), "same")
        else:
            V = convolve(V, g.reshape([1, w + 1, 1]), "same")

        if d == 3:
            V = convolve(V, dg.reshape([w + 1, 1, 1]), "same")
        else:
            V = convolve(V, g.reshape([w + 1, 1, 1]), "same")

    elif len(V.shape) == 2:
        if d == 1:
            V = convolve(V, dg.reshape([1, w + 1]), "same")
        else:
            V = convolve(V, g.reshape([1, w + 1]), "same")

        if d == 2:
            V = convolve(V, dg.reshape([w + 1, 1]), "same")
        else:
            V = convolve(V, g.reshape([w + 1, 1]), "same")

    return V


def r_to_size(r_interval, size_type):
    """Convert radius interval to size (Area/Volume) interval.

    Args:
      r_interval (tuple[float]): radius interval.
      size_type (int): segmentation type (according to pygpseq.const).

    Returns:
      tuple(float): size (volume of area) interval.
    """

    if const.SEG_3D == size_type:
        # Volume interval
        o_interval = (4 / float(3)) * np.pi
        o_interval *= np.power(np.array(r_interval), 3)
    else:
        # Area interval
        o_interval = np.pi * np.square(np.array(r_interval))

    return o_interval


def rotate3d(coords, theta, axis):
    """
    Rotate coordinates around an axis.

    Args:
        coords (nd.array): coordinate table.
        theta (float): rotation angle.
        axis (int): rotation axis, axis order matches coordinate columns.

    Returns:
        tuple: rotated coordinates.
    """

    rotation_mat = False

    if axis == 0:  # X axis rotation
        rotation_mat = np.matrix(
            [
                [1, 0, 0],
                [0, np.cos(theta), -np.sin(theta)],
                [0, np.sin(theta), np.cos(theta)],
            ]
        )

    if axis == 1:  # Y axis rotation
        rotation_mat = np.matrix(
            [
                [np.cos(theta), 0, np.sin(theta)],
                [0, 1, 0],
                [-np.sin(theta), 0, np.cos(theta)],
            ]
        )

    if axis == 2:  # Z axis rotation
        rotation_mat = np.matrix(
            [
                [np.cos(theta), -np.sin(theta), 0],
                [np.sin(theta), np.cos(theta), 0],
                [0, 0, 1],
            ]
        )

    if axis < 0 or axis > 2:  # Unrecognized axis
        return ()

    transformed_mat = rotation_mat * coords
    a, b, c = transformed_mat.A

    return (a, b, c)


def round_unicode(n, nsig):
    """Round operation on unicode number in scientific notation.

    Args:
      n (unicode): number in scientific notation.
      nsig (int): number of significant digits.

    Returns:
      unicode: rounded unicode number in scientific notation.
    """

    # Convert unicode to string
    n = str(n.replace(u"\u2212", "-"))

    # Split on the exponent
    n = n.split("e") if "e" in n else [n]
    # Round the float part
    n[0] = str(round(float(n[0]), nsig))

    # Re-join with the exponent and return
    return unicode("e".join(n))


def smooth_gaussian(x, y, sigma_smooth=None, nbins=None):
    """Smoothen a curve.

    Args:
      x (numeric): x coordinates.
      y (numeric): y coordinates.
      nbins (int): curve precision (opt, def: 500).
      sigma_smooth (float): smoothing factor (opt, def: 0.01).

    Returns:
      np.array: smoothened curve.
    """

    # SET PARAMS ===============================================================

    if nbins is None:
        nbins = 500

    if sigma_smooth is None:
        sigma_smooth = 0.1

    # SMOOTHEN =================================================================

    # Evenly sampled domain
    xs = np.linspace(0, max(x), nbins)

    # Initialize output
    ynum = np.zeros(len(xs))
    ysum = np.zeros(len(xs))

    # Weighted moving average
    for i in [i for i in range(len(x)) if not np.isnan(y[i])]:
        norm = get_norm_pdf(x[i], sigma_smooth, xs)
        ynum += norm
        ysum += norm * y[i]

    # Output
    return ysum / ynum


def smooth_sparse_gaussian(
    x, y, nbins=None, sigma_smooth=None, rescale_sigma=None, **kwargs
):
    """Produce a smooth approximation of sparse data.
    Basically a smoothened binned distribution.

    Args:
      x (float): x coordinates.
      y (float): y coordinates.
      nbins (int): curve precision (opt, def: 200).
      sigma_smooth (float): smoothing factor (opt, def: 0.01).
      rescale_sigma (bool): whether to multiply sigma_smooth to max(x).

    Returns:
      dict: various metrics profiles (mean, median, mode, std).
    """

    if nbins is None:
        nbins = 200
    if sigma_smooth is None:
        sigma_smooth = 0.01
    if rescale_sigma is None:
        rescale_sigma = True

    if rescale_sigma:
        sigma_smooth *= max(x)

    # Bin data
    data = binned_profile(x, y, nbins)

    # Prepare output
    out = {"x": data["breaks"].tolist(), "n": data["n"].tolist()}

    # Smoothen profiles
    for field in ["mean", "median", "mode", "std", "max"]:
        out[field + "_raw"] = data[field]
        out[field] = smooth_gaussian(data["breaks"], data[field], sigma_smooth, nbins)

    # Output
    return out


def wilcox_sets(df, groupkey, setkey):
    """Perform Wilcoxon-Mann-Whitney U test on the provided list of sets.

    Args:
      df (pandas.DataFrame): data frame with distributions to be compared.
      groupkey (string): df column of set labels.
      setkey (string): df column with distributions to be compared.

    Returns:
      np.array: dataset with WMW U-test p-value of compared distribution.

    Examples:
      >>> dd = [('condition', 'S100'), ('y', 'float')]
      >>> p = np.array([('c1', 1), ('c2', 2)], dtype = dd)
      >>> p = pd.DataFrame(p)
      >>> print(wilcox_sets(p, 'condition', 'y'))
      [('y', 'c2', 'c1', 1.0, '')]

    """

    # Identify sets
    set_names = [c for c in set(df[groupkey])]
    n_sets = len(set_names)
    grouped = df.groupby(groupkey)

    # Initialize output
    dtype_definition = [
        ("field", "S100"),
        ("i", "S100"),
        ("j", "S100"),
        ("p", "float"),
        ("sig", "S5"),
    ]
    p_vals = np.zeros((n_sets * (n_sets - 1) / 2,), dtype=dtype_definition)

    # Cycle counter
    c = 0
    for i in range(n_sets):
        for j in range(i + 1, n_sets):
            # Run Wilcoxon-Mann-Whitney U test
            p = stats.mannwhitneyu(
                grouped.get_group(set_names[i])[setkey],
                grouped.get_group(set_names[j])[setkey],
                alternative="two-sided",
            ).pvalue

            # Significance string
            sig = ""
            if p <= 0.0001:
                sig = "***"
            elif p <= 0.001:
                sig = "**"
            elif p <= 0.01:
                sig = "*"
            elif p <= 0.05:
                sig = "."

            # Append result
            p_vals[c] = np.array(
                (setkey, set_names[i], set_names[j], p, sig), dtype=dtype_definition
            )

            # Increase cycle counter
            c += 1

    # Output
    return p_vals


# END ==========================================================================

################################################################################

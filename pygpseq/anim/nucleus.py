# -*- coding: utf-8 -*-

"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
@description: contains Nucleus wrapper.
"""

# DEPENDENCIES =================================================================

import numpy as np
from scipy import ndimage as ndi
from scipy.ndimage.interpolation import shift
from scipy.ndimage.measurements import center_of_mass
from scipy.ndimage.morphology import distance_transform_edt
import skimage.io as io
from skimage.measure import label, mesh_surface_area
from skimage.measure import marching_cubes_lewiner as marching_cubes
import warnings

from pygpseq import const
from pygpseq.tools import distance as dist, io as iot, image as imt
from pygpseq.tools import stat as stt, string as st, vector as vt

# CLASSES ======================================================================


class Nucleus(iot.IOinterface):
    """Nucleus wrapper.

    Attributes:
            __version__ (string): package version.
            s (int): series id (1-indexed).
            n (int): nucleus id (1-indexed).
            box (tuple(float)): nucleus bounding box corner coordinates.
            aspect (tuple[float]): relative/absolute px/vx size.
            dna_bg (float): dna channel estimated background.
            sig_bg (float): signal channel estimated background.
            flat_size (int): nucleus size in Z projection.
            size (int): number of px/vx in the nucleus.
            unit (string): px or vx.
            surf (float): nucleus mesh surface.
            sumI (float): sum of the intensity of every px/vx in the nucleus.
            flat_sumI(float): sum of intensity over sum projection.
            meanI (float): mean of the intensity of every px/vx in the nucleus.
            shape (float): nucleus shape descriptor.
            thr (float): global threshold used to identify the nucleus.
    """

    __version__ = const.VERSION
    c = 0
    s = 0
    n = 0
    box = ()
    box_origin = ()
    box_sides = ()
    box_mass_center = ()
    aspect = (1, 1, 1)
    shift = np.array((0, 0, 0))
    dna_bg = 0
    sig_bg = 0
    flat_size = 0
    size = 0
    unit = ""
    surf = 0
    sumI = 0
    flat_sumI = 0
    meanI = 0
    shape = 0
    thr = 0

    def __init__(
        self,
        logpath,
        n,
        series_id,
        mask,
        i,
        thr,
        offset,
        aspect,
        dna_bg,
        sig_bg,
        calc_n_surface=None,
        cond_name=None,
        **kwargs
    ):
        """Run IOinterface __init__ method.

        Args:
        logpath (string): path to the log file.
        n (int): nucleus id (1-indexed).
        series_id (int): series id (1-indexed).
        mask (numpy.array): binary image.
        i (numpy.array): image.
        thr (uint16): threshold obtained with Otsu's method.
        offset (tuple[int]): dimensions box/square offset.
        aspect (tuple[float]): pixel/voxel dimension proportion.
        dna_bg (uint16): median background for DNA channel.
        sig_bg (uint16): median background for Signal channel.
        calc_n_surface (bool): True to calculate the nucleus mesh surface.
                                                         Optional, defaults to True.
        cname (str): condition name.
        **kwargs
        """

        # parent class __init__
        super(Nucleus, self).__init__(path=logpath, append=True)

        # Default values
        if calc_n_surface is None:
            calc_n_surface = True

        # Store parameters locally
        self.c = "%s" % cond_name if type(None) != type(cond_name) else ""
        self.s = series_id
        self.n = n
        self.box = self.get_bounding_box(mask, offset)
        self.thr = thr
        self.dna_bg = dna_bg
        self.sig_bg = sig_bg
        self.aspect = aspect

        # Apply box selection to the image
        i = imt.apply_box(i, self.box)
        mask = imt.apply_box(mask, self.box)
        if "sigMask" in kwargs.keys():
            sigMask = imt.apply_box(kwargs["sigMask"], self.box)

            # Select largest object only
            L = label(sigMask)
            if L.max() > 1:
                sizes = imt.get_objects_xysize(L)
                sigMask = L == sizes.index(max(sizes)) + 1

                com = np.array(center_of_mass(mask))
                sig_com = np.array(center_of_mass(sigMask))
                self.shift = com - sig_com
            elif L.max() == 0:
                msg = "Segmentation failed in signal channel,"
                msg += " no shift correction [%d.%d]." % (self.s, self.n)
                self.printout(msg, -1)

        # Nuclear measurements
        self.size = mask.sum()
        self.flat_size = mask.max(0).sum() if len(i.shape) == 3 else self.size
        self.unit = imt.get_unit(i.shape)
        self.sumI = i[mask == 1].sum()
        self.meanI = self.sumI / self.size

        flat_mask = imt.mk_z_projection(mask, const.MAX_PROJ)
        self.flat_sumI = imt.mk_z_projection(i, const.SUM_PROJ)
        self.flat_sumI = self.flat_sumI[flat_mask == 1].sum()

        self.shape = imt.describe_shape(mask, self.aspect)
        if len(mask.shape) == 3 and calc_n_surface:
            self.surf = imt.calc_surface(mask, self.aspect)
        else:
            self.surf = self.size

        self.box_origin = np.array([c[0] + 1 for c in self.box])
        self.box_sides = np.array([np.diff(c) for c in self.box])
        self.box_mass_center = center_of_mass(mask)

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

    def check_box_offset(self, shape, offset=None):
        """Check bounding box offset.

        Note:
        If no offset is specified, it defaults to 0. If only one offset is
        specified, it will be used for every dimension. If the number of
        offsets specified does not match the number of dimensions, onlyt the
        first will be used for every dimension.

        Args:
        shape (tuple[int]): image shape.
        offset (tuple[int]): bounding box offset in px/vx [Z Y X].

        Returns:
        tuple[int]: corrected box offset.
        """

        if offset is None:
            offset = 0

        # Make offset into a list
        if type([]) != type(offset):
            offset = list(offset)

        # Identify the offset for every dimension
        if len(offset) != len(shape):
            offset = [offset[0] for d in shape]

        # Output
        return offset

    def export(self, **kwargs):
        """Export nuclear data."""

        # Set output suffix
        if "suffix" not in kwargs.keys():
            suffix = ""
        else:
            suffix = st.add_leading_dot(kwargs["suffix"])

        # Get nuclear data
        data, dp, vp, log = self.get_data(**kwargs)

        # Export as csv file
        out_fname = kwargs["series_name"] + ".nucleus" + str(self.n)
        out_fname += suffix + ".csv"
        np.savetxt(
            kwargs["out_dir"] + fname,
            data,
            header=",".join([h for h in data.dtype.names]),
            delimiter=",",
            comments="",
        )

    def get_2d_bounding_box(self, mask, offset=None):
        """Return the bounding box (2d) of the object in mask.

        Note:
        An offset can be specified for each dimension. If no offset is
        specified, it defaults to 0. If only one offset is specified, it
        will be used for every dimension. If the number of offsets specified
        does not match the number of dimensions, onlyt the first will be
        used for every dimension.

        Args:
        mask (np.array): thresholded image.
        offset (tuple[int]): bounding box offset in px/vx [Z Y X].

        Returns:
        list[int]: bounding square corner coordinates.
        """

        # Check provided offset
        offset = self.check_box_offset(mask.shape, offset)

        # Binarize mask if it is not
        mask = mask.astype("bool").astype("uint8")

        box = []

        # Y-side boundaries
        vy = mask.max(1).tolist()
        vy_min = max(0, vy.index(1) - offset[1])
        vy.reverse()
        vy_max = min(mask.shape[1] - 1, len(vy) - vy.index(1) - 1 + offset[1])
        box.append((vy_min, vy_max))

        # X-side boundaries
        vx = mask.max(0).tolist()
        vx_min = max(0, vx.index(1) - offset[0])
        vx.reverse()
        vx_max = min(mask.shape[0] - 1, len(vx) - vx.index(1) - 1 + offset[0])
        box.append((vx_min, vx_max))

        return box

    def get_3d_bounding_box(self, mask, offset=None):
        """Return the bounding box (3d) of the object in mask.

        Note:
        An offset can be specified for each dimension. If no offset is
        specified, it defaults to 0. If only one offset is specified, it
        will be used for every dimension. If the number of offsets specified
        does not match the number of dimensions, onlyt the first will be
        used for every dimension.

        Args:
        mask (np.array): thresholded image.
        offset (tuple[int]): bounding box offset in px/vx [Z Y X].

        Returns:
        list[int]: bounding box corner coordinates.
        """

        # Check provided offset
        offset = self.check_box_offset(mask.shape, offset)

        # Binarize mask if it is not
        mask = mask.astype("bool").astype("uint8")

        # Retrieve 2D bounding box
        box = [()]
        box.extend(self.get_2d_bounding_box(mask.max(0), offset[1:2]))

        # Z-side boundaries
        vz = mask.max(1).max(1).tolist()
        vz_min = max(0, vz.index(1) - offset[0])
        vz.reverse()
        vz_max = min(mask.shape[0] - 1, len(vz) - vz.index(1) - 1 + offset[0])
        box[0] = (vz_min, vz_max)

        return box

    def get_bounding_box(self, mask, offset=None):
        """Return the bounding box (2d or 3d) of the object in mask.

        Note:
        An offset can be specified for each dimension. If no offset is
        specified, it defaults to 0. If only one offset is specified, it
        will be used for every dimension. If the number of offsets
        specified does not match the number of dimensions, onlyt the first
        will be used for every dimension.

        Args:
        mask (np.array): thresholded image.
        offset (tuple[int]): bounding box offset in px/vx [Z Y X].

        Returns:
        list[int]: bounding box corner coordinates.
        """

        if len(mask.shape) == 2:
            return self.get_2d_bounding_box(mask, offset)
        elif len(mask.shape) == 3:
            return self.get_3d_bounding_box(mask, offset)

    def get_data(
        self, dna_ch, sig_ch, an_type, aspect, debugging, part_n_erosion, **kwargs
    ):
        """Get nuclear data.

        Args:
        dna_ch (np.array): image (dimensionality based on an_type).
        sig_ch (np.array): image (dimensionality based on an_type).
        an_type (int): analysis type according to pygpseq.const.
        aspect (tuple[float]): pixel/voxel dimension proportion.
        debugging (bool): True for debugging mode.
        part_n_erosion (float): partial nucleus erosion distance threshold.
        **kwargs

        Returns:
                tuple: nuclear data, density profile, volume profile and log string.
        """

        # SET PARAMS ===========================================================

        # Set output suffix
        if not "suffix" in kwargs.keys():
            suffix = ""
        else:
            suffix = st.add_leading_dot(kwargs["suffix"])

        # Set plotting
        if not "plotting" in kwargs.keys():
            kwargs["plotting"] = True

        # RETRIEVE DATA ========================================================

        # Start log
        log = ""

        # Apply box selection to channels
        dna = imt.apply_box(imt.slice_k_d_img(dna_ch, len(self.box)), self.box)

        if 0 != np.sum(self.shift):
            log += self.printout(
                "Shifting signal channel: %s" % (self.shift.round(3).tolist()), 3
            )
            shifted = shift(sig_ch, self.shift, mode="wrap")
            sig_ch = shifted.astype(sig_ch.dtype)
        sig = imt.apply_box(imt.slice_k_d_img(sig_ch, len(self.box)), self.box)

        # Produce or select mask
        if not "mask" in kwargs.keys():
            bi = Binarize(path=self.logpath, append=True, **kwargs)
            bi.verbose = self.verbose
            mask, thr, tmp_log = bi.run(dna.copy())
        else:
            mask = imt.apply_box(kwargs["mask"], self.box)

        # Select largest object only
        L = label(mask)
        if 1 < L.max():
            sizes = imt.get_objects_xysize(L)
            mask = L == sizes.index(max(sizes)) + 1
        elif 0 == L.max():
            msg = "Found empty nucleus"
            msg += " [%d.%d]." % (self.s, self.n)
            self.printout(msg, -1)

        # Apply mask to boxes
        dna[mask == 0] = 0
        sig[mask == 0] = 0

        if const.AN_MID == an_type and 3 == len(mask.shape):
            # Identify middle section
            if "mid_type" in kwargs.keys():
                mid = imt.get_mid_section_idx(dna, mask, kwargs["mid_type"])
            else:
                mid = imt.get_mid_section_idx(dna, mask)

            # Select only mid-section
            mask = mask[mid, :, :]
            dna = dna[mid, :, :]
            sig = sig[mid, :, :]

        # Perform distance transform
        laminD, centrD = dist.calc_nuclear_distances(kwargs["dist_type"], mask, aspect)

        # Export single-nucleus images in debugging mode
        if debugging:
            fname = kwargs["out_dir"] + const.OUTDIR_DEBUG
            fname += "s" + str(self.s) + "n" + str(self.n)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                if kwargs["plotting"]:
                    io.imsave("%s%s%s.tif" % (fname, self.c, suffix), mask.astype("u4"))
                if kwargs["plotting"]:
                    io.imsave(
                        "%s%s.laminD%s.tif" % (fname, self.c, suffix),
                        laminD.astype(np.uint32),
                    )
                if kwargs["plotting"]:
                    io.imsave(
                        "%s%s.centrD%s.tif" % (fname, self.c, suffix),
                        centrD.astype(np.uint32),
                    )
                if kwargs["plotting"]:
                    io.imsave(
                        "%s%s.dna%s.tif" % (fname, self.c, suffix),
                        dna.astype(np.uint32),
                    )
                if kwargs["plotting"]:
                    io.imsave(
                        "%s%s.sig%s.tif" % (fname, self.c, suffix),
                        sig.astype(np.uint32),
                    )

        # Select pixels for partial 3D nuclear analysis
        sm = np.zeros(mask.shape, dtype="u4")
        # Convert image into a list
        mask_flat = mask.reshape([np.prod(mask.shape)])
        mask_flat = mask_flat.tolist()
        mask_flat = [i for i in range(len(mask_flat)) if 1 == mask_flat[i]]

        # Prepare output
        data = np.zeros(len(mask_flat), dtype=const.DTYPE_NUCLEAR_DATA)

        # Flatten data for export
        data["dna"] = vt.flatten_and_select(dna, mask_flat)
        data["sig"] = vt.flatten_and_select(sig, mask_flat)
        data["lamin_d"] = vt.flatten_and_select(laminD, mask_flat)
        data["centr_d"] = vt.flatten_and_select(centrD, mask_flat)
        data["part"] = vt.flatten_and_select(sm, mask_flat)
        data["n"] = [self.n for i in data["dna"]]

        # Remove background
        data["dna"] = np.array(data["dna"])
        data["dna"][data["dna"] < self.dna_bg] = self.dna_bg
        data["dna"] = np.array(data["dna"]) - self.dna_bg
        data["sig"] = np.array(data["sig"])
        data["sig"][data["sig"] < self.sig_bg] = self.sig_bg
        data["sig"] = np.array(data["sig"]) - self.sig_bg

        # Add normalized distance
        laminD_norm = dist.normalize_nuclear_distance(
            kwargs["dist_type"], laminD, centrD
        )
        data["lamin_dnorm"] = vt.flatten_and_select(laminD_norm, mask_flat)

        # Prepare density profile
        density_profile, volume_profile = self.calc_density_profile(
            data["dna"], data["lamin_dnorm"], kwargs["nbins"]
        )

        # Output
        return (data, density_profile, volume_profile, log)

    def get_summary(self):
        """Get nuclear summary."""

        # Output
        data = [
            self.s,
            self.n,
            self.flat_size,
            self.size,
            self.surf,
            self.sumI,
            self.flat_sumI,
            self.meanI,
            self.shape,
        ]
        data.extend([c[0] + 1 for c in self.box])
        data.extend([c[1] + 1 for c in self.box])
        data.extend([x + 1 for x in self.box_mass_center])

        if 3 == len(self.box):
            DTYPE_NUCLEAR_SUMMARY = const.DTYPE_NUCLEAR_SUMMARY_3D
        else:
            DTYPE_NUCLEAR_SUMMARY = const.DTYPE_NUCLEAR_SUMMARY_2D

        data = np.array(tuple(data), dtype=DTYPE_NUCLEAR_SUMMARY)

        return data

    def calc_density_profile(self, dna, dnorm, nbins=200):
        """Build nucleus density profile.

        Args:
                dna (np.ndarray): single-voxel intensity array.
                dnorm (np.ndarray): single voxel normalized lamin distance array.
                nbins (int): number of bins over normalized lamin distance.
        """

        density_profile = [self.c, self.s, self.n]
        volume_profile = [self.c, self.s, self.n]

        # Prepare denominators
        M = dna.shape[0]  # Nuclear voxel count
        sumI = np.nansum(dna)  # Nuclear voxel intensity sum
        denom = sumI / M  # Average voxel intensity

        # Calculate for each bin
        breaks = np.linspace(0, 1, nbins + 1)
        for i in range(1, len(breaks)):

            # Identify voxels in the bin
            layerN = dnorm >= breaks[i - 1] if i == 1 else dnorm > breaks[i - 1]
            layerN = np.logical_and(layerN, dnorm <= breaks[i])
            Nvx = np.nansum(layerN.astype("i"))
            volume_profile.append(Nvx)

            if 0 == Nvx:
                density_profile.append(np.nan)
            else:
                # Build profile
                numer = np.nansum(dna[layerN]) / Nvx
                density_profile.append(numer / denom)

        return (np.array(density_profile), np.array(volume_profile))


# END ==========================================================================

################################################################################

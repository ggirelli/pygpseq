# -*- coding: utf-8 -*-

"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
@description: contains Analyzer wrapper, alongside all analysis-related
              parameters and methods.
"""

# DEPENDENCIES =================================================================

from pygpseq import const
from pygpseq.tools import path as pt
from pygpseq.tools import io as iot

# CLASSES ======================================================================


class Analyzer(iot.IOinterface):
    """GPSeq image Analyzer wrapper.

    Args:
      __version__ (string): class-bound attribute.
      dna_names (list[string]): dna channel names.
      sig_names (list[string]): signal channel names.
      adaptive_neighbourhood (int): adaptive threshold neighbourhood square
                                    side. If even, it's increased by one.
                                    If <= 1, turned off.
      aspect (tuple[float]): relative/absolute px/vx size.
      umes (str): unit of measure for the aspect.
      radius_interval (tuple[float]): nuclear radius interval.
      min_z_size (float): minimum Z size of the nucleus. If > 1, the ceil is
                          taken as the minimum number of slices. If < 1, the
                          float is taken as the fraction of the stack.
      seg_type (int): segmentation type according to pygpseq.const.
      an_type (int): analysis type according to pygpseq.const.
      mid_type (int): mid-section definition according to pygpseq.const.
      nsf (tuple[int]): nuclear selection features according to pygpseq.const.
                        Uses up to three features. The first two features are
                        used in the summary scattersplot (i.e., order matters).
      offset (tuple[int]): bounding box offset in px/vx [Z Y X].
      part_n_erosion float: partial nucleus erosion distance threshold.
      calc_n_surface (bool): True to calculate the nuclei mesh surface.
      sigma_density (float): sigma for smoothing.
      sigma_smooth (float): sigma for density calculation.
      nbins (int): number of bins (precision) for profile calculation.
      do_clear_Z_borders (bool): True to clear Z borders.
      rescale_deconvolved (bool): True to rescale deconvolved images.
      correctCA (bool): True to correct for chromatic aberrations (TODO).
      normalize_distance (bool): True to use relative distance from lamina.
      cdescr (dict): dictionary with better condition descriptions. The keys
                     are the condition subfolder names, the values the
                     descriptions.
    """

    __version__ = const.VERSION
    dna_names = ("dapi",)
    sig_names = ("tmr", "cy5")
    adaptive_neighbourhood = 101
    aspect = (1.0, 1.0, 1.0)
    umes = "nm"
    radius_interval = (10.0, float("inf"))
    min_z_size = 0.25
    seg_type = const.SEG_SUM_PROJ
    mask_folder = None
    mask2d_folder = None
    mask_prefix = "mask_"
    labeled = False
    compressed = False
    an_type = const.AN_SUM_PROJ
    mid_type = const.MID_SEC_DEFAULT
    dist_type = const.LD_DEFAULT
    nsf = (const.NSEL_SIZE, const.NSEL_SUMI, const.NSEL_SHAPE)
    offset = (0, 5, 5)
    part_n_erosion = 0.5
    calc_n_surface = False
    sigma_density = 0.1
    sigma_smooth = 0.1
    nbins = 200
    do_clear_Z_borders = False
    do_fill_holes = True
    correct_shift = False
    rescale_deconvolved = False
    correctCA = False
    normalize_distance = True
    cdescr = {}

    def __init__(self):
        """Run IOinterface __init__ method."""
        super(Analyzer, self).__init__()

    def __setattr__(self, name, value):
        """Check the attribute and set it."""
        self.check_attr(name, value)
        return super(Analyzer, self).__setattr__(name, value)

    def check_attr(self, name, value):
        """Run attribute format and value asserts.

        Args:
          name (string): attribute name
          value: attribute value
        """

        # Default answer
        checked = True

        if name in ["dna_names", "sig_names"]:
            assert_msg = '"%s" must be a non-empty tuple of strings.' % name
            assert type(()) == type(value), assert_msg
            assert 0 != len(value), assert_msg
            assert all(type("") == type(s) for s in value), assert_msg

        elif name in ["aspect", "radius_interval"]:
            assert_msg = '"%s" must be a non-empty tuple of floats.' % name
            assert type(()) == type(value), assert_msg
            assert 0 != len(value), assert_msg
            assert all(type(0.0) == type(s) for s in value), assert_msg

        elif name == "min_z_size":
            assert_msg = '"%s" must be a float lower than 1 or' % name
            assert_msg += " an integer greater than 1."
            assert type(value) in [type(0), type(0.0)], assert_msg
            if type(0) == type(value):
                assert value >= 1, assert_msg
            elif type(0.0) == type(value):
                assert value <= 1 and value >= 0, assert_msg

        elif name == "seg_type":
            # Check that it is one of the allowed constants
            seg_types = [const.SEG_SUM_PROJ, const.SEG_MAX_PROJ, const.SEG_3D]
            assert_msg = '"%s" must be one of the following values: ' % name
            assert_msg += str(seg_types)
            assert value in seg_types, assert_msg

        elif name == "an_type":
            # Check that it is one of the allowed constants
            an_types = [const.AN_SUM_PROJ, const.AN_MAX_PROJ, const.AN_3D, const.AN_MID]
            assert_msg = '"%s" must be one of the following values: ' % name
            assert_msg += str(an_types)
            assert value in an_types, assert_msg

        elif name == "nsf":
            assert_msg = '"%s" must be a tuple of the following values: %s' % (
                name,
                str(range(len(const.NSEL_FIELDS))),
            )
            assert type(()) == type(value), assert_msg
            assertc = all(v in range(len(const.NSEL_FIELDS)) for v in value)
            assert assertc, assert_msg

        elif name == "offset":
            assert_msg = '"%s" must be a non-empty tuple of integers.' % name
            assert type(()) == type(value), assert_msg
            assert 0 != len(value), assert_msg
            assert all(type(0) == type(s) for s in value), assert_msg

        elif name in ["sigma_smooth", "sigma_density"]:
            assert_msg = '"%s" must be a positive float.'
            assert type(0.0) == type(value), assert_msg
            assert 0 < value, assert_msg

        elif name in ["rm_z_tips", "rescale_deconvolved", "correctCA"]:
            assert_msg = '"%s" must be a boolean.'
            assert type(True) == type(value), assert_msg

        elif name == "cdescr":
            assert_msg = '"%s" must be a dictionary with string values.'
            assert type({}) == type(value), assert_msg
            assertc = all(type("") == type(v) for v in value.values())
            assert assertc, assert_msg

    def check_anseg_types(self):
        """Check seg_type and an_type."""

        # 2D segmentation does not allow 3D analysis
        no3d_cond = self.an_type == const.AN_3D
        no3d_cond = no3d_cond and self.seg_type != const.SEG_3D
        if no3d_cond:
            # Revert analysis to default
            msg = "3D analysis is not available for 2D segmentation.\n"
            msg += "Using sum z-projection analysis instead...\n"
            self.printout(msg, -1)
            self.an_type = const.AN_SUM_PROJ

        # 3D segmentation does not allow 2D analysis
        no3d_cond = self.an_type not in [const.AN_3D, const.AN_MID]
        no3d_cond = no3d_cond and self.seg_type == const.SEG_3D
        if no3d_cond:
            # Revert analysis to default
            msg = "3D segmentation is not available for 2D analysis.\n"
            msg += "Using sum z-projection segmentation instead...\n"
            self.printout(msg, -1)
            self.seg_type = const.SEG_SUM_PROJ

        # 2D segmentation does not allow mid-section analysis
        nomid_cond = self.an_type == const.AN_MID
        nomid_cond = nomid_cond and self.seg_type != const.SEG_3D
        if nomid_cond:
            # Revert analysis to default
            msg = "Mid-section analysis is not available for 2D segmentation.\n"
            msg += "Using sum z-projection analysis instead...\n"
            self.printout(msg, -1)
            self.an_type = const.SEG_SUM_PROJ


# END ==========================================================================

################################################################################

# -*- coding: utf-8 -*-

'''
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
@description: contains Analyzer wrapper, alongside all analysis-related
              parameters and methods.
'''

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
    dna_names = ('dapi',)
    sig_names = ('tmr', 'cy5')
    adaptive_neighbourhood = 101
    aspect = (1., 1., 1.)
    umes = "nm"
    radius_interval = (10., float('inf'))
    min_z_size = .25
    seg_type = const.SEG_SUM_PROJ
    an_type = const.AN_SUM_PROJ
    mid_type = const.MID_SEC_DEFAULT
    nsf = (const.NSEL_SIZE, const.NSEL_SUMI, const.NSEL_SHAPE)
    offset = (0, 5, 5)
    part_n_erosion = .5
    calc_n_surface = False
    sigma_density = .1
    sigma_smooth = .1
    nbins = 200
    do_clear_Z_borders = False
    rescale_deconvolved = False
    correctCA = False
    normalize_distance = True
    cdescr = {}

    def __init__(self):
        """Run IOinterface __init__ method. """
        super(Analyzer, self).__init__()

    def __setattr__(self, name, value):
        """Check the attribute and set it. """

        # Check the attribute
        check = self.check_attr(name, value)

        if True == check:
            # Set the attribute
            return(super(Analyzer, self).__setattr__(name, value))
        else:
            # Don't set the attribute
            return(None)

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
        checked = True

        if name in ['dna_names', 'sig_names']:
            # Require tuple
            if not type(()) == type(value):
                checked = False
            # Require non-empty tuple
            elif 0 == len(value):
                checked = False
            else:
                # Require tuple of strings
                types = [type('') == type(s) for s in value]
                types = all(types)
                if not types:
                    checked = False

            if not checked:
                msg = '"' + name + '" must be a non-empty tuple of strings.\n'
                msg += 'Keeping previous value. [' + str(self[name]) + ']'
                self.printout(msg, -1)
        elif name in ['aspect', 'radius_interval']:
            # Require tuple
            if not type(()) == type(value):
                checked = False
            # Require non-empty tuple
            elif 0 == len(value):
                checked = False
            else:
                # Require tuple of floats
                types = [type(1.) == type(s) for s in value]
                types = all(types)
                if not types:
                    checked = False

            if not checked:
                msg = '"' + name + '" must be a non-empty tuple of floats.\n'
                msg += 'Keeping previous value. [' + str(self[name]) + ']'
                self.printout(msg, -1)
        elif 'min_z_size' == name:
            # Allow float or integer
            if not type(value) in [type(0), type(0.)]:
                checked = False

            if not checked:
                msg = '"' + name + '" must be a float lower than 1 or'
                msg += ' an integer greater than 1.\n'
                msg += 'Keeping previous value. [' + str(self[name]) + ']'
                self.printout(msg, -1)
        elif 'seg_type' == name:
            # Check that it is one of the allowed constants
            seg_types = [const.SEG_SUM_PROJ, const.SEG_MAX_PROJ, const.SEG_3D]
            if not value in seg_types:
                checked = False

            if not checked:
                msg = '"' + name + '" must be one of the following values:\n'
                msg += str(seg_types) + '\n'
                msg += 'Keeping previous value. [' + str(self[name]) + ']'
                self.printout(msg, -1)
        elif 'an_type' == name:
            # Check that it is one of the allowed constants
            an_types = [const.AN_SUM_PROJ, const.AN_MAX_PROJ,
                const.AN_3D, const.AN_MID]
            if not value in an_types:
                checked = False

            if not checked:
                msg = '"' + name + '" must be one of the following values:\n'
                msg += str(an_types) + '\n'
                msg += 'Keeping previous value. [' + str(self[name]) + ']'
                self.printout(msg, -1)
        elif 'nsf' == name:
            # Require tuple
            if not type(()) == type(value):
                checked = False
            # Check that it is one of the allowed constants
            elif not all([v in range(len(const.NSEL_FIELDS)) for v in value]):
                    checked = False

            if not checked:
                msg = '"' + name + '" must be a tuple'
                msg += ' of the following values:\n'
                msg += str(range(len(const.NSEL_FIELDS))) + '\n'
                msg += 'Keeping previous value [' + str(self.nsf) + '].'
                self.printout(msg, -1)
        elif 'offset' == name:
            # Require tuple
            if not type(()) == type(value):
                checked = False
            # Require non-empty tuple
            elif 0 == len(value):
                checked = False
            else:
                # Require tuple of integers
                types = [type(0) == type(s) for s in value]
                types = all(types)
                if not types:
                    checked = False

            if not checked:
                msg = '"' + name + '" must be a non-empty tuple of integers.\n'
                msg += 'Keeping previous value. [' + str(self[name]) + ']'
                self.printout(msg, -1)
        elif name in ["sigma_smooth", "sigma_density"]:
            # Require float
            if not type(1.) == type(value):
                checked = False
            # Require positive float
            elif 0 > value:
                    checked = False

            if not checked:
                msg = '"' + name + '" must be a positive float.\n'
                msg += 'Keeping previous value. [' + str(self[name]) + ']'
                self.printout(msg, -1)
        elif name in ['rm_z_tips', 'rescale_deconvolved', 'correctCA']:
            # Require boolean
            if not type(True) == type(value):
                checked = False
            
                msg = '"' + name + '" must be a Boolean.\n'
                msg += 'Keeping previous value. [' + str(self[name]) + ']'
                self.printout(msg, -1)
        elif 'cdescr' == name:
            # Require a dictionary
            if not type({}) == type(value):
                checked = False
            else:
                # Require a dictionary with string values
                types = [type('') == type(v) for v in value.values()]
                types = all(types)
                if not types:
                    checked = False

            if not checked:
                msg = '"' + name + '" must be a dictionary with string values.\n'
                msg += 'Keeping previous value. [' + str(self.cdescr) + ']'
                self.printout(msg, -1)

        # Output
        return(checked)

    def check_anseg_types(self):
        """Check seg_type and an_type. """

        # 2D segmentation does not allow 3D analysis
        no3d_cond = self.an_type == const.AN_3D
        no3d_cond = no3d_cond and self.seg_type != const.SEG_3D
        if no3d_cond:
            # Revert analysis to default
            msg = '3D analysis is not available for 2D segmentation.\n'
            msg += 'Using sum z-projection analysis instead...\n'
            self.printout(msg, -1)
            self.an_type = const.AN_SUM_PROJ

        # 3D segmentation does not allow 2D analysis
        no3d_cond = not self.an_type in [const.AN_3D, const.AN_MID]
        no3d_cond = no3d_cond and self.seg_type == const.SEG_3D
        if no3d_cond:
            # Revert analysis to default
            msg = '3D segmentation is not available for 2D analysis.\n'
            msg += 'Using sum z-projection segmentation instead...\n'
            self.printout(msg, -1)
            self.seg_type = const.SEG_SUM_PROJ

        # 2D segmentation does not allow mid-section analysis
        nomid_cond = self.an_type == const.AN_MID
        nomid_cond = nomid_cond and self.seg_type != const.SEG_3D
        if nomid_cond:
            # Revert analysis to default
            msg = 'Mid-section analysis is not available for 2D segmentation.\n'
            msg += 'Using sum z-projection analysis instead...\n'
            self.printout(msg, -1)
            self.an_type = const.SEG_SUM_PROJ

# END ==========================================================================

################################################################################

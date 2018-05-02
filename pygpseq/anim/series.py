# -*- coding: utf-8 -*-

'''
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
@description: contains Series wrapper, which in turn contains Nucleus.
'''

# DEPENDENCIES =================================================================

import math
import os

import matplotlib.pyplot as plt
import numpy as np
from skimage.measure import label

from pygpseq import const

from pygpseq.tools.binarize import Binarize
from pygpseq.tools import io as iot
from pygpseq.tools import image as imt
from pygpseq.tools import plot
from pygpseq.tools import stat as stt
from pygpseq.tools import string as st
from pygpseq.tools import vector as vt

from pygpseq.anim.nucleus import Nucleus

# CLASSES ======================================================================

class Series(iot.IOinterface):
    """Series (Field of View, i.e., two-channel image) wrapper.

    Attributes:
      __version__ (string): package version.
      n (int): series id (1-indexed).
      name (string): series name.
      nuclei (list[pygpseq.wraps.Nuclei]): nuclei list.
      basedir (string): series folder path.
      dna_bg (float): estimated dna channel background.
      sig_bg (float): estimated signal channel background.
      flist (list): series file info.
    """

    __version__ = const.VERSION
    c = None
    n = 0
    name = ''
    nuclei = []
    basedir = '.'
    dna_bg = None
    sig_bg = None
    filist = []

    def __init__(self, ds, condition = None, **kwargs):
        """Run IOinterface __init__ method.

        Args:
          ds (dict): series information list.
          condition (pyGPSeq.wraps.Condition): condition wrapper (opt).
        """

        # If required, inherit from `condition` wrap
        if None != condition:
            logpath = condition.logpath
            super(Series, self).__init__(path = logpath, append = True)
            self.basedir = condition.path
            self.c = condition.name
        else:
            super(Series, self).__init__()
        
        # Save input parameters
        self.name = ds[0]
        self.filist = ds[1]
        self.n = ds[2]

    def __getitem__(self, key):
        """ Allow get item. """
        if key in dir(self):
            return(getattr(self, key))
        else:
            return(None)

    def __setitem__(self, key, value):
        """ Allow set item. """
        if key in dir(self):
            self.__setattr__(key, value)

    def adjust_options(self, read_only_dna = None, log = None, **kwargs):
        """Adjust options to be passed to the Nucleus class.

        Args:
          dna_names (tuple[string]): dna channel names.
          sig_names (tuple[string]): signal channel names.
          an_type (pyGPSeq.const): analysis type.

        Returns:
          dict: adds the following kwargs:
            series_name (string): series wrap name.
            basedir (string): series wrap base directory.
            dna_ch (numpy.array): image (dimensionality based on an_type).
            sig_ch (numpy.array): image (dimensionality based on an_type).
        """

        # Start log
        if None == log: log = ''

        # Only work on dna channel
        if None == read_only_dna:
            read_only_dna = False

        # Add necessary options
        kwargs['series_name'] = self.name
        kwargs['basedir'] = self.basedir

        # Read DNA channel
        kwargs['dna_ch'], log = self.get_channel(kwargs['dna_names'],
            log, **kwargs)
        if not read_only_dna:
            kwargs['sig_ch'], log = self.get_channel(kwargs['sig_names'],
                log, **kwargs)

        # Output
        return((kwargs, log))

    def export_nuclei(self, **kwargs):
        """Export current series nuclei. """

        # Set output suffix
        if not 'suffix' in kwargs.keys():
            suffix = ''
        else:
            suffix = st.add_leading_dot(kwargs['suffix'])

        # Add necessary options
        self.printout('Current series: "' + self.name + '"...', 1)
        kwargs, log = self.adjust_options(**kwargs)

        # Export nuclei
        [n.export(**kwargs) for n in self.nuclei]
        
        # Produce log
        log = np.zeros(len(self.nuclei), dtype = const.DTYPE_NUCLEAR_SUMMARY)
        for l in [n.get_summary(**kwargs) for n in self.nuclei]:
            # Append nuclear data to the series log
            summary = [self.n]
            summary.extend(l)
            log[i, :] = summary

        # Export series log
        np.savetxt(kwargs['out_dir'] + self.name + '.summary' + suffix + '.csv',
            log, delimiter = ',', comments = '',
            header = ",".join([h for h in log.dtype.names]))

        return(log)

    def find_channel(self, channel_names):
        """Return the first channel to correspond to channel_names. """

        # Fix the param type
        if type(str()) == type(channel_names):
            channel_names = [channel_names]

        # Cycle through the available channels
        for cname in channel_names:

            # Identify the requested channel
            idx = self.find_channel_id(cname)

            # Return the channel
            if -1 != idx:
                return([i for i in self.filist.items()][idx])

        # Return empty dictionary if no matching channel is found
        return({})

    def find_channel_id(self, channel_name):
        """Return the id of the channel file with the specified name. """

        # Retrieve available channel names
        names = self.get_channel_names()

        if 0 != names.count(channel_name):
            # Return matching channel id
            return(names.index(channel_name))
        else:
            # Return -1 if no matching channel is found
            return(-1)

    def find_nuclei(self, **kwargs):
        """Segment current series.

        Args:
          **kwargs
            dna_names (tuple[string]): dna channel names.
            cond_name (string): condition wrapper name.
            seg_type (pyGPSeq.const): segmentation type.
            rm_z_tips (bool): remove nuclei touching the tips of the stack.
            radius_interval (tuple[float]): allowed nuclear radius interval.
            offset (tuple[int]): dimensions box/square offset.
            aspect (tuple[float]): pixel/voxel dimension proportion.

        Returns:
          tuple: series current instance and log string.
        """

        # Set output suffix
        if not 'suffix' in kwargs.keys():
            suffix = ''
        else:
            suffix = st.add_leading_dot(kwargs['suffix'])

        # Check plotting
        if not 'plotting' in kwargs.keys():
            kwargs['plotting'] = True

        log = ""
        log += self.printout('Current series: "' + self.name + '"...', 1)

        # Read images
        kwargs, alog = self.adjust_options(read_only_dna = False, **kwargs)
        log += alog

        # Extract from kwargs
        seg_type = kwargs['seg_type']
        dna_ch = kwargs['dna_ch']
        sig_ch = kwargs['sig_ch']

        # Make new channel copy
        i = dna_ch.copy()

        # Produce a mask
        bi = Binarize(path = kwargs['logpath'], append = True, **kwargs)
        bi.verbose = self.verbose
        mask, thr, tmp_log = bi.run(i)
        log += tmp_log

        # Estimate background 
        if None == self.dna_bg:
            self.dna_bg = imt.estimate_background(dna_ch, mask, seg_type)
        kwargs['dna_bg'] = self.dna_bg
        if None == self.sig_bg:
            self.sig_bg = imt.estimate_background(sig_ch, mask, seg_type)
        kwargs['sig_bg'] = self.sig_bg
        log += self.printout('Estimating background:', 2)
        log += self.printout('DNA channel: ' + str(kwargs['dna_bg']), 3)
        log += self.printout('Signal channel: ' + str(kwargs['sig_bg']), 3)

        # Filter object size
        mask, tmp_log = bi.filter_obj_XY_size(mask)
        log += tmp_log
        mask, tmp_log = bi.filter_obj_Z_size(mask)
        log += tmp_log

        # Save mask
        log += self.printout('Saving series object mask...', 2)
        L = label(mask)

        # Plot
        fig = plt.figure()
        if 3 == len(mask.shape):
            plt.imshow(L.max(0).astype('u4'))
        else:
            plt.imshow(L.astype('u4'))
        plt.gca().get_xaxis().set_visible(False)
        plt.gca().get_yaxis().set_visible(False)
        plot.set_font_size(kwargs['font_size'])

        title = 'Nuclei in "' + kwargs['cond_name'] + '", ' + str(self.name)
        title += ' [' + str(L.max()) + ' objects]'
        plt.title(title)

        # Export as png
        fname = kwargs['out_dir'] + const.OUTDIR_MASK + kwargs['cond_name']
        fname += '.' + self.name + '.mask' + suffix + '.png'
        if kwargs['plotting']: plot.export(fname, 'png')

        # Close plot figure
        plt.close(fig)
        
        # Initialize nuclei
        log += self.printout('Bounding ' + str(L.max()) + ' nuclei...', 2)
        kwargs['logpath'] = self.logpath
        kwargs['i'] = i
        kwargs['thr'] = thr
        kwargs['series_id'] = self.n
        kwargs['cond_name'] = self.c
        seq = range(1, L.max() + 1)
        self.nuclei = [Nucleus(n = n, mask = L == n, **kwargs) for n in seq]

        return((self, log))

    def get_c(self):
        """Return number of channels in the series. """
        return(len(self.filist))

    def get_channel(self, ch_name, log = None, **kwargs):
        """Read the series specified channel.

        Args:
          ch_name (string): channel name.
          log (string): log string.
          **kwargs
        
        Returns:
          tuple: channel image and log string.
        """

        # Start log (used when verbosity is off)
        if None == log: log = ""
        log += self.printout('Reading channel "' + str(ch_name) + '"...', 2)

        # Read channel
        f = self.find_channel(ch_name)
        imch = imt.read_tiff(os.path.join(self.basedir, f[0]))
        imch = imt.slice_k_d_img(imch, 3)

        # Deconvolved images correction
        if 'rescale_deconvolved' in kwargs.keys():
            if kwargs['rescale_deconvolved']:
                # Get DNA scaling factor and rescale
                sf = imt.get_rescaling_factor(f[0], **kwargs)
                imch = (imch / sf).astype('float')
                msg = 'Rescaling "' + f[0] + '" [' + str(sf) + ']...'
                log += self.printout(msg, 3)

        # Make Z-projection
        if kwargs['an_type'] in [const.AN_SUM_PROJ, const.AN_MAX_PROJ]:
            msg = 'Generating Z-projection [' + str(kwargs['an_type']) + ']...'
            log += self.printout(msg, 3)
            if 2 != len(imch.shape):
                imch = imt.mk_z_projection(imch, kwargs['an_type'])

        # Prepare output
        return((imch, log))

    def get_channel_names(self, channel_field = None):
        """Return the names of the channels in the series. """
        if None == channel_field:
            channel_field = const.REG_CHANNEL_NAME
        return([c[channel_field] for c in self.filist.values()])

    def get_nuclei_data(self, nuclei_ids, **kwargs):
        """Retrieve a single nucleus from the current series. """

        # Read channel images
        kwargs, log = self.adjust_options(**kwargs)

        # Re-build mask
        bi = Binarize(path = self.logpath, append = True, **kwargs)
        bi.verbose = self.verbose
        mask, thr, tmp_log = bi.run(kwargs['dna_ch'].copy())
        log += tmp_log

        # Empty nuclear data array
        data = []
        densp = []
        for nucleus_id in nuclei_ids:
            # Select nucleus
            n = self.nuclei[nucleus_id -1]

            # Setup nucleus instance verbosity
            if not self.verbose:
                n.verbose = False

            # Retrieve nuclear data
            ndata, ndensp, nlog = n.get_data(mask = mask, **kwargs)

            # Update log and save nuclear data
            log += nlog
            data.append(ndata)
            densp.append(ndensp)

        return((data, densp, log))

    def propagate_attr(self, key):
        """Propagate attribute current value to every nucleus. """
        for i in range(len(self.nuclei)):
            self.nuclei[i][key] = self[key]

# END ==========================================================================

################################################################################

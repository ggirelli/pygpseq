# -*- coding: utf-8 -*-

'''
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
@description: IO management library.
'''

# DEPENDENCIES =================================================================

from datetime import datetime
import numpy as np
import os, sys, tempfile
from time import time
import warnings

from pygpseq import const

from pygpseq.tools import path as pt, string as st

# FUNCTIONS ====================================================================

class IOinterface(object):
    """Interface with input-output methods.

    Attributes:
      logpath (string): path to log file.
      verbose (bool): True to be verbose.
    """

    logpath = ''
    verbose = True

    def __init__(self, **kwargs):
        """
        Args:
          path (string): path to the log file
          append (bool): whether to append to existing log (opt)
        """

        super(IOinterface, self).__init__()

        # Append to existing log?
        if not 'append' in kwargs.keys():
            append = False
        else:
            append = kwargs['append']

        # Select path for log file
        if not 'path' in kwargs.keys():
            curpath = pt.add_trailing_slash(tempfile.gettempdir())
            curpath += 'pyGPSeq/log/'
            now = datetime.fromtimestamp(time()).strftime('%Y-%m-%d_%H:%M:%S')
            curpath += now + '.log'
        else:
            curpath = kwargs['path']

        # Create log file directory if missing
        self.check_log_dir(curpath)

        if not append:
            # Edit path if logfile exists already
            c = 1
            while os.path.exists(curpath):
                fname, fext = os.path.splitext(curpath)
                fname = fname + ' (' + str(c) + ')'
                curpath = fname + fext
                c += 1

        self.logpath = curpath
    
    def check_log_dir(self, path = None):
        """Create log file directory if missing.

        Args:
          path (string): optional log path.

        Returns:
          None: creates the folder that will contain the log, if missing.
        """

        if None == path:
            path = self.logpath

        dpath = os.path.dirname(path)
        if not os.path.isdir(dpath):
            os.makedirs(dpath)

    def gen_log_name(self):
        """Generate logfile name. """
        now = datetime.fromtimestamp(time()).strftime('%Y-%m-%d_%H:%M:%S')
        name = now + '_pyGPSeq.log'
        return(name)

    def log(self, s):
        """Write to logfile.

        Args:
          s (string): formatted message string.

        Returns:
          None: writes to logfile.
        """

        # Create log file directory if missing
        self.check_log_dir()

        # Write to logfile
        f = open(self.logpath, 'a')
        f.write(s)
        f.close()

    def printout(self, s, lvl):
        """Output for user.

        Args:
          s (string): message string.
          lvl (int): message level.

        Returns:
          string: formatted string.
        """

        if self.verbose:
            s = printout(s, lvl)
            self.log(s)
        return(printout(s, lvl, verbose = False))

def merge_profiles(profiles):
    """Formats the profiles into a single table.

    Args:
      profiles (list): list of profile dictionaries.

    Returns:
      np.array: merged profiles.
    """

    nsteps = len(profiles[0]['ratio']['x'])
    nrows = nsteps * len(profiles)

    out = np.zeros((nrows,), dtype = const.DTYPE_PROFILE_EXPORT)

    c = 0
    for profile in profiles:
        subcsv = np.zeros((nsteps,), dtype = const.DTYPE_PROFILE_EXPORT)
        subcsv['condition'] = [profile['condition'] for i in range(nsteps)]
        subcsv['x'] = profile['ratio']['x']

        for k1 in ['dna', 'sig', 'ratio']:
            for k2 in ['mean', 'median', 'mode', 'std']:
                subcsv[k1 + '_' + k2] = profile[k1][k2]
                subcsv[k1 + '_' + k2 + '_raw'] = profile[k1][k2 + '_raw']

        subcsv['n'] = [profile['n'] for i in range(nsteps)]
        out[c:(c+nsteps),] = subcsv
        c += nsteps

    return(out)

def merge_summaries(sums):
    """Formats the nuclear summaries into a single table.

    Args:
      sums (list): list of nuclear summaries, np.array.

    Returns:
      np.array: merged summaries.
    """

    nrows = sum([st.shape[0] for st in sums])
    ncols = len(sums[0].dtype)

    out = np.zeros((nrows,), dtype = const.DTYPE_NDATA_EXPORT)

    c = 0
    crow = 0
    for summary in sums:
        out[out.dtype.names[0]][crow:(crow + summary.shape[0])] = c + 1
        for name in summary.dtype.names:
            out[name][crow:(crow + summary.shape[0])] = summary[name]
        c += 1
        crow += summary.shape[0]

    return(out)

def printout(s, lvl, verbose = True, canAbort = True):
    """Log to shell.

    Args:
      s (string): message string.
      lvl (int): message level.
      verbose (bool): True to display formatted message.
      canAbort (bool): whether to abort if lvl==-2.
    """
    
    # Only a string can be logged
    if type(str()) != type(s):
        return()

    # Add level-based prefix
    if -2 == lvl:
        if not canAbort:
            s = '\n~~ ERROR ~~ ლ(ಠ益ಠლ)\n%s' % s
        if canAbort:
            print("\n~~ ERROR ~~ ლ(ಠ益ಠლ)\n%s\nTerminated.\n" % s)
            raise Exception(s)
    elif -1 == lvl:
        s = 'WARNING: %s' % s
    elif 0 == lvl:
        s = ' %s' % s
    elif 1 == lvl:
        s = '  · %s' % s
    elif 2 == lvl:
        s = '    > %s' % s
    elif 3 == lvl:
        s = '    >> %s' % s
    elif 4 <= lvl:
        s = '    >>> %s' % s

    # Log
    if verbose: print(s)
    return(st.add_trailing_new_line(s))

# END ==========================================================================

################################################################################

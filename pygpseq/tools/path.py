# -*- coding: utf-8 -*-

"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
@description: path string management library.
"""

# DEPENDENCIES =================================================================

import os
import re

from pygpseq import const

# FUNCTIONS ====================================================================


def add_extension(path, ext):
    """Add extension to path.
    If the provided path is missing the extension, add it.

    Args:
      path (string)
      ext (string): extension.

    Returns:
      string: path + extension.
    """

    # Add lading dot to extension, if needed
    ext = add_leading_dot(ext)

    # Check presence of the extension
    if not path.endswith(ext):
        path += ext

    # Output
    return path


def add_leading_dot(s):
    """Add leading dot."""
    if s[0] != ".":
        s = "." + s
    return s


def add_trailing_slash(s):
    """Add trailing slash."""
    if not s.endswith("/"):
        s = s + "/"
    return s


def select_folders(path, ext):
    """Select subdirectories containing files with the given extension.

    Args:
      path (string): base directory path.
      ext (string): extensions to look for.

    Returns:
      list: list of selected path's subdirectories containing file *.ext.
    """

    # Check input params
    path = add_trailing_slash(os.path.abspath(path))
    ext = add_leading_dot(ext)

    # Return an empty list if the provided path is not a directory
    if not os.path.isdir(path):
        return []
    # Retreive subdirectory list
    sdirs = [add_trailing_slash(path + x) for x in next(os.walk(path))[1]]

    # Will contain the selected subdirectories
    selected = []

    for sdir in sdirs:

        # Retreive file list
        flist = os.listdir(sdir)

        # Look for files with the proper extension
        if len(select_files(sdir, ext)) != 0:
            selected.append(sdir)

    selected.sort()

    return selected


def select_files(path, ext):
    """Select the files with the proper .ext in the provided path.

    Args:
      path (string): base directory path.
      ext (string): extensions to look for.

    Returns:
      list: list of selected path's files *.ext.
    """

    # Check input params
    path = add_trailing_slash(os.path.abspath(path))
    ext = add_leading_dot(ext)

    # Return an empty list if the provided path is not a directory
    if not os.path.isdir(path):
        return []
    # Retreive file list
    flist = os.listdir(path)

    # Will contain the selected files
    selected = []

    for f in flist:
        if os.path.isfile(path + f):
            # Retreive file extension
            fname, fext = os.path.splitext(f)

            # Compare file extension with .ext
            if fext == ext:
                selected.append(f)

    selected.sort()

    return selected


def select_series(flist, reg, series_field=None):
    """Group a list of files by series based on the provided regexp.

    Args:
      flist (list): list of file paths.
      reg (string): regular expression to identify series.
      series_field (string): name of the series ID field in the regexp.

    Returns:
      list: series channel paths with separated fields in a dictionary.
    """

    # Default series string field name
    if series_field is None:
        series_field = const.REG_SERIES_ID

    if type(str()) == type(reg):
        # Compile regexp
        reg = re.compile(reg)

    # Output dictionary
    d = {}

    # Select files matching the regexp
    for f in flist:

        # Search for matches
        m = re.search(reg, f)

        if m != None:
            # Identify the series string
            scur = m.groupdict()[series_field]

            # Save in the output dictionary
            if scur in d:
                # Add channel to existing series
                d[scur][f] = m.groupdict()
            else:
                # Add new series
                d[scur] = {f: m.groupdict()}

    return d


# END ==========================================================================

################################################################################

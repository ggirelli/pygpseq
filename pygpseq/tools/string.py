# -*- coding: utf-8 -*-

"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
@description: String management library.
"""

# FUNCTIONS ====================================================================


def add_leading_dot(s, not_empty=None):
    """Adds leading dot.

    Args:
      s (string): string to recieve leading dot.
      not_empty (bool): whether to add leading dot to empty string.
    """

    if not_empty is None:
        not_empty = True
    if len(s) == 0 and not_empty:
        return s
    if not s.startswith("."):
        s = "." + s
    return s


def add_trailing_new_line(s):
    """Adds trailing new line."""
    if s[len(s) - 1] != "\n":
        s += "\n"
    return s


# END ==========================================================================

################################################################################

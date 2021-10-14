# -*- coding: utf-8 -*-

"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
@description: vectors management library.
"""

# DEPENDENCIES =================================================================

import numpy as np

# FUNCTIONS ====================================================================


def flatten_and_select(v, s):
    """Flatten the array v and select rows based on s.

    Args:
      v (np.array)
      s (list): list of v rows indexes.
    """

    v = v.reshape((np.prod(v.shape))).tolist()
    v = [v[i] for i in s]
    return v


def merge_nparrays(npal):
    """Merge a list of numpy arrays into one.
    The numpy arrays MUST have the same dtype definition.

    Args:
      npal (list[np.array]): list of numpy arrays.

    Returns:
      np.array: a merged numpy array.
    """

    # Check matching dtype definitions
    if any(npal[i].dtype != npal[0].dtype for i in range(len(npal))):
        return npal

    # Count total rows
    nrows = sum(a.shape[0] for a in npal)

    # Initialize output
    merged = np.zeros((nrows,), dtype=npal[0].dtype)

    # Row pointer
    c = 0

    # Cycle through the arrays
    for i in range(len(npal)):
        # Get current array size
        nrows = npal[i].shape[0]

        # Save current array in the merged output
        merged[c : (c + nrows)] = npal[i]

        # Increment row pointer
        c += nrows

    # Output
    return merged


def rm_from_mask(L, torm):
    """Remove elements from a mask.

    Args:
      L (np.array[int]): labelled objects.
      torm (list): list of objects indexes (to remove).
    """

    if len(torm) <= L.max() - len(torm):
        # Update list of objects to be discarded
        torm = [e + 1 for e in torm]

        # Identify which objects to discard
        rm_mask = np.vectorize(lambda x: x in torm)(L)

    else:
        # Select objects to be kept
        tokeep = [e + 1 for e in range(L.max()) if e not in torm]

        # Identify which objects to discard
        rm_mask = np.vectorize(lambda x: x not in tokeep)(L)

    # Discard and re-label
    L[rm_mask] = 0
    # Output
    return L > 0


def uniquec(l):
    """Count the instances of the uniqued integers in l.

    Args:
      l (list[int]): list of integers.

    Returns:
      list[tuple]: list of (n, count(n)) for every n in unique(l).
    """

    # Possible integer values
    possible = range(max(l) + 1)

    # Count elements occurrences
    counts = [0 for i in possible]
    for n in l:
        counts[n] += 1

    # Return tupled counts
    return [(i, counts[i]) for i in possible if counts[i]]


# END ==========================================================================

################################################################################

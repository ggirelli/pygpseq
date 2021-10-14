# -*- coding: utf-8 -*-

"""
@author: Gabriele Girelli
@email: gigi.ga90@gmail.com
@description: pygpseq package.
"""

# DEPENDENCIES =================================================================

import matplotlib

matplotlib.use("ps")

__all__ = [
    "jinja2",
    "joblib",
    "matplotlib",
    "numpy",
    "pandas",
    "pickle",
    "scipy",
    "skimage",
    "weasyprint",
]

import os
import time

from joblib import Parallel, delayed

from pygpseq import anim, const, fish, tools

# END ==========================================================================

################################################################################

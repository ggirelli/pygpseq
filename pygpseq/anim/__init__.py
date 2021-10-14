# -*- coding: utf-8 -*-

"""
@author: Gabriele Girelli
@email: gigi.ga90@gmail.com
@module: pygpseq.anim
@description: wrappers for GPSeq image dataset standard analysis.
"""

# DEPENDENCIES =================================================================

__all__ = [
    "joblib",
    "logging",
    "matplotlib",
    "multiprocessing",
    "numpy",
    "pandas",
    "scipy",
    "skimage",
    "weasyprint",
]

from pygpseq.anim.main import Main
from pygpseq.anim.analysis import Analyzer
from pygpseq.anim.nucleus import Nucleus
from pygpseq.anim.series import Series
from pygpseq.anim.condition import Condition

# END ==========================================================================

################################################################################

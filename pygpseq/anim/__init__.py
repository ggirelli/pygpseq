# -*- coding: utf-8 -*-

'''
@author: Gabriele Girelli
@email: gigi.ga90@gmail.com
@module: pygpseq.anim
@description: wrappers for GPSeq image dataset standard analysis.
'''

# DEPENDENCIES =================================================================

__all__ = ['joblib', 'logging', 'matplotlib', 'multiprocessing', 'numpy',
	'pandas', 'scipy', 'skimage', 'weasyprint']

from .main import Main
from .analysis import Analyzer
from .nucleus import Nucleus
from .series import Series
from .condition import Condition

# END ==========================================================================

################################################################################

# -*- coding: utf-8 -*-

'''
@author: Gabriele Girelli
@email: gigi.ga90@gmail.com
@description: pygpseq package.
'''

# DEPENDENCIES =================================================================

__all__ = ['jinja2', 'joblib', 'matplotlib', 'numpy', 'pandas', 'pickle',
	'scipy', 'skimage', 'weasyprint']

import os
import time

from joblib import Parallel, delayed

from . import const

from .tools.binarize import Binarize
from .tools import image as imt
from .tools import io as iot
from .tools import path as pt
from .tools import string as st

from .anim.main import Main
from .anim.analysis import Analyzer
from .anim.nucleus import Nucleus
from .anim.series import Series
from .anim.condition import Condition

# END ==========================================================================

################################################################################
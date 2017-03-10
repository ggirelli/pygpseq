#!/usr/bin/python
# -*- coding: utf-8 -*-

""" Test of single pyGPSeq run. """

import os
import pkg_resources
import tempfile
from context import pygpseq as gp

gpi = gp.Main(ncores = 1)

gpi.skip = []

gpi.basedir = os.path.abspath('./tests')
gpi.outdir = tempfile.gettempdir() + '/pyGPSeq_test/'
gpi.logpath = tempfile.gettempdir() + '/pyGPSeq_test/log'

gpi.aspect = [1., 1., 1.]

gpi.reg = '^(?P<' + gp.const.REG_CHANNEL_NAME + '>[^/]*)'
gpi.reg += '\.(?P<' + gp.const.REG_CHANNEL_ID + '>channel[0-9]+)'
gpi.reg += '\.(?P<' + gp.const.REG_SERIES_ID + '>series[0-9]+)'
gpi.reg += '(?P<' + gp.const.REG_EXT + '>(_cmle)?\.tif)$'

gpi.rescale_deconvolved = False
gpi.cdescr['test'] = 'test'

# 3D test
gpi.skip = [4,5]
gpi.printout('## 3D test', 0)
gpi.seg_type = gp.const.SEG_3D
gpi.an_type = gp.const.AN_3D
gpi = gpi.run()

# Mid test
gpi.skip = [1,2,4,5]
gpi.printout('## Mid test', 0)
gpi.seg_type = gp.const.SEG_3D
gpi.an_type = gp.const.AN_MID
gpi = gpi.run()

# Max test
gpi.skip = [1,4,5]
gpi.printout('## Max test', 0)
gpi.seg_type = gp.const.SEG_SUM_PROJ
gpi.an_type = gp.const.AN_SUM_PROJ
gpi = gpi.run()

# Mix test
gpi.skip = [1,2,4,5]
gpi.printout('## Mix test [sum/max]', 0)
gpi.seg_type = gp.const.SEG_SUM_PROJ
gpi.an_type = gp.const.AN_MAX_PROJ
gpi = gpi.run()

# Sum test
gpi.skip = [1,4,5]
gpi.printout('## Sum test', 0)
gpi.seg_type = gp.const.SEG_MAX_PROJ
gpi.an_type = gp.const.AN_MAX_PROJ
gpi = gpi.run()

# Mix test
gpi.skip = [1,2,4,5]
gpi.printout('## Mix test [max/sum]', 0)
gpi.seg_type = gp.const.SEG_MAX_PROJ
gpi.an_type = gp.const.AN_SUM_PROJ
gpi = gpi.run()

# Final plots
gpi.skip = [1, 2, 3]
gpi.printout('## Final plot test', 0)
gpi.seg_type = gp.const.SEG_SUM_PROJ
gpi.an_type = gp.const.AN_MAX_PROJ
gpi = gpi.run()

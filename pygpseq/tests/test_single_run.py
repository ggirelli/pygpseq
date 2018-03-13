""" Test a single pyGPSeq run. """

from unittest import TestCase

import os
import tempfile
import pygpseq as gp

class TestSingleGP(TestCase):
	def test_init(self):
		gpi = gp.Main(ncores = 1)
		self.assertTrue(True)

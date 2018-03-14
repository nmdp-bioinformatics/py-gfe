#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#    pygfe pyGFE.
#    Copyright (c) 2017 Be The Match operated by National Marrow Donor Program. All Rights Reserved.
#
#    This library is free software; you can redistribute it and/or modify it
#    under the terms of the GNU Lesser General Public License as published
#    by the Free Software Foundation; either version 3 of the License, or (at
#    your option) any later version.
#
#    This library is distributed in the hope that it will be useful, but WITHOUT
#    ANY WARRANTY; with out even the implied warranty of MERCHANTABILITY or
#    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
#    License for more details.
#
#    You should have received a copy of the GNU Lesser General Public License
#    along with this library;  if not, write to the Free Software Foundation,
#    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA.
#
#    > http://www.fsf.org/licensing/licenses/lgpl.html
#    > http://www.opensource.org/licenses/lgpl-license.php
#



"""
test_pygfe
----------------------------------

Tests for `pygfe` module.
"""


import sys
import unittest

from pygfe.gfe import GFE
import os


class TestPygfe(unittest.TestCase):

    def setUp(self):
        self.data_dir = os.path.dirname(__file__) + "/resources"
        pass

    def tearDown(self):
        pass

    def test_000_pygfe(self):
        pygfe = GFE()
        self.assertIsInstance(pygfe, GFE)
        self.assertGreater(len(pygfe.structures), 1)
        self.assertTrue('HLA-A' in pygfe.structures)
        self.assertFalse('HLA-Z' in pygfe.structures)
        pass

    def test_001_pygfe_load(self):
        pygfe = GFE(verbose=True, load_features=True)
        self.assertIsInstance(pygfe, GFE)
        self.assertGreater(len(pygfe.structures), 1)
        self.assertGreater(len(pygfe.all_feats), 1)
        self.assertTrue('HLA-A' in pygfe.structures)
        self.assertFalse('HLA-Z' in pygfe.structures)
        pass







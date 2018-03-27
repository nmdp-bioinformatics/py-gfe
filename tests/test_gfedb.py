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

from pygfe.gfedb import GfeDB
import os

from py2neo import Graph

neo4jpass = 'gfedb'
if os.getenv("NEO4JPASS"):
    neo4jpass = os.getenv("NEO4JPASS")

neo4juser = 'neo4j'
if os.getenv("NEO4JUSER"):
    neo4juser = os.getenv("NEO4JUSER")

neo4jurl = "http://neo4j.b12x.org:80"
if os.getenv("NEO4JURL"):
    neo4jurl = os.getenv("NEO4JURL")


class TestGfeDB(unittest.TestCase):

    def setUp(self):
        self.data_dir = os.path.dirname(__file__) + "/resources"
        pass

    def tearDown(self):
        pass

    def test_000_pygfe(self):
        self.graph = Graph(neo4jurl, user=neo4juser, password=neo4jpass,
                           bolt=False)
        self.gfedb = GfeDB(self.graph)
        self.assertIsInstance(self.gfedb, GfeDB)
        pass

    # def test_001_pygfe_load(self):
    #     pygfe = pyGFE(verbose=True, load_features=True)
    #     self.assertIsInstance(pygfe, pyGFE)
    #     self.assertGreater(len(pygfe.structures), 1)
    #     self.assertGreater(len(pygfe.all_feats), 1)
    #     self.assertTrue('HLA-A' in pygfe.structures)
    #     self.assertFalse('HLA-Z' in pygfe.structures)
    #     pass







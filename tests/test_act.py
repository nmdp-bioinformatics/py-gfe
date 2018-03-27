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

import os
from Bio import SeqIO
from pygfe.act import ACT
from py2neo import Graph
from pygfe.gfedb import GfeDB
from pygfe.gfe import GFE
from pygfe.models.allele_call import AlleleCall

from seqann import BioSeqAnn
from BioSQL import BioSeqDatabase

import pymysql


neo4jpass = 'gfedb'
if os.getenv("NEO4JPASS"):
    neo4jpass = os.getenv("NEO4JPASS")

neo4juser = 'neo4j'
if os.getenv("NEO4JUSER"):
    neo4juser = os.getenv("NEO4JUSER")

neo4jurl = "http://neo4j.b12x.org:80"
if os.getenv("NEO4JURL"):
    neo4jurl = os.getenv("NEO4JURL")

biosqlpass = ''
if os.getenv("BIOSQLPASS"):
    biosqlpass = os.getenv("BIOSQLPASS")

biosqluser = ''
if os.getenv("BIOSQLUSER"):
    biosqluser = os.getenv("BIOSQLUSER")

biosqlhost = "localhost"
if os.getenv("BIOSQLHOST"):
    biosqlhost = os.getenv("BIOSQLHOST")

biosqldb = "bioseqdb"
if os.getenv("BIOSQLDB"):
    biosqldb = os.getenv("BIOSQLDB")

biosqlport = 3306
if os.getenv("BIOSQLPORT"):
    biosqlport = os.getenv("BIOSQLPORT")


def conn():
    try:
        # print(biosqlpass, biosqluser, biosqlhost, biosqldb, biosqlport, sep="\t")
        conn = pymysql.connect(host=biosqlhost,
                               port=biosqlport, user=biosqluser,
                               passwd=biosqlpass, db=biosqldb)
        conn.close()
        return True
    except Exception as e:
        print("Exception while checking MYSQL Connection:" + str(e))
        return False


class TestAct(unittest.TestCase):

    def setUp(self):
        self.data_dir = os.path.dirname(__file__) + "/resources"

        pass

    def tearDown(self):
        pass

    # def test_000_act(self):
    #     graph = Graph(neo4jurl, user=neo4juser, password=neo4jpass,
    #                   bolt=False)
    #     if conn():
    #         server = BioSeqDatabase.open_database(driver="pymysql",
    #                                               user=biosqluser,
    #                                               passwd=biosqlpass,
    #                                               host=biosqlhost,
    #                                               db=biosqldb)
    #         print("Server found!", file=sys.stderr)
    #         seqann = BioSeqAnn(server=server)
    #     else:
    #         print("No Server found!", file=sys.stderr)
    #         seqann = BioSeqAnn()
    #     act = ACT(gfedb=GfeDB(graph),
    #               seqann=seqann,
    #               gfe=GFE())

    #     self.assertIsInstance(act, ACT)
    #     pass

    # def test_001_type_sequence(self):
    #     act = ACT(self.graph)
    #     seqs = list(SeqIO.parse(self.data_dir + "/known_A.fasta", "fasta"))
    #     typing = act.type_hla("HLA-A", str(seqs[0].seq), '')
    #     print(type(typing))
    #     self.assertIsInstance(act, ACT)
    #     pass

    # def test_002_type_unknown(self):
    #     seqs = list(SeqIO.parse(self.data_dir + "/unknown_A.fasta", "fasta"))
    #     typing = self.act.type_from_seq("HLA-A", str(seqs[0].seq))
    #     self.assertIsInstance(self.act, ACT)
    #     self.assertIsInstance(typing, AlleleCall)
    #     pass

    def test_002_type_unknown(self):
        #test_id = act.get_alleleid("A*01:01:01:01")
        graph = Graph(neo4jurl, user=neo4juser, password=neo4jpass,
                      bolt=False)
        if conn():
            server = BioSeqDatabase.open_database(driver="pymysql",
                                                  user=biosqluser,
                                                  passwd=biosqlpass,
                                                  host=biosqlhost,
                                                  db=biosqldb)
            print("Server found!", file=sys.stderr)
            seqann = BioSeqAnn(server=server)
        else:
            print("No Server found!", file=sys.stderr)
            seqann = BioSeqAnn()
        act = ACT(gfedb=GfeDB(graph),
                  seqann=seqann,
                  gfe=GFE())

        seqs = list(SeqIO.parse(self.data_dir + "/unknown_A.fasta", "fasta"))
        typing = act.type_from_seq("HLA-A", str(seqs[0].seq))
        print(typing)
        print("TEST 2")
        print(type(typing))
        self.assertIsInstance(act, ACT)
        # self.assertGreater(len(pygfe.structures), 1)
        # self.assertTrue('HLA-A' in pygfe.structures)
        # self.assertFalse('HLA-Z' in pygfe.structures)
        pass

    # def test_003_pygfe(self):
    #     act = ACT(self.graph)
    #     #test_id = act.get_alleleid("A*01:01:01:01")
    #     seqs = list(SeqIO.parse(self.data_dir + "/test_seq.fasta", "fasta"))
    #     lookups = act.type_hla("HLA-A", str(seqs[0].seq), '')
    #     #print(lookups)
    #     #lookups = act.sequence_lookup("HLA-A", str(seqs[0].seq))
    #     self.assertIsInstance(act, ACT)
    #     # self.assertGreater(len(pygfe.structures), 1)
    #     # self.assertTrue('HLA-A' in pygfe.structures)
    #     # self.assertFalse('HLA-Z' in pygfe.structures)
    #     pass
    # def test_001_pygfe_load(self):
    #     pygfe = pyGFE(verbose=True, load_features=True)
    #     self.assertIsInstance(pygfe, pyGFE)
    #     self.assertGreater(len(pygfe.structures), 1)
    #     self.assertGreater(len(pygfe.all_feats), 1)
    #     self.assertTrue('HLA-A' in pygfe.structures)
    #     self.assertFalse('HLA-Z' in pygfe.structures)
    #     pass







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
#from pygfe.models.allele_call import AlleleCall
from pygfe.pygfe import pyGFE
import os

from Bio import SeqIO
#from pygfe.act import ACT
from py2neo import Graph
#from pygfe.gfedb import GfeDB
#from pygfe.gfe import GFE
from pygfe.models.error import Error
#from pygfe.models.feature import Feature
from pygfe.models.typing import Typing
from pygfe.models.feature import Feature
from pygfe.models.seqdiff import Seqdiff

from seqann.sequence_annotation import BioSeqAnn
from BioSQL import BioSeqDatabase

import pymysql
import pandas as pd
from pygfe.cypher import all_feats
from pygfe.cypher import all_gfe2hla
import pickle
import time

neo4jpass = 'gfedb'
if os.getenv("NEO4JPASS"):
    neo4jpass = os.getenv("NEO4JPASS")

neo4juser = 'neo4j'
if os.getenv("NEO4JUSER"):
    neo4juser = os.getenv("NEO4JUSER")

neo4jurl = "http://neo4j.b12x.org:80"
if os.getenv("NEO4JURL"):
    neo4jurl = os.getenv("NEO4JURL")

biosqlpass = "my-secret-pw"
if os.getenv("BIOSQLPASS"):
    biosqlpass = os.getenv("BIOSQLPASS")

biosqluser = 'root'
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

import logging
logging.basicConfig(format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
                    datefmt='%m/%d/%Y %I:%M:%S %p',
                    level=logging.INFO)


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


class TestPygfe(unittest.TestCase):

    def setUp(self):
        self.data_dir = os.path.dirname(__file__) + "/resources"
        pass

    def tearDown(self):
        pass

    def test_000_pygfe(self):
        graph = Graph(neo4jurl, user=neo4juser, password=neo4jpass,
                      bolt=False)
        #if conn():
        server = BioSeqDatabase.open_database(driver="pymysql",
                                              user=biosqluser,
                                              passwd=biosqlpass,
                                              host=biosqlhost,
                                              db=biosqldb)
        seqann = BioSeqAnn(server=server, verbose=False)
        seqann = "X"
        #else:
        #    print
        #    seqann = BioSeqAnn()
        pygfe = pyGFE(graph=graph,
                      seqann=seqann,
                      load_features=False,
                      verbose=False,
                      load_all=True,
                      loci=["HLA-A"])
        self.assertIsInstance(pygfe, pyGFE)
        seqs = list(SeqIO.parse(self.data_dir + "/unknown_A.fasta", "fasta"))
        typing = pygfe.type_from_seq("HLA-A", str(seqs[1].seq))
        #self.assertEqual(typing.gfe, 'HLA-Aw770-1-1-1-1-1-1-1-1-1-1-1-1-1-1-1-4')
        self.assertEqual(typing.hla, 'HLA-A*01:01:01:01')
        self.assertEqual(typing.status, "novel")
        self.assertIsInstance(typing, Typing)
        pass

    def test_001_pygfe(self):
        graph = Graph(neo4jurl, user=neo4juser, password=neo4jpass,
                      bolt=False)
        #if conn():

        pickle_file1 = "unique_db-feats.pickle"
        pickle_file2 = "feature-service.pickle"
        pickle_gfe2feat = "gfe2feat.pickle"
        pickle_file3 = "gfe2hla.pickle"
        pickle_file4 = "seq2hla.pickle"

        with open(pickle_gfe2feat, 'rb') as handle1:
            gfe_feats = pickle.load(handle1)

        with open(pickle_file1, 'rb') as handle1:
            feats = pickle.load(handle1)

        with open(pickle_file2, 'rb') as handle2:
            cached_feats = pickle.load(handle2)

        with open(pickle_file3, 'rb') as handle3:
            gfe2hla = pickle.load(handle3)

        with open(pickle_file4, 'rb') as handle:
            seq2hla = pickle.load(handle)

        seqann = BioSeqAnn(verbose=False, cached_features=cached_feats, align=True)

        pygfe = pyGFE(graph=graph,
                      seqann=seqann,
                      gfe_feats=gfe_feats,
                      gfe2hla=gfe2hla,
                      seq2hla=seq2hla,
                      features=feats,
                      verbose=False)
        self.assertIsInstance(pygfe, pyGFE)
        seqs = list(SeqIO.parse(self.data_dir + "/unknown_A.fasta", "fasta"))
        typing = pygfe.type_from_seq("HLA-A", str(seqs[1].seq))
        print(typing)
        #self.assertEqual(typing.gfe, 'HLA-Aw770-1-1-1-1-1-1-1-1-1-1-1-1-1-1-1-4')
        self.assertEqual(typing.hla, 'HLA-A*01:01:01:01')
        self.assertEqual(typing.status, "novel")
        self.assertIsInstance(typing, Typing)
        pass

    def test_002_pygfe(self):
        graph = Graph(neo4jurl, user=neo4juser, password=neo4jpass,
                      bolt=False)
        #if conn():

        pickle_file1 = "unique_db-feats.pickle"
        pickle_file2 = "feature-service.pickle"
        pickle_gfe2feat = "gfe2feat.pickle"
        #pickle_file3 = "gfe2hla.pickle"
        pickle_file4 = "seq2hla.pickle"

        with open(pickle_gfe2feat, 'rb') as handle1:
            gfe_feats = pickle.load(handle1)

        with open(pickle_file1, 'rb') as handle1:
            feats = pickle.load(handle1)

        with open(pickle_file2, 'rb') as handle2:
            cached_feats = pickle.load(handle2)

        # with open(pickle_file3, 'rb') as handle3:
        #     gfe2hla = pickle.load(handle3)

        with open(pickle_file4, 'rb') as handle:
            seq2hla = pickle.load(handle)

        seqann = BioSeqAnn(verbose=False, cached_features=cached_feats, align=True)

        pygfe = pyGFE(graph=graph,
                      seqann=seqann,
                      gfe_feats=gfe_feats,
                      seq2hla=seq2hla,
                      features=feats,
                      verbose=False)
        self.assertIsInstance(pygfe, pyGFE)
        seqs = list(SeqIO.parse(self.data_dir + "/unknown_A.fasta", "fasta"))
        typing = pygfe.type_from_seq("HLA-A", str(seqs[1].seq))
        print(typing)
        #self.assertEqual(typing.gfe, 'HLA-Aw770-1-1-1-1-1-1-1-1-1-1-1-1-1-1-1-4')
        self.assertEqual(typing.hla, 'HLA-A*01:01:01:01')
        self.assertEqual(typing.status, "novel")
        self.assertIsInstance(typing, Typing)
        pass


    def test_005_picklefiles(self):
        graph = Graph("http://ec2-34-207-175-160.compute-1.amazonaws.com:80", user=neo4juser, password=neo4jpass,
                      bolt=False)
        #if conn():
        server = BioSeqDatabase.open_database(driver="pymysql",
                                              user=biosqluser,
                                              passwd=biosqlpass,
                                              host=biosqlhost,
                                              db=biosqldb,
                                              port=biosqlport)
        seqann = BioSeqAnn(server=server, verbose=False)

        gfe = GFE()
        #cached_feats = gfe.all_feats
        # print("Finished loading cached_feats")
        # pickle_service = "feature-service.pickle"
        # with open(pickle_service, 'wb') as handle2:
        #     pickle.dump(cached_feats, handle2, protocol=pickle.HIGHEST_PROTOCOL)

        feat_df = pd.DataFrame(graph.data(all_feats()))
        feat_df['ID'] = feat_df.apply(lambda row: ":".join([row['DB'],
                                                            row['LOC'],
                                                            str(row['RANK']),
                                                            row['TERM'],
                                                            row['SEQ']]),
                                      axis=1)
        feats = feat_df[['ID', 'ACCESSION']].set_index('ID').to_dict()['ACCESSION']

        print("Finished loading feats")
        pickle_feats = "unique_db-feats.pickle"
        with open(pickle_feats, 'wb') as handle1:
            pickle.dump(feats, handle1, protocol=pickle.HIGHEST_PROTOCOL)

        gfedb = GfeDB(graph=graph, persist=False, verbose=False)
        act = ACT(gfedb=gfedb, seqann=seqann, load_gfe2hla=True,
                  load_gfe2feat=True,
                  load_seq2hla=True, gfe=gfe)

        print("Finished loading all!!")

        gfe2hla = act.gfe2hla
        seq2hla = act.seq2hla
        gfe2feat = act.gfe_feats

        pickle_gfe2feat = "gfe2feat.pickle"
        with open(pickle_gfe2feat, 'wb') as handle5:
            pickle.dump(gfe2feat, handle5, protocol=pickle.HIGHEST_PROTOCOL)

        pickle_gfe2hla = "gfe2hla.pickle"
        with open(pickle_gfe2hla, 'wb') as handle3:
            pickle.dump(gfe2hla, handle3, protocol=pickle.HIGHEST_PROTOCOL)

        pickle_seq2hla = "seq2hla.pickle"
        with open(pickle_seq2hla, 'wb') as handle4:
            pickle.dump(seq2hla, handle4, protocol=pickle.HIGHEST_PROTOCOL)

        pass

    def test_003_loader1(self):
        graph = Graph(neo4jurl, user=neo4juser, password=neo4jpass,
                      bolt=False)
        #if conn():
        server = BioSeqDatabase.open_database(driver="pymysql",
                                              user=biosqluser,
                                              passwd=biosqlpass,
                                              host=biosqlhost,
                                              db=biosqldb,
                                              port=3307)
        seqann = BioSeqAnn(server=server, verbose=True)

        pickle_file1 = "unique_db-feats.pickle"
        pickle_file2 = "feature-service.pickle"
        pickle_gfe2feat = "gfe2feat.pickle"
        pickle_file3 = "gfe2hla.pickle"
        pickle_file4 = "seq2hla.pickle"

        with open(pickle_gfe2feat, 'rb') as handle1:
            gfe_feats = pickle.load(handle1)

        with open(pickle_file1, 'rb') as handle1:
            feats = pickle.load(handle1)

        with open(pickle_file2, 'rb') as handle2:
            cached_feats = pickle.load(handle2)

        with open(pickle_file3, 'rb') as handle3:
            gfe2hla = pickle.load(handle3)

        with open(pickle_file4, 'rb') as handle:
            seq2hla = pickle.load(handle)

        start = time.time()
        pygfe = pyGFE(graph=graph,
                      seqann=seqann,
                      load_features=False,
                      verbose=True,
                      features=feats,
                      seq2hla=seq2hla,
                      gfe2hla=gfe2hla,
                      gfe_feats=gfe_feats,
                      cached_features=cached_feats,
                      loci=["HLA-A"])
        self.assertIsInstance(pygfe, pyGFE)
        seqs = list(SeqIO.parse(self.data_dir + "/unknown_A.fasta", "fasta"))
        typing1 = pygfe.type_from_seq("HLA-A", str(seqs[1].seq), "3.31.0")
        typing2 = pygfe.type_from_seq("HLA-A", str(seqs[1].seq), "3.20.0")
        end = time.time()
        time_taken = end - start
        self.assertEqual(typing1.hla, 'HLA-A*01:01:01:01')
        self.assertEqual(typing1.status, "novel")
        self.assertIsInstance(typing1, Typing)
        self.assertEqual(typing2.hla, 'HLA-A*01:01:01:01')
        self.assertEqual(typing2.status, "novel")
        self.assertIsInstance(typing2, Typing)
        pass

    def test_003_loader2(self):
        start = time.time()
        graph = Graph(neo4jurl, user=neo4juser, password=neo4jpass,
                      bolt=False)
        #if conn():
        server = BioSeqDatabase.open_database(driver="pymysql",
                                              user=biosqluser,
                                              passwd=biosqlpass,
                                              host=biosqlhost,
                                              db=biosqldb)
        seqann = BioSeqAnn(server=server, verbose=False)

        pickle_file1 = "unique_db-feats.pickle"
        pickle_file2 = "feature-service.pickle"
        pickle_gfe2feat = "gfe2feat.pickle"
        pickle_file3 = "gfe2hla.pickle"
        pickle_file4 = "seq2hla.pickle"
        with open(pickle_gfe2feat, 'rb') as handle1:
            gfe_feats = pickle.load(handle1)

        with open(pickle_file1, 'rb') as handle1:
            feats = pickle.load(handle1)

        with open(pickle_file2, 'rb') as handle2:
            cached_feats = pickle.load(handle2)

        with open(pickle_file3, 'rb') as handle3:
            gfe2hla = pickle.load(handle3)

        with open(pickle_file4, 'rb') as handle:
            seq2hla = pickle.load(handle)

        pygfe = pyGFE(graph=graph,
                      seqann=seqann,
                      load_features=False,
                      verbose=False,
                      features=feats,
                      seq2hla=seq2hla,
                      gfe2hla=gfe2hla,
                      gfe_feats=gfe_feats,
                      cached_features=cached_feats,
                      loci=["HLA-A"])
        self.assertIsInstance(pygfe, pyGFE)
        seqs = list(SeqIO.parse(self.data_dir + "/known_A.fasta", "fasta"))
        typing1 = pygfe.type_from_seq("HLA-A", str(seqs[0].seq), "3.20.0")
        typing2 = pygfe.type_from_seq("HLA-A", str(seqs[0].seq), "3.31.0")
        end = time.time()
        time_taken = end - start
        print("TIME TAKEN: " + str(time_taken))
        self.assertEqual(typing2.hla, 'HLA-A*01:01:01:01')
        self.assertEqual(typing2.status, "documented")
        self.assertIsInstance(typing2, Typing)
        self.assertEqual(typing1.hla, 'HLA-A*01:01:01:01')
        self.assertEqual(typing1.status, "documented")
        self.assertIsInstance(typing1, Typing)
        pass

    def test_004_loader3(self):
        start = time.time()
        graph = Graph(neo4jurl, user=neo4juser, password=neo4jpass,
                      bolt=False)
        #if conn():
        server = BioSeqDatabase.open_database(driver="pymysql",
                                              user=biosqluser,
                                              passwd=biosqlpass,
                                              host=biosqlhost,
                                              db=biosqldb, port=3307)
        seqann = BioSeqAnn(server=server, verbose=True)
        pygfe = pyGFE(graph=graph,
                      seqann=seqann,
                      verbose=False,
                      load_features=False,
                      load_gfe2hla=True,
                      load_seq2hla=True,
                      load_gfe2feat=True,
                      loci=["HLA-A"])
        self.assertIsInstance(pygfe, pyGFE)
        seqs = list(SeqIO.parse(self.data_dir + "/known_A.fasta", "fasta"))
        #typing1 = pygfe.type_from_seq("HLA-A", str(seqs[0].seq), "3.20.0")
        typing2 = pygfe.type_from_seq("HLA-A", str(seqs[0].seq), "3.31.0")
        server.close()
        end = time.time()
        time_taken = end - start
        print("TIME TAKEN: " + str(time_taken))
        self.assertEqual(typing2.hla, 'HLA-A*01:01:01:01')
        self.assertEqual(typing2.status, "documented")
        self.assertIsInstance(typing2, Typing)
        # self.assertEqual(typing1.hla, 'HLA-A*01:01:01:01')
        # self.assertEqual(typing1.status, "documented")
        # self.assertIsInstance(typing1, Typing)
        pass

    def test_001_load_features(self):
        graph = Graph(neo4jurl, user=neo4juser, password=neo4jpass,
                      bolt=False)
        #if conn():
        server = BioSeqDatabase.open_database(driver="pymysql",
                                              user=biosqluser,
                                              passwd=biosqlpass,
                                              host=biosqlhost,
                                              db=biosqldb,
                                              port=3307)
        seqann = BioSeqAnn(server=server)
        #else:
        #    seqann = BioSeqAnn()
        pygfe = pyGFE(graph=graph,
                      seqann=seqann,
                      verbose=True,
                      load_features=True,
                      verbosity=2,
                      loci=["HLA-A"])
        self.assertIsInstance(pygfe, pyGFE)
        self.assertGreater(len(pygfe.gfe.structures), 1)
        self.assertGreater(len(pygfe.gfe.all_feats), 1)
        self.assertTrue('HLA-A' in pygfe.gfe.structures)
        self.assertFalse('HLA-Z' in pygfe.gfe.structures)
        pass

    def test_004_hml(self):
        #start = time.time()
        graph = Graph(neo4jurl, user=neo4juser, password=neo4jpass,
                      bolt=False)
        #if conn():
        server = BioSeqDatabase.open_database(driver="pymysql",
                                              user=biosqluser,
                                              passwd=biosqlpass,
                                              host=biosqlhost,
                                              db=biosqldb,
                                              port=3307)
        seqann = BioSeqAnn(server=server, dbversion="3190", verbose=True)

        pickle_file1 = "unique_db-feats.pickle"
        pickle_file2 = "feature-service.pickle"
        pickle_gfe2feat = "gfe2feat.pickle"
        pickle_file3 = "gfe2hla.pickle"
        pickle_file4 = "seq2hla.pickle"
        with open(pickle_gfe2feat, 'rb') as handle1:
            gfe_feats = pickle.load(handle1)

        with open(pickle_file1, 'rb') as handle1:
            feats = pickle.load(handle1)

        with open(pickle_file2, 'rb') as handle2:
            cached_feats = pickle.load(handle2)

        with open(pickle_file3, 'rb') as handle3:
            gfe2hla = pickle.load(handle3)

        with open(pickle_file4, 'rb') as handle:
            seq2hla = pickle.load(handle)

        pygfe = pyGFE(graph=graph,
                      seqann=seqann,
                      load_features=False,
                      verbose=True,
                      features=feats,
                      seq2hla=seq2hla,
                      gfe2hla=gfe2hla,
                      gfe_feats=gfe_feats,
                      cached_features=cached_feats,
                      loci=["HLA-C"])
        self.assertIsInstance(pygfe, pyGFE)
        seqs = list(SeqIO.parse(self.data_dir + "/hml_fail.fasta", "fasta"))
        typing1 = pygfe.type_from_seq("HLA-C", str(seqs[1].seq), "3.19.0")
        #typing2 = pygfe.type_from_seq("HLA-DRB1", str(seqs[0].seq), "3.31.0")
        #typing2 = pygfe.type_from_seq("HLA-DRB1", str(seqs[0].seq), "3.31.0")
        #end = time.time()
        #time_taken = end - start
        print(typing1)
        #print("=====")
        #print(typing2)
        # self.assertEqual(typing2.hla, 'HLA-A*01:01:01:01')
        # self.assertEqual(typing2.status, "documented")
        #self.assertIsInstance(typing2, Typing)
        # self.assertEqual(typing1.hla, 'HLA-A*01:01:01:01')
        # self.assertEqual(typing1.status, "documented")
        self.assertIsInstance(typing1, Typing)
        pass

    def test_005_A(self):
        #start = time.time()
        graph = Graph(neo4jurl, user=neo4juser, password=neo4jpass,
                      bolt=False)
        #if conn():
        server = BioSeqDatabase.open_database(driver="pymysql",
                                              user=biosqluser,
                                              passwd=biosqlpass,
                                              host=biosqlhost,
                                              db=biosqldb,
                                              port=3307)
        seqann = BioSeqAnn(server=server, dbversion="3200", verbose=True)

        pickle_file1 = "unique_db-feats.pickle"
        pickle_file2 = "feature-service.pickle"
        pickle_gfe2feat = "gfe2feat.pickle"
        pickle_file3 = "gfe2hla.pickle"
        pickle_file4 = "seq2hla.pickle"
        with open(pickle_gfe2feat, 'rb') as handle1:
            gfe_feats = pickle.load(handle1)

        with open(pickle_file1, 'rb') as handle1:
            feats = pickle.load(handle1)

        with open(pickle_file2, 'rb') as handle2:
            cached_feats = pickle.load(handle2)

        with open(pickle_file3, 'rb') as handle3:
            gfe2hla = pickle.load(handle3)

        with open(pickle_file4, 'rb') as handle:
            seq2hla = pickle.load(handle)

        pygfe = pyGFE(graph=graph,
                      seqann=seqann,
                      load_features=False,
                      verbose=True,
                      features=feats,
                      seq2hla=seq2hla,
                      gfe2hla=gfe2hla,
                      gfe_feats=gfe_feats,
                      cached_features=cached_feats,
                      loci=["HLA-DQB1"])
        self.assertIsInstance(pygfe, pyGFE)
        seqs = list(SeqIO.parse(self.data_dir + "/A_fail.fasta", "fasta"))
        typing1 = pygfe.type_from_seq("HLA-DQB1", str(seqs[1].seq), "3.20.0")
        #typing2 = pygfe.type_from_seq("HLA-DRB1", str(seqs[0].seq), "3.31.0")
        #typing2 = pygfe.type_from_seq("HLA-DRB1", str(seqs[0].seq), "3.31.0")
        #end = time.time()
        #time_taken = end - start
        print(typing1)
        #print("=====")
        #print(typing2)
        # self.assertEqual(typing2.hla, 'HLA-A*01:01:01:01')
        # self.assertEqual(typing2.status, "documented")
        #self.assertIsInstance(typing2, Typing)
        # self.assertEqual(typing1.hla, 'HLA-A*01:01:01:01')
        # self.assertEqual(typing1.status, "documented")
        self.assertIsInstance(typing1, Typing)
        pass

    def test_006_align(self):

        graph = Graph(neo4jurl, user=neo4juser, password=neo4jpass,
                      bolt=False)
        #if conn():
        server = BioSeqDatabase.open_database(driver="pymysql",
                                              user=biosqluser,
                                              passwd=biosqlpass,
                                              host=biosqlhost,
                                              db=biosqldb,
                                              port=3307)
        seqann = BioSeqAnn(align=True, server=server, dbversion="3310", verbose=True)

        pickle_file1 = "unique_db-feats.pickle"
        pickle_file2 = "feature-service.pickle"
        pickle_gfe2feat = "gfe2feat.pickle"
        pickle_file3 = "gfe2hla.pickle"
        pickle_file4 = "seq2hla.pickle"
        with open(pickle_gfe2feat, 'rb') as handle1:
            gfe_feats = pickle.load(handle1)

        with open(pickle_file1, 'rb') as handle1:
            feats = pickle.load(handle1)

        with open(pickle_file2, 'rb') as handle2:
            cached_feats = pickle.load(handle2)

        with open(pickle_file3, 'rb') as handle3:
            gfe2hla = pickle.load(handle3)

        with open(pickle_file4, 'rb') as handle:
            seq2hla = pickle.load(handle)

        pygfe = pyGFE(graph=graph,
                      seqann=seqann,
                      load_features=False,
                      verbose=True,
                      features=feats,
                      seq2hla=seq2hla,
                      gfe2hla=gfe2hla,
                      gfe_feats=gfe_feats,
                      cached_features=cached_feats,
                      loci=["HLA-A"])
        self.assertIsInstance(pygfe, pyGFE)
        seqs = list(SeqIO.parse(self.data_dir + "/align_tests.fasta", "fasta"))
        typing1 = pygfe.type_from_seq("HLA-A", str(seqs[0].seq), "3.31.0")
        typing2 = pygfe.type_from_seq("HLA-A", str(seqs[1].seq), "3.31.0")
        typing3 = pygfe.type_from_seq("HLA-A", str(seqs[2].seq), "3.31.0")
        typing4 = pygfe.type_from_seq("HLA-A", str(seqs[3].seq), "3.31.0")
        self.assertEqual(typing1.hla, 'HLA-A*02:01:01:12')
        self.assertEqual(typing2.hla, 'HLA-A*02:01:01:12')
        self.assertEqual(typing3.hla, 'HLA-A*02:01:01:12')
        self.assertEqual(typing4.hla, 'HLA-A*02:01:01:12')
        #end = time.time()
        #time_taken = end - start
        #print(typing1)
        #print(typing1.aligned.keys())
        #print(typing1.novel_features)
        #difss = pygfe.hla_seqdiff("HLA-A","3.31.0","HLA-A*01:01:01:01","HLA-A*01:01:01:07")

        #self.assertIsInstance(typing1, Typing)
        pass








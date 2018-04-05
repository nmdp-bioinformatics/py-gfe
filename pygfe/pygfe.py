'''
Created on Feb 8, 2017

@author: mhalagan
'''

from seqann import BioSeqAnn
from pygfe.feature_client.models.feature import Feature

from pygfe.cypher import sequence_search
from pygfe.cypher import gfe_search
from pygfe.cypher import similar_gfe_classI
from pygfe.cypher import similar_gfe_classII
from pygfe.cypher import gfe_hla
from pygfe.cypher import hla_Ggroups
from pygfe.cypher import hla_gfe
from pygfe.cypher import hla_ars
from pygfe.cypher import gfe_ars
from pygfe.cypher import get_sequence
from pygfe.cypher import similar_kir
from pygfe.cypher import get_features
# from pygfe.cypher import ref_query
from pygfe.cypher import search_hla_features

from pygfe.cypher import hla_alleleid
from pygfe.cypher import gfe_alleleid
from pygfe.cypher import fullseqid
from pygfe.cypher import seqid
from pygfe.cypher import search_feature
from pygfe.cypher import persisted_query

from pygfe.cypher import groups_classI
from pygfe.cypher import groups_classII

from pygfe.models.error import Error
from pygfe.models.feature import Feature
from pygfe.models.typing import Typing
# from pygfe.models.gfe_call import GfeCall
# from pygfe.models.gfe_typing import GfeTyping
# from pygfe.models.allele_call import AlleleCall
# from pygfe.models.feature_call import FeatureCall
# from pygfe.models.typing_status import TypingStatus
# from pygfe.models.ars_call import ArsCall
# from pygfe.models.persisted import Persisted
# from pygfe.models.persisted_data import PersistedData

from py2neo import Node, Relationship
import pandas as pa

import os
import glob
import re
import json

import sys
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from BioSQL import BioSeqDatabase
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import OrderedDict
from Bio.Alphabet import IUPAC
from py2neo import Graph

from typing import Dict

flatten = lambda l: [item for sublist in l for item in sublist]
is_gfe = lambda x: True if re.search("\d+-\d+-\d+", x) else False
is_kir = lambda x: True if re.search("KIR", x) else False
is_classII = lambda x: True if re.search("HLA-D", x) else False
is_classI = lambda x: True if re.search("HLA-\Dw", x) else False
lc = lambda x: x.lower() if not re.search("UTR", x) else x.lower().replace("utr", "UTR")

from pandas import DataFrame
from pygfe.gfe import GFE
from pygfe.act import ACT
from pygfe.gfedb import GfeDB
from pygfe.cypher import all_gfe2hla
from pygfe.cypher import all_seq2hla
from pygfe.graph_search import GraphSearch
import logging

logging.basicConfig(format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
                    datefmt='%m/%d/%Y %I:%M:%S %p',
                    level=logging.INFO)


class pyGFE(ACT, GraphSearch):
    '''
    classdocs
    '''
    def __init__(self, url="http://feature.nmdp-bioinformatics.org",
                 loci=['HLA-A', 'HLA-B', 'HLA-C',
                       'HLA-DRB1', 'HLA-DQB1', 'HLA-DRB4',
                       'HLA-DRB5', 'HLA-DPB1', 'HLA-DPA1',
                       'HLA-DQA1', 'HLA-DRB3'],
                 load_features=False, store_features=False,
                 graph: Graph=None, persist: bool=False,
                 seqann: BioSeqAnn=None,
                 features: Dict=None,
                 cached_features: Dict=None,
                 verbose: bool=False,
                 pid: str="NA",
                 gfe2hla: Dict=None,
                 gfe_feats: DataFrame=None,
                 seq2hla: DataFrame=None,
                 verbosity=1):
        '''
        Constructor
        '''
        # TODO: Add catch if seqann or graph aren't defined
        self.gfe = GFE(store_features=store_features,
                       load_features=load_features,
                       cached_features=cached_features,
                       url=url,
                       loci=loci,
                       verbose=verbose,
                       verbosity=verbosity)
        self.gfedb = GfeDB(graph=graph, persist=persist, verbose=verbose)
        self.logger = logging.getLogger("Logger." + __name__)
        self.logname = "ID {:<10} - ".format(str(pid))
        self.seqann = seqann
        self.features = features
        self.gfe2hla = gfe2hla
        self.seq2hla = seq2hla
        self.gfe_feats = gfe_feats
        self.verbose = verbose


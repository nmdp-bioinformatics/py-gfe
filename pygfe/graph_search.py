'''
Created on Feb 8, 2017

@author: mhalagan
'''


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
import pymysql

from pygfe.gfedb import GfeDB

flatten = lambda l: [item for sublist in l for item in sublist]
is_gfe = lambda x: True if re.search("\d+-\d+-\d+", x) else False
is_kir = lambda x: True if re.search("KIR", x) else False
is_classII = lambda x: True if re.search("HLA-D", x) else False
is_classI = lambda x: True if re.search("HLA-\Dw", x) else False
lc = lambda x: x.lower() if not re.search("UTR", x) else x.lower().replace("utr", "UTR")


class GraphSearch(object):
    '''
    Example:

        >>> from Bio import SeqIO
        >>> from BioSQL import BioSeqDatabase
        >>> from seqann.sequence_annotation import BioSeqAnn
        >>> from pygfe.pygfe import pyGFE
        >>> seq_file = 'test_dq.fasta'
        >>> gfe = pyGFE()
        >>> server = BioSeqDatabase.open_database(driver="pymysql", user="root",
        ...                                       passwd="", host="localhost",
        ...                                      db="bioseqdb")
        >>> seqann = BioSeqAnn(server=server)
        >>> seq_rec = list(SeqIO.parse(seq_file, 'fasta'))[0]
        >>> annotation = seqann.annotate(seq_rec, "HLA-DQB1")
        >>> gfe = gfe.get_gfe(annotation, "HLA-DQB1")
        >>> print(gfe)
        HLA-DQB1w0-4-0-141-0-12-0-4-0-0-0-0-0
    '''

    def __init__(self, gfedb: GfeDB=None):
        '''
        Constructor
        '''
        self.gfedb = gfedb

    def get_alleleid(self, typing):
        """
        creates GFE from HLA sequence and locus

        :param locus: string containing HLA locus.
        :param sequence: string containing sequence data.

        :return: GFEobject.
        """
        if is_gfe(typing):
            gfe_q = gfe_alleleid(typing)
            gfeid = pa.DataFrame(self.graph.data(gfe_q))
            if not gfeid.empty:
                return gfeid['ID'][0]
            else:
                self.gfedb.max_alleleid += 1
                return self.gfedb.max_alleleid
        else:
            hla_q = hla_alleleid(typing)
            hlaid = pa.DataFrame(self.graph.data(hla_q))
            if not hlaid.empty:
                return hlaid['ID'][0]
            else:
                self.gfedb.max_alleleid += 1
                return self.gfedb.max_alleleid

    # TODO: Make this come from the
    #       pygfe cache instead of making
    #       a bunch of calls to the GFE DB
    def unique_features(self, features):
        """
        creates GFE from HLA sequence and locus

        :param locus: string containing HLA locus.
        :param sequence: string containing sequence data.

        :return: GFEobject.
        """
        unique = []
        for feat in features:
            feat_q = search_feature(feat.term.upper(), feat.rank, feat.sequence)
            seq_features = pa.DataFrame(self.graph.data(feat_q))
            if seq_features.empty:
                unique.append(feat)
        return unique

    def sequence_lookup(self, locus, sequence):
        """
        Looks up sequence from

        :param locus: string containing HLA locus.
        :param sequence: string containing sequence data.

        :return: GFEobject.
        """
        lookup_query = sequence_search(locus, sequence)
        sequence_data = pa.DataFrame(self.graph.data(lookup_query))
        if not sequence_data.empty:
            features = list()
            gfe = list(set([x for x in sequence_data["GFE"]]))
            hla = list(set([x for x in sequence_data["HLA"]]))
            seq_features = pa.DataFrame(self.graph.data(get_features(gfe[0])))
            for i in range(0, len(seq_features['term'])):
                feature = Feature(accession=seq_features['accession'][i], rank=seq_features['rank'][i], sequence=seq_features['sequence'][i], term=lc(seq_features['term'][i]))
                features.append(feature)

            typing = Typing(hla=hla[0], related_gfe=[GfeTyping(gfe=gfe[0], shares=features, features_shared=len(features))])
            return [typing, gfe[0], features]
        else:
            return

    def gfe_lookup(self, gfe, features):
        """
        creates GFE from HLA sequence and locus

        :param locus: string containing HLA locus.
        :param sequence: string containing sequence data.

        :return: GFEobject.
        """
        gfe_data = pa.DataFrame(self.graph.data(gfe_search(gfe)))
        if not gfe_data.empty:
            typing_list = list()
            for hla in gfe_data["HLA"]:
                typing_list.append(Typing(hla=hla, related_gfe=[GfeTyping(gfe=gfe, shares=features, features_shared=len(features))]))
            return typing_list
        else:
            return

    def get_class(self, gfe):
        """
        creates GFE from HLA sequence and locus

        :param locus: string containing HLA locus.
        :param sequence: string containing sequence data.

        :return: GFEobject.
        """
        [loc, accessions] = gfe.split("w")
        if loc in self.loci:
            return self.loci[loc]
        return

    def get_gfe_call(self, hla):
        """
        creates GFE from HLA sequence and locus

        :param locus: string containing HLA locus.
        :param sequence: string containing sequence data.

        :return: GFEobject.
        """
        gfe = self.get_gfe(hla)
        gfe_c = GfeCall(hla=gfe, act_version=self.version,
                        gfedb_version='0.0.2')
        return gfe_c

    def get_gfe(self, hla):
        """
        creates GFE from HLA sequence and locus

        :param locus: string containing HLA locus.
        :param sequence: string containing sequence data.

        :return: GFEobject.
        """
        gfe_query = hla_gfe(hla)
        gfe_data = pa.DataFrame(self.graph.data(gfe_query))
        if not gfe_data.empty:
            gfe = [x for x in gfe_data["GFE"]]
            return gfe
        else:
            return

    def get_hla(self, gfe):
        """
        creates GFE from HLA sequence and locus

        :param locus: string containing HLA locus.
        :param sequence: string containing sequence data.

        :return: GFEobject.
        """
        hla_query = gfe_hla(gfe)
        hla_data = pa.DataFrame(self.graph.data(hla_query))
        if not len(hla_data["HLA"][0]) == 0:
            hla_alleles = flatten(hla_data["HLA"])
            return hla_alleles
        else:
            return

    def sequence(self, seq_type, allele):
        """
        creates GFE from HLA sequence and locus

        :param locus: string containing HLA locus.
        :param sequence: string containing sequence data.

        :return: GFEobject.
        """
        seq_query = get_sequence(seq_type, allele)
        seq_data = pa.DataFrame(self.graph.data(seq_query))
        if not seq_data.empty:
            seq = [x for x in seq_data["SEQ"]]
            return seq
        else:
            return

    def get_gfe_features(self, hlas, features):

        gfe_features = {}
        for hla in hlas:
            gfes = self.get_gfe(hla)
            for gfe in gfes:
                gfe_parts = self.breakup_gfe(gfe)
                for feat in features:
                    f = lc(feat)
                    if not f in gfe_features:
                        gfe_features.update({f: [gfe_parts[f]]})
                    else:
                        fa = gfe_features[f]
                        if not gfe_parts[f] in fa:
                            fa.append(gfe_parts[f])
                            gfe_features[f] = fa
        return gfe_features

    def search_features(self, hlas, features):
        """
        creates GFE from HLA sequence and locus

        :param locus: string containing HLA locus.
        :param sequence: string containing sequence data.

        :return: GFEobject.
        """
        [locus, allele] = hlas[0].split("*")
        search_feats = self.get_gfe_features(hlas, features)
        feat_query = search_hla_features("HLA-A", search_feats)
        feat_data = pa.DataFrame(self.graph.data(feat_query))
        fc = FeatureCall(alleles=hlas, features_searched=features,
                         act_version=self.version, gfedb_version='0.0.2')
        found_hla = {}
        if not feat_data.empty:
            for i in range(0, len(feat_data['HLA'])):
                hla = feat_data['HLA'][i]
                gfe = feat_data['GFE'][i]
                if hla not in found_hla:
                    hla_typing = Typing(hla=hla,
                                        related_gfe=[GfeTyping(gfe=gfe)])
                    found_hla.update({hla: hla_typing})
                else:
                    hla_typing = found_hla[hla]
                    hla_typing.related_gfe.append(GfeTyping(gfe=gfe))
                    found_hla.update({hla: hla_typing})
            fc.matched = list(found_hla.values())
            return fc
        else:
            return fc



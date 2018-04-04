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
import pandas as pd
from py2neo import Graph
#import swagger_client
# from swagger_client.rest import ApiException
# from swagger_client.api_client import ApiClient
import os
import glob
import re
import json

from pygfe.util import get_structures

flatten = lambda l: [item for sublist in l for item in sublist]
is_gfe = lambda x: True if re.search("\d+-\d+-\d+", x) else False
is_kir = lambda x: True if re.search("KIR", x) else False
is_classII = lambda x: True if re.search("HLA-D", x) else False
is_classI = lambda x: True if re.search("HLA-\Dw", x) else False
lc = lambda x: x.lower() if not re.search("UTR", x) else x.lower().replace("utr", "UTR")


# gfedb = GFEDB(url=url, pass=)
# act = ACT(gfedb=gfedb, gfe=GFE, seqann=seqann)
# seqsearch = SeqSearch(gfedb=gfedb)

# TODO: Add version to GFE DB that can be retrieved
#       in a query
# TODO: Change to model
class GfeDB(object):
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
    def __init__(self, graph: Graph=None,
                 persist: bool=False,
                 verbose: bool=False):
        '''
        Constructor
        '''
        self.persist = persist
        self.graph = graph

        # TODO: pull version from graph
        self.version = ''
        self.structures = get_structures()

        # max_hlaid = pd.DataFrame(self.graph.data('MATCH(hla:IMGT_HLA) RETURN max(hla.alleleId) AS ID'))
        # max_gfeid = pd.DataFrame(self.graph.data('MATCH(gfe:GFE) RETURN max(gfe.alleleId) AS ID'))
        # max_fullseqid = pd.DataFrame(self.graph.data('MATCH(seq:SEQUENCE) RETURN max(seq.sequenceId) AS ID'))
        # max_sequenceid = pd.DataFrame(self.graph.data('MATCH(feat:FEATURE) RETURN max(feat.sequenceId) AS ID'))

        #self.max_alleleid = max(max_hlaid['ID'][0], max_gfeid['ID'][0])
        #self.max_seqid = max(max_fullseqid['ID'][0], max_sequenceid['ID'][0])

        # TODO: change out
        self.admin = ''
        if os.getenv("NEO4JADMIN"):
            self.admin = os.getenv("NEO4JADMIN")

    def get_seqid(self, sequence, seqtype, rank):
        """
        creates GFE from HLA sequence and locus

        :param locus: string containing HLA locus.
        :param sequence: string containing sequence data.

        :return: GFEobject.
        """
        if seqtype == "SEQUENCE":
            fullseq_q = fullseqid(sequence)
            fullid = pd.DataFrame(self.graph.data(fullseq_q))
            if not fullid.empty:
                return fullid['ID'][0]
            else:
                self.max_seqid += 1
                return self.max_seqid
        else:
            seq_q = seqid(sequence, seqtype, rank)
            seq_id = pd.DataFrame(self.graph.data(seq_q))
            if not seq_id.empty:
                return seq_id['ID'][0]
            else:
                self.max_seqid += 1
                return self.max_seqid

    def get_alleleid(self, typing: str) -> int:
        """
        creates GFE from HLA sequence and locus

        :param locus: string containing HLA locus.
        :param sequence: string containing sequence data.

        :return: GFEobject.
        """
        if is_gfe(typing):
            gfe_q = gfe_alleleid(typing)
            gfeid = pd.DataFrame(self.graph.data(gfe_q))
            if not gfeid.empty:
                return gfeid['ID'][0]
            else:
                self.max_alleleid += 1
                return self.max_alleleid
        else:
            hla_q = hla_alleleid(typing)
            hlaid = pd.DataFrame(self.graph.data(hla_q))
            if not hlaid.empty:
                return hlaid['ID'][0]
            else:
                self.max_alleleid += 1
                return self.max_alleleid

    def persist_typing(self, typing: str):
        """
        creates GFE from HLA sequence and locus

        :param locus: string containing HLA locus.
        :param sequence: string containing sequence data.

        :return: GFEobject.
        """
        [loc, accesions] = typing.gfe.split("w")

        gfe_alleleid = self.get_alleleid(typing.gfe)
        gfe_node = Node("GFE", name=str(typing.gfe),
                        imgtdb=str(typing.imgtdb_version),
                        locus=str(loc),
                        alleleId=int(gfe_alleleid),
                        status=str("persisted"))

        full_seq = typing.full_gene.sequence
        full_seqid = self.get_seqid(full_seq, "SEQUENCE", 0)
        seq_node = Node("SEQUENCE", name="Sequence", rank=int(0),
                        sequenceId=int(full_seqid),
                        sequence=str(full_seq),
                        length=len(full_seq),
                        nuc=[s for s in list(full_seq)],
                        status=str("persisted"))

        hla_nodes = []
        relationships = []
        for hla_typing in typing.typing:
            hla_alleleid = self.get_alleleid(hla_typing.hla)
            hla_node = Node("IMGT", name=hla_typing.hla,
                            imgtdb=typing.imgtdb_version, locus=loc,
                            alleleId=int(hla_alleleid))
            gfe_rel = Relationship(hla_node, "HAS_GFE", gfe_node,
                                   status="persisted")
            seqgfe_rel = Relationship(hla_node, "HAS_FEATURE", seq_node,
                                      status="persisted")
            seqhla_rel = Relationship(gfe_node, "HAS_FEATURE", seq_node,
                                      status="persisted")
            hla_nodes.append(hla_node)
            relationships.append(gfe_rel)
            relationships.append(seqgfe_rel)
            relationships.append(seqhla_rel)

        for feat in typing.typing_status.novel_features:
            feat_seqid = self.get_seqid(feat.sequence,
                                        feat.term.upper(), feat.rank)
            feature_node = Node("FEATURE", name=str(feat.term.upper()),
                                rank=str(feat.rank),
                                sequenceId=int(feat_seqid),
                                sequence=str(feat.sequence),
                                length=len(feat.sequence),
                                nuc=[s for s in list(feat.sequence)],
                                status=str("persisted"))
            gfefeat_rel = Relationship(gfe_node, "HAS_FEATURE",
                                       feature_node, accession=str(feat.accession),
                                       status="persisted")
            relationships.append(gfefeat_rel)
            for hla_n in hla_nodes:
                hlafeat_rel = Relationship(hla_n, "HAS_FEATURE",
                                           feature_node,
                                           accession=str(feat.accession),
                                           status="persisted")
                relationships.append(hlafeat_rel)

        for feat in typing.features:
            if not feat in typing.typing_status.novel_features:
                feat_seqid = self.get_seqid(feat.sequence,
                                            feat.term.upper(), feat.rank)
                feature_node = Node("FEATURE", name=str(feat.term.upper()),
                                    rank=str(feat.rank),
                                    sequenceId=int(feat_seqid),
                                    sequence=str(feat.sequence),
                                    length=str(len(feat.sequence)),
                                    nuc=[s for s in list(feat.sequence)])
                gfefeat_rel = Relationship(gfe_node, "HAS_FEATURE",
                                           feature_node, accession=str(feat.accession),
                                           status="persisted")
                relationships.append(gfefeat_rel)
                for hla_n in hla_nodes:
                    hlafeat_rel = Relationship(hla_n, "HAS_FEATURE",
                                               feature_node,
                                               accession=str(feat.accession),
                                               status="persisted")
                    relationships.append(hlafeat_rel)

        tx = self.graph.begin()
        for rel in relationships:
            tx.merge(rel)
        tx.commit()

    def get_persisted(self):
        """
        creates GFE from HLA sequence and locus

        :param locus: string containing HLA locus.
        :param sequence: string containing sequence data.

        :return: GFEobject.
        """
        # TODO: Change model of persisted
        persisted = Persisted(act_version=self.version, gfedb_version='0.0.2')
        per_data = pd.DataFrame(self.graph.data(persisted_query()))
        persisted_a = []
        if not per_data.empty:
            for i in range(0, len(per_data['HLA'])):
                per = PersistedData(hla=per_data['HLA'][i], gfe=per_data['GFE'][i],
                                    term=per_data['TERM'][i], rank=per_data['RANK'][i],
                                    accession=per_data['ACCESSION'][i],
                                    sequence=per_data['SEQUENCE'][i])
                persisted_a.append(per)

        persisted.persisted_data = persisted_a
        return persisted




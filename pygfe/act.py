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
from pygfe.models.gfe_call import GfeCall
from pygfe.models.gfe_typing import GfeTyping
from pygfe.models.allele_call import AlleleCall
from pygfe.models.feature_call import FeatureCall
from pygfe.models.typing_status import TypingStatus
from pygfe.models.ars_call import ArsCall
from pygfe.models.persisted import Persisted
from pygfe.models.persisted_data import PersistedData

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

from pygfe.gfe import GFE
from pygfe.gfedb import GfeDB

flatten = lambda l: [item for sublist in l for item in sublist]
is_gfe = lambda x: True if re.search("\d+-\d+-\d+", x) else False
is_kir = lambda x: True if re.search("KIR", x) else False
is_classII = lambda x: True if re.search("HLA-D", x) else False
is_classI = lambda x: True if re.search("HLA-\Dw", x) else False
lc = lambda x: x.lower() if not re.search("UTR", x) else x.lower().replace("utr", "UTR")


class ACT(object):
    '''
    classdocs
    '''

    def __init__(self, gfedb: GfeDB=None,
                 seqann: BioSeqAnn=None,
                 gfe: GFE=None,
                 verbose: bool=False):
        '''
        Constructor
        '''
        self.gfe = gfe
        self.seqann = seqann
        self.gfedb = gfedb

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
            seq_features = pa.DataFrame(self.gfedb.graph.data(feat_q))
            if seq_features.empty:
                unique.append(feat)
        return unique

    def type_from_seq(self, locus, sequence):
        """
        creates GFE from HLA sequence and locus

        :param locus: string containing HLA locus.
        :param sequence: string containing sequence data.

        :return: GFEobject.
        """
        ac_object = AlleleCall()

        # TODO: ADD VERSION
        ac_object.act_version = ''
        ac_object.typing_status = TypingStatus()
        sequence = sequence.upper()
        sequence_typing = self.sequence_lookup(locus, sequence)
        if sequence_typing:
            ac_object.typing_status.status = "documented"
            ac_object.typing = [sequence_typing[0]]
            ac_object.gfe = sequence_typing[1]
            ac_object.features = sequence_typing[2]
            ac_object.full_gene = Feature(rank="1", sequence=sequence, term="gene")
            return ac_object
        else:
            gfe_o = self.gfe_create(locus, sequence)
            ac_object.gfe = gfe_o['gfe']

            ac_object.features = [Feature(accession=f.accession, rank=f.rank, sequence=f.sequence, term=f.term) for f in gfe_o['structure']]
            ac_object.typing_status.novel_features = self.unique_features(ac_object.features)
            related_gfe = self.gfe_lookup(ac_object.gfe, ac_object.features)

            if(len(ac_object.typing_status.novel_features) != 0):
                ac_object.typing_status.status = "novel"
            else:
                ac_object.typing_status.status = "novel_combination"

            if related_gfe:
                ac_object.typing = related_gfe
            else:
                ac_object.typing = self.find_similar(ac_object.gfe, ac_object.features)

            return ac_object

    def type_from_gfe(self, type_gfe):
        """
        creates GFE from HLA sequence and locus

        :param locus: string containing HLA locus.
        :param sequence: string containing sequence data.

        :return: GFEobject.
        """
        ac_object = AlleleCall()
        ac_object.gfe = type_gfe
        ac_object.act_version = self.version
        ac_object.typing_status = TypingStatus()
        seq_features = pa.DataFrame(self.gfedb.graph.data(get_features(type_gfe)))
        if not seq_features.empty:
            ac_object.typing_status.status = "documented"
            features = list()
            for i in range(0, len(seq_features['term'])):
                feature = Feature(accession=seq_features['accession'][i], rank=seq_features['rank'][i], sequence=seq_features['sequence'][i], term=lc(seq_features['term'][i]))
                features.append(feature)
            seq_o = self.gfe.get_sequence(type_gfe)
            ac_object.full_gene = Feature(rank="1", sequence=seq_o.sequence, term="gene")
            ac_object.features = features
            related_gfe = self.gfe_lookup(type_gfe, ac_object.features)
            ac_object.typing = related_gfe
        else:
            seq_o = self.gfe.get_sequence(type_gfe)
            ac_object.full_gene = Feature(rank="1", sequence=seq_o.sequence, term="gene")
            ac_object.features = [Feature(accession=f.accession, rank=f.rank, sequence=f.sequence, term=lc(f.term)) for f in seq_o.structure]
            ac_object.typing_status.novel_features = self.unique_features(ac_object.features)
            ac_object.typing = self.find_similar(ac_object.gfe, ac_object.features)

            if(len(ac_object.typing_status.novel_features) != 0):
                ac_object.typing_status.status = "novel"
            else:
                ac_object.typing_status.status = "novel_combination"

            if self.admin == self.user and self.persist:
                self.persist_typing(ac_object)
        return ac_object

    def sequence_lookup(self, locus, sequence):
        """
        Looks up sequence from

        :param locus: string containing HLA locus.
        :param sequence: string containing sequence data.

        :return: GFEobject.
        """
        lookup_query = sequence_search(locus, sequence)
        sequence_data = pa.DataFrame(self.gfedb.graph.data(lookup_query))
        if not sequence_data.empty:
            features = list()
            gfe = list(set([x for x in sequence_data["GFE"]]))
            hla = list(set([x for x in sequence_data["HLA"]]))
            seq_features = pa.DataFrame(self.gfedb.graph.data(get_features(gfe[0])))
            for i in range(0, len(seq_features['term'])):
                feature = Feature(accession=seq_features['accession'][i], rank=seq_features['rank'][i], sequence=seq_features['sequence'][i], term=lc(seq_features['term'][i]))
                features.append(feature)

            typing = Typing(hla=hla[0], related_gfe=[GfeTyping(gfe=gfe[0], shares=features, features_shared=len(features))])
            return [typing, gfe[0], features]
        else:
            return

    def gfe_create(self, locus, sequence):
        """
        creates GFE from HLA sequence and locus

        :param locus: string containing HLA locus.
        :param sequence: string containing sequence data.

        :return: GFEobject.
        """
        seq_rec = SeqRecord(Seq(sequence, IUPAC.unambiguous_dna), id="GFE")
        annotation = self.seqann.annotate(seq_rec, locus)
        features, gfe = self.gfe.get_gfe(annotation, locus)
        return {'gfe': gfe, 'structure': features}

    def gfe_lookup(self, gfe, features):
        """
        creates GFE from HLA sequence and locus

        :param locus: string containing HLA locus.
        :param sequence: string containing sequence data.

        :return: GFEobject.
        """
        gfe_data = pa.DataFrame(self.gfedb.graph.data(gfe_search(gfe)))
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

    def find_similar(self, gfe, features):
        """
        creates GFE from HLA sequence and locus

        :param locus: string containing HLA locus.
        :param sequence: string containing sequence data.

        :return: GFEobject.
        """
        if is_classI(gfe):
            gfe_dict = self.breakup_gfe(gfe)
            [locus, feature_accessions] = gfe.split("w")
            cypher = similar_gfe_classI(gfe, gfe_dict["EXON-2"],
                                        gfe_dict["EXON-3"])
            similar_data = pa.DataFrame(self.gfedb.graph.data(cypher))
            return self.create_typing(similar_data, gfe, features)
        elif is_classII(gfe):
            gfe_dict = self.breakup_gfe(gfe)
            [locus, feature_accessions] = gfe.split("w")
            cypher = similar_gfe_classII(gfe, gfe_dict["EXON-2"])
            similar_data = pa.DataFrame(self.gfedb.graph.data(cypher))
            return self.create_typing(similar_data, gfe, features)
        elif is_kir(gfe):
            return self.find_gfe_kir(gfe, features)
        else:
            return

    def create_typing(self, similar_data, gfe, features):
        """
        creates GFE from HLA sequence and locus

        :param locus: string containing HLA locus.
        :param sequence: string containing sequence data.

        :return: GFEobject.
        """
        if not similar_data.empty:

            gfe_sim = {sim_gfe: self.calcDiff(sim_gfe, gfe)
                       for sim_gfe in similar_data['GFE']}

            hla_sim = {similar_data['GFE'][i]: similar_data['HLA'][i]
                       for i in range(0, len(similar_data['GFE']))}

            found_hla = {}
            max_val = max(gfe_sim.values())
            max_gfes = [[g, hla_sim[g]] for g in gfe_sim.keys() if gfe_sim[g] == max_val]

            for gfes_hla in max_gfes:
                gfes, hla = gfes_hla[0], gfes_hla[1]
                if hla not in found_hla:
                    matched_features = self.matching_features(gfe, gfes, self.map_structures(features))
                    hla_typing = Typing(hla=hla, related_gfe=[GfeTyping(gfe=gfes, features_shared=len(matched_features))])
                    found_hla.update({hla: hla_typing})
                else:
                    hla_typing = found_hla[hla]
                    matched_features = self.matching_features(gfe, gfes, self.map_structures(features))
                    hla_typing.related_gfe.append(GfeTyping(gfe=gfes, features_shared=len(matched_features)))
                    found_hla.update({hla: hla_typing})

            typing_list = list(found_hla.values())
            return typing_list
        else:
            return list()

    def find_gfe_kir(self, gfe, features):
        """
        creates GFE from HLA sequence and locus

        :param locus: string containing HLA locus.
        :param sequence: string containing sequence data.

        :return: GFEobject.
        """
        [locus, feature_accessions] = gfe.split("w")
        cypher = similar_kir(locus)
        similar_data = pa.DataFrame(self.gfedb.graph.data(cypher))
        return self.create_typing(similar_data, gfe, features)

    def gfe2hla(self, gfe):
        """
        creates GFE from HLA sequence and locus

        :param locus: string containing HLA locus.
        :param sequence: string containing sequence data.

        :return: GFEobject.
        """
        cypher = gfe_hla(gfe)
        hlacypher = pa.DataFrame(self.gfedb.graph.data(cypher))
        hla_alleles = set()

        if not len(hlacypher["HLA"][0]) == 0:
            hla_alleles = flatten(hlacypher["HLA"])
            return hla_alleles
        else:
            return

    def map_structures(self, gfe_structs):
        """
        creates GFE from HLA sequence and locus

        :param locus: string containing HLA locus.
        :param sequence: string containing sequence data.

        :return: GFEobject.
        """
        feat_map = {}
        for feat in gfe_structs:
            feat_name = "-".join([feat.term.upper(), str(feat.rank)])
            feat_map.update({feat_name: feat.sequence})

        return feat_map

    def matching_features(self, gfe1, gfe2, structures):
        """
        creates GFE from HLA sequence and locus

        :param locus: string containing HLA locus.
        :param sequence: string containing sequence data.

        :return: GFEobject.
        """

        # TODO: Update looping
        gfe_parts1 = self.breakup_gfe(gfe1)
        gfe_parts2 = self.breakup_gfe(gfe2)
        feat_list = list()
        for feat in gfe_parts1:
            if feat in gfe_parts2:
                if gfe_parts1[feat] == gfe_parts2[feat]:
                    feat_term, feat_rank = feat.split('-')
                    shared_feat = Feature(term=feat_term.upper(), rank=feat_rank, sequence=structures[feat], accession=gfe_parts1[feat])
                    feat_list.append(shared_feat)

        return(feat_list)

    def breakup_gfe(self, gfe):
        """
        creates GFE from HLA sequence and locus

        :param locus: string containing HLA locus.
        :param sequence: string containing sequence data.

        :return: GFEobject.
        """
        [locus, feature_accessions] = gfe.split("w")
        accessions = feature_accessions.split("-")
        i = 0
        features = {}
        for feature_rank in self.gfedb.structures[locus]:
            accession = accessions[i]
            features.update({feature_rank: accession})
            i += 1
        return(features)

    def calcDiff(self, gfe1, gfe2):
        """
        creates GFE from HLA sequence and locus

        :param locus: string containing HLA locus.
        :param sequence: string containing sequence data.

        :return: GFEobject.
        """
        count = 0
        gfe1_parts = gfe1.split("-")
        gfe2_parts = gfe2.split("-")
        for i in range(0, len(gfe1_parts)):
            if gfe1_parts[i] == gfe2_parts[i]:
                count += 1

        return count

    # TODO: MERGE with get_gfe
    def get_gfe_call(self, hla):
        """
        creates GFE from HLA sequence and locus

        :param locus: string containing HLA locus.
        :param sequence: string containing sequence data.

        :return: GFEobject.
        """
        gfe = self.get_gfe(hla)
        gfe_c = GfeCall(hla=gfe,
                        act_version=self.version,
                        gfedb_version='0.0.2')
        return gfe_c

    # TODO: MERGE with get_gfe_call
    def get_gfe(self, hla):
        """
        creates GFE from HLA sequence and locus

        :param locus: string containing HLA locus.
        :param sequence: string containing sequence data.

        :return: GFEobject.
        """
        gfe_query = hla_gfe(hla)
        gfe_data = pa.DataFrame(self.gfedb.graph.data(gfe_query))
        if not gfe_data.empty:
            gfe = [x for x in gfe_data["GFE"]]
            return gfe
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
        seq_data = pa.DataFrame(self.gfedb.graph.data(seq_query))
        if not seq_data.empty:
            seq = [x for x in seq_data["SEQ"]]
            return seq
        else:
            return

    def typing_to_bioseq(self, typing):
        """
        creates GFE from HLA sequence and locus

        :param locus: string containing HLA locus.
        :param sequence: string containing sequence data.

        :return: GFEobject.
        """

        # TODO: Add more annotation and qualifiers
        # TODO: Use the full sequence accession as the ID
        sequence = typing.full_gene.sequence
        seqid = 'GFE'
        comments = []
        description = "GFE " + typing.gfe
        comments.append("Typing Status: " + typing.typing_status.status)
        if hasattr(typing.typing_status, 'novel_features'):
            uniq = " ".join(["-".join([feat.term, str(feat.rank)]) for feat in typing.typing_status.novel_features])
            comments.append("Novel features: " + uniq)
        
        allele_typed = "/".join([typ.hla for typ in typing.typing])
        comments.append("Allele Call: " + allele_typed)
        comments.append("")
        comments.append("Typed with ACT Service " + self.version)
        if hasattr(typing, 'gfe_version'):
            comments.append("Annotated with GFE Service " + typing.gfe_version)

        if hasattr(typing.full_gene, 'accession'):
            seqid = 'GFEw' + str(typing.full_gene.accession)
        seqrecord = SeqRecord(Seq(sequence, IUPAC.unambiguous_dna), id=seqid, description=description)
        source_feature = SeqFeature(FeatureLocation(0, len(str(seqrecord.seq))), type="source", strand=1)
        seqrecord.annotations["sequence_version"] = 1
        seqrecord.annotations["molecule_type"] = "DNA"
        seqrecord.annotations["data_file_division"] = "HUM"
        if hasattr(typing.full_gene, 'accession'):
            seqrecord.annotations["accessions"] = [seqid]

        seqrecord.annotations["comment"] = comments
        seqrecord.features.append(source_feature)
        seqrecord.features[0].qualifiers = OrderedDict([('organism', ['Homo sapiens']), ('mol_type', ['genomic DNA']), ('db_xref', ['taxon:9606'])])

        for feat in typing.features:
            start_pos = sequence.find(feat.sequence)
            if len(seqrecord.features) == 1:
                start_pos = 0

            feat_type = feat.term
            if feat_type == "five_prime_UTR" or feat_type == "three_prime_UTR":
                feat_type = "UTR"
            end_pos = start_pos + len(feat.sequence)
            seq_feature = SeqFeature(FeatureLocation(start_pos, end_pos), type=feat_type, strand=1)
            
            if feat.term == 'exon' or feat.term == 'intron':
                seq_feature.qualifiers = OrderedDict([('number', [feat.rank])])
            seqrecord.features.append(seq_feature)

        return seqrecord




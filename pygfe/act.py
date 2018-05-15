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
#from pygfe.models.feature import Feature
from pygfe.models.typing import Typing
from pygfe.models.feature import Feature
from pygfe.models.seqdiff import Seqdiff

from py2neo import Node, Relationship
import pandas as pa

import time
import os
import glob
import re
import json

from typing import Dict
import seqann

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
from pandas import DataFrame

import pygfe
from pygfe.gfe import GFE
from pygfe.gfedb import GfeDB

from pygfe.cypher import all_gfe2hla
from pygfe.cypher import all_seq2hla
from pygfe.cypher import all_gfe2feats
from pygfe.cypher import hla_seqdiff

import logging

flatten = lambda l: [item for sublist in l for item in sublist]
is_gfe = lambda x: True if re.search("\d+-\d+-\d+", x) else False
is_kir = lambda x: True if re.search("KIR", x) else False
is_classII = lambda x: True if re.search("HLA-D", x) else False
is_classI = lambda x: True if re.search("HLA-\Dw", x) else False
lc = lambda x: x.lower() if not re.search("UTR", x) else x.lower().replace("utr", "UTR")
isutr = lambda f: True if re.search("UTR", f) else False

logging.basicConfig(format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
                    datefmt='%m/%d/%Y %I:%M:%S %p',
                    level=logging.INFO)


class ACT(object):
    '''
    classdocs
    '''
    def __init__(self, gfedb: GfeDB=None,
                 seqann: BioSeqAnn=None,
                 gfe: GFE=None,
                 features: Dict=None,
                 gfe2hla: Dict=None,
                 gfe_feats: DataFrame=None,
                 seq2hla: DataFrame=None,
                 load_gfe2hla: bool=False,
                 load_seq2hla: bool=False,
                 load_gfe2feat: bool=False,
                 pid: str="NA",
                 verbose: bool=False):
        '''
        Constructor
        '''
        self.gfe_feats = gfe_feats
        self.gfe = gfe
        self.seqann = seqann
        self.gfedb = gfedb
        self.features = features
        self.gfe2hla = gfe2hla
        self.seq2hla = seq2hla
        self.logger = logging.getLogger("Logger." + __name__)
        self.logname = "ID {:<10} - ".format(str(pid))

        if load_gfe2feat:
            self.gfe_feats = pa.DataFrame(self.gfedb.graph.data(all_gfe2feats()))

        if load_seq2hla:
            self.seq2hla = pa.DataFrame(self.gfedb.graph.data(all_seq2hla()))

        if load_gfe2hla:
            tmp_gfe = {}
            gfehla_df = pa.DataFrame(self.gfedb.graph.data(all_gfe2hla()))
            for loc in gfehla_df['LOC'].unique().tolist():
                if re.search("HLA-\D$", loc):
                    loc_df = gfehla_df.loc[gfehla_df['LOC'] == loc]
                    loc1 = self.gfe.structures[loc]['exon-2']
                    loc2 = self.gfe.structures[loc]['exon-3']
                    loc_df['EXON23'] = loc_df['GFE'].apply(lambda gfe: "-".join([gfe.split("-")[loc1],gfe.split("-")[loc2]]))
                    tmp_gfe.update({loc: loc_df})

                if is_classII(loc):
                    loc_df = gfehla_df.loc[gfehla_df['LOC'] == loc]
                    loc1 = self.gfe.structures[loc]['exon-2']
                    loc_df['EXON2'] = loc_df['GFE'].apply(lambda gfe:gfe.split("-")[loc1])
                    tmp_gfe.update({loc: loc_df})
            self.gfe2hla = tmp_gfe

    def hla_seqdiff(self, locus, imgtdb_version, typ1, typ2):
        diff_data = pa.DataFrame(self.gfedb.graph.data(hla_seqdiff(locus, imgtdb_version, typ1, typ2)))
        diffs = []
        for i in range(0, len(diff_data)):
            index = diff_data['POS'][i]
            diff_feat = self.seqann.refdata.align_coordinates[locus.split("-")[1]][index]
            ref_seq = diff_data['VAR2'][i]
            in_seq = diff_data['VAR1'][i]
            relative_location = index + self.seqann.refdata.location[locus]
            diff = Seqdiff(term=diff_feat, location=str(relative_location), ref=ref_seq, inseq=in_seq)
            diffs.append(diff)
        return diffs

    # TODO: Make this come from the
    #       pygfe cache instead of making
    #       a bunch of calls to the GFE DB
    def unique_features(self, features, locus, imgtdb_version):
        """
        creates GFE from HLA sequence and locus

        :param locus: string containing HLA locus.
        :param sequence: string containing sequence data.
        :return: GFEobject.
        """
        unique = []
        if self.features:
            for feat in features:
                rank = feat.rank if not isutr(feat.term) else 1
                feat_str = ":".join([imgtdb_version, locus, str(rank), feat.term.upper(), feat.sequence])
                if not feat_str in self.features:
                    unique.append(feat)
        else:
            for feat in features:
                feat_q = search_feature(feat.term.upper(),
                                        feat.rank,
                                        feat.sequence)
                seq_features = pa.DataFrame(self.gfedb.graph.data(feat_q))
                if seq_features.empty:
                    unique.append(feat)
        return unique

    def type_from_seq(self, locus, sequence, imgtdb_version):
        """
        creates GFE from HLA sequence and locus

        :param locus: string containing HLA locus.
        :param sequence: string containing sequence data.

        :return: GFEobject.
        """
        # TODO: Add full gene accession
        ac_object = Typing()
        ac_object.imgtdb_version = "".join(imgtdb_version.split("."))
        ac_object.pygfe_version = pygfe.__version__
        ac_object.gfedb_version = '0.0.2'
        ac_object.seqann_version = seqann.__version__
        sequence = sequence.upper()
        sequence_typing = self.sequence_lookup(locus, sequence, ac_object.imgtdb_version)
        if sequence_typing:
            ac_object.status = "documented"
            ac_object.hla = sequence_typing[0]
            ac_object.gfe = sequence_typing[1]
            ac_object.closest_gfe = sequence_typing[1]
            ac_object.features = sequence_typing[2]

            if self.verbose:
                self.logger.info(self.logname + locus + " sequence documented for " + imgtdb_version + " | " + ac_object.gfe + " = " + ac_object.hla)

            return ac_object
        else:
            # time GFE creation
            time_start = time.time()
            gfe_o = self.gfe_create(locus, sequence)
            if not 'annotation' in gfe_o:
                self.logger.error(self.logname + "Failed to create annotation!!")
                ac_object.gfe = "NA"
                ac_object.status = "Failed_Annotation"
                ac_object.features = []
                ac_object.closest_gfe = "NA"
                ac_object.hla = "NA"
                return ac_object

            if self.verbose:
                time_taken = int(time.time() - time_start)
                self.logger.info(self.logname + " gfe_create time for " + locus + " "
                             + imgtdb_version + " = " + str(time_taken) + " minutes")
            annotation = gfe_o['annotation']
            ac_object.gfe = gfe_o['gfe']
            ac_object.features = [Feature(accession=f.accession,
                                          rank=f.rank,
                                          sequence=f.sequence,
                                          term=f.term) for f in gfe_o['structure']]
            novel_features = self.unique_features(ac_object.features, locus, ac_object.imgtdb_version)
            if(len(novel_features) != 0):
                if self.verbose:
                    self.logger.info(self.logname + " # novel features = "
                                 + str(len(novel_features)))
                ac_object.novel_features = novel_features
                ac_object.status = "novel"
            else:
                self.logger.info(self.logname + " novel combination")
                ac_object.status = "novel_combination"

            similar_results = self.find_similar(ac_object.gfe,
                                                ac_object.features,
                                                ac_object.imgtdb_version)
            if similar_results:
                ac_object.hla = similar_results[0]
                ac_object.closest_gfe = similar_results[1]
                if self.seqann.align:
                    if self.verbose:
                        self.logger.info(self.logname + " finding sequence differences")
                    ac_object.seqdiff = self.diff_seq(similar_results[0],
                                                      annotation)
                    ac_object.differences = len(ac_object.seqdiff)
            else:
                ac_object.hla = "NA"
                ac_object.closest_gfe = "NA"
                if self.verbose:
                    self.logger.warn(self.logname + " No allele call made!")
            return ac_object

    def rename_feat(self, f):
        if isutr(f):
            return f
        else:
            return "_".join(f.split("-"))

    def diff_seq(self, ref_allele, annotation):
        locus = ref_allele.split("*")[0]
        feat_order = list(self.gfe.struct_order[locus].keys())
        feat_order.sort()
        db = self.seqann.refdata.dbversion
        dbv = ".".join([list(db)[0], "".join(list(db)[1:3]), list(db)[3]])
        aligned_seq = list("".join([annotation.aligned[self.rename_feat(self.gfe.struct_order[locus][i])] for i in feat_order]))
        seq_l = ",".join(["".join(["\"", s, "\""]) for s in aligned_seq])
        query = "MATCH (a:IMGT_HLA)-[r:HAS_ALIGNMENT]-(n:GEN_ALIGN) " \
                + "WHERE a.locus = \"" + locus + "\" " \
                + "AND a.name = \"" + ref_allele + "\" " \
                + "AND r.imgt_release = \"" + dbv + "\" " \
                + "WITH [" + seq_l + "] AS RSEQ,n " \
                + "WITH [x in range(0,size(n.seq)-1) | x] AS ind,RSEQ,n " \
                + "UNWIND ind AS number " \
                + "WITH DISTINCT number,RSEQ,n " \
                + "WHERE NOT(n.seq[number] = RSEQ[number]) "\
                + "RETURN number, n.seq[number] AS REF, RSEQ[number] AS SEQ"
        diff_data = pa.DataFrame(self.gfedb.graph.data(query))
        diffs = []
        for i in range(0, len(diff_data)):
            index = diff_data['number'][i]
            diff_feat = self.seqann.refdata.align_coordinates[locus.split("-")[1]][index]
            ref_seq = diff_data['REF'][i]
            in_seq = diff_data['SEQ'][i]
            relative_location = index + self.seqann.refdata.location[locus]
            # TODO: Change rank!
            diff = Seqdiff(term=diff_feat, location=str(relative_location), ref=ref_seq, inseq=in_seq)
            diffs.append(diff)
        return diffs

    def type_from_gfe(self, type_gfe):
        """
        creates GFE from HLA sequence and locus

        :param locus: string containing HLA locus.
        :param sequence: string containing sequence data.

        :return: GFEobject.
        """
        ac_object = Typing()
        # ac_object.gfe = type_gfe
        # ac_object.act_version = self.version
        # ac_object.typing_status = TypingStatus()
        # seq_features = pa.DataFrame(self.gfedb.graph.data(get_features(type_gfe)))
        # if not seq_features.empty:
        #     ac_object.typing_status.status = "documented"
        #     features = list()
        #     for i in range(0, len(seq_features['term'])):
        #         feature = Feature(accession=seq_features['accession'][i], rank=seq_features['rank'][i], sequence=seq_features['sequence'][i], term=lc(seq_features['term'][i]))
        #         features.append(feature)
        #     seq_o = self.gfe.get_sequence(type_gfe)
        #     ac_object.full_gene = Feature(rank="1", sequence=seq_o.sequence, term="gene")
        #     ac_object.features = features
        #     related_gfe = self.gfe_lookup(type_gfe, ac_object.features)
        #     ac_object.typing = related_gfe
        # else:
        #     seq_o = self.gfe.get_sequence(type_gfe)
        #     ac_object.full_gene = Feature(rank="1", sequence=seq_o.sequence, term="gene")
        #     ac_object.features = [Feature(accession=f.accession, rank=f.rank, sequence=f.sequence, term=lc(f.term)) for f in seq_o.structure]
        #     ac_object.typing_status.novel_features = self.unique_features(ac_object.features)
        #     ac_object.typing = self.find_similar(ac_object.gfe, ac_object.features)

        #     if(len(ac_object.typing_status.novel_features) != 0):
        #         ac_object.typing_status.status = "novel"
        #     else:
        #         ac_object.typing_status.status = "novel_combination"

        #     if self.admin == self.user and self.persist:
        #         self.persist_typing(ac_object)
        return ac_object

    def sequence_lookup(self, locus, sequence, imgtdb_version):
        """
        Looks up sequence from

        :param locus: string containing HLA locus.
        :param sequence: string containing sequence data.

        :return: GFEobject.
        """
        # TODO: add dbversion!
        if not self.seq2hla.empty:
            df = self.seq2hla[(self.seq2hla['DB'] == imgtdb_version) & (self.seq2hla['LOC'] == locus) & (self.seq2hla['SEQ'] == sequence)]
            if not df.empty:
                hla = df['HLA'].tolist()[0]
                gfe = df['GFE'].tolist()[0]
                features = []
                feats = self.gfe_feats[(self.gfe_feats['GFE'] == gfe) &
                                          (self.gfe_feats['DB'] == imgtdb_version)]['FEATS'].tolist()[0]
                if feats and len(feats) > 0:
                    for feat in feats:
                        feature = Feature(accession=feat['accession'],
                                          rank=feat['rank'],
                                          sequence=feat['sequence'],
                                          term=lc(feat['term']))
                        features.append(feature)
                else:
                    seq_features = pa.DataFrame(self.gfedb.graph.data(get_features(gfe)))
                    for i in range(0, len(seq_features['term'])):
                        feature = Feature(accession=seq_features['accession'][i],
                                          rank=seq_features['rank'][i],
                                          sequence=seq_features['sequence'][i],
                                          term=lc(seq_features['term'][i]))
                        features.append(feature)
                return [hla, gfe, features]
        else:
            lookup_query = sequence_search(locus, sequence)
            sequence_data = pa.DataFrame(self.gfedb.graph.data(lookup_query))
            if not sequence_data.empty:
                features = list()
                gfe = list(set([x for x in sequence_data["GFE"]]))
                hla = list(set([x for x in sequence_data["HLA"]]))
                seq_features = pa.DataFrame(self.gfedb.graph.data(get_features(gfe[0])))
                for i in range(0, len(seq_features['term'])):
                    feature = Feature(accession=seq_features['accession'][i],
                                      rank=seq_features['rank'][i],
                                      sequence=seq_features['sequence'][i],
                                      term=lc(seq_features['term'][i]))
                    features.append(feature)
                #typing = Typing(hla=hla[0], related_gfe=[GfeTyping(gfe=gfe[0], shares=features, features_shared=len(features))])
                return [hla[0], gfe[0], features]
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
        if not hasattr(annotation, 'annotation'):
            return annotation
        else:
            features, gfe = self.gfe.get_gfe(annotation, locus)
            return {'gfe': gfe, 'structure': features, 'annotation': annotation}

    def gfe_lookup(self, gfe, features):
        """
        creates GFE from HLA sequence and locus

        :param locus: string containing HLA locus.
        :param sequence: string containing sequence data.

        :return: GFEobject.
        """
        #print("WTF gfe_lookup")
        gfe_data = pa.DataFrame(self.gfedb.graph.data(gfe_search(gfe)))
        if not gfe_data.empty:
            return ["/".join(gfe_data["HLA"]), gfe]
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

    def find_similar(self, gfe, features, imgtdb_version):
        """
        creates GFE from HLA sequence and locus

        :param locus: string containing HLA locus.
        :param sequence: string containing sequence data.

        :return: GFEobject.
        """
        if is_classI(gfe):
            gfe_dict = self.breakup_gfe(gfe)
            [locus, feature_accessions] = gfe.split("w")
            if self.gfe2hla:
                exon23 = "-".join([gfe_dict["EXON-2"], gfe_dict["EXON-3"]])
                df = self.gfe2hla[locus][(self.gfe2hla[locus]['EXON23'] == exon23) & (self.gfe2hla[locus]['DB'] == imgtdb_version)][['GFE', 'HLA']]\
                         .reset_index()
                return self.create_typing(df, gfe, features)
            else:
                cypher = similar_gfe_classI(gfe, gfe_dict["EXON-2"],
                                            gfe_dict["EXON-3"],
                                            imgtdb_version)
                similar_data = pa.DataFrame(self.gfedb.graph.data(cypher))
                return self.create_typing(similar_data, gfe, features)
        elif is_classII(gfe):
            gfe_dict = self.breakup_gfe(gfe)
            [locus, feature_accessions] = gfe.split("w")
            if self.gfe2hla:
                exon2 = gfe_dict["EXON-2"]
                df = self.gfe2hla[locus][(self.gfe2hla[locus]['EXON2'] == exon2) & (self.gfe2hla[locus]['DB'] == imgtdb_version)][['GFE', 'HLA']]\
                         .reset_index()
                return self.create_typing(df, gfe, features)                
            else:
                cypher = similar_gfe_classII(gfe, gfe_dict["EXON-2"],
                                             imgtdb_version)
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
            gfe_sim = {sim_gfe: self.calcSim(sim_gfe, gfe)
                       for sim_gfe in similar_data['GFE']}
            gfe_diff = {sim_gfe: self.calcDiff(sim_gfe, gfe)
                       for sim_gfe in similar_data['GFE']}
            hla_sim = {similar_data['GFE'][i]: similar_data['HLA'][i]
                       for i in range(0, len(similar_data['GFE']))}
            max_val = max(gfe_sim.values())
            min_val = min(gfe_diff.values())

            # TODO: Add more complex logic here about 
            max_gfes = [[g, hla_sim[g]] for g in gfe_sim.keys() if gfe_sim[g] == max_val]
            min_gfes = [[g, hla_sim[g]] for g in gfe_sim.keys() if gfe_diff[g] == min_val]
            diff_max = [g[0] for g in max_gfes if gfe_diff[g[0]] > 0]
            if len(diff_max) == len(max_gfes) and len(min_gfes) == 1 and min_val == 0:
                max_gfes = min_gfes
            hla_list = []
            gfe_list = []
            for gfes_hla in max_gfes:
                gfes, hla = gfes_hla[0], gfes_hla[1]
                hla_list.append(hla)
                gfe_list.append(gfes)
            return ["/".join(hla_list), "/".join(gfe_list)]
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
                    if not feat.upper == "FIVE_PRIME_UTR" and not feat == "THREE_PRIME_UTR":
                        feat_term, feat_rank = feat.split('_')
                        shared_feat = Feature(term=feat_term.upper(), rank=feat_rank, sequence=structures[feat], accession=gfe_parts1[feat])
                        feat_list.append(shared_feat)
                    else:
                        shared_feat = Feature(term=feat, rank=1, sequence=structures["-".join([feat, str(1)])], accession=gfe_parts1["-".join([feat, str(1)])])
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
            if feature_rank.upper() == "FIVE_PRIME_UTR" or feature_rank.upper() == "THREE_PRIME_UTR":
                features.update({"_".join([feature_rank, str(1)]).upper(): accession})
            else:
                features.update({feature_rank.upper(): accession})
            i += 1
        return(features)

    def calcSim(self, gfe1, gfe2):
        """
        creates GFE from HLA sequence and locus

        :param locus: string containing HLA locus.
        :param sequence: string containing sequence data.

        :return: GFEobject.
        """
        count = 0
        gfe1_parts = gfe1.split("w")[1].split("-")
        gfe2_parts = gfe2.split("w")[1].split("-")
        for i in range(0, len(gfe1_parts)):
            if gfe1_parts[i] == gfe2_parts[i] and \
                    not str(gfe1_parts[i]) == str(0):
                count += 1

        return count

    def calcDiff(self, gfe1, gfe2):
        """
        creates GFE from HLA sequence and locus

        :param locus: string containing HLA locus.
        :param sequence: string containing sequence data.

        :return: GFEobject.
        """
        count = 0
        gfe1_parts = gfe1.split("w")[1].split("-")
        gfe2_parts = gfe2.split("w")[1].split("-")
        for i in range(0, len(gfe1_parts)):
            if gfe1_parts[i] != gfe2_parts[i] and \
                    not str(gfe1_parts[i]) == str(0):
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
        #gfe = self.get_gfe(hla)
        # gfe_c = GfeCall(hla=gfe,
        #                 act_version=self.version,
        #                 gfedb_version='0.0.2')
        return ''

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
        comments.append("Typed with GFE Service " + self.version)
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




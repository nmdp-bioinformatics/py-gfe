'''
Created on Feb 8, 2017

@author: mhalagan
'''

import pygfe
from pygfe.cypher import sequence_search, all_db_imgt
from pygfe.cypher import similar_gfe_classI
from pygfe.cypher import similar_gfe_classII
from pygfe.cypher import similar_kir
from pygfe.cypher import get_features
from pygfe.cypher import search_feature
from pygfe.models.feature import Feature
from pygfe.models.typing import Typing
from pygfe.models.error import Error
from pygfe.models.seqdiff import Seqdiff
from pygfe.cypher import all_gfe2hla
from pygfe.cypher import all_seq2hla
from pygfe.cypher import all_gfe2feats
from pygfe.util import get_structures

import seqann
from seqann.blast_cmd import get_locus
from seqann.util import checkseq
from seqann.util import isutr
from seqann.util import is_classII
from seqann.util import is_kir
from seqann import BioSeqAnn

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord

from pandas import DataFrame

from py2neo import Graph

from typing import List
from typing import Dict
from typing import Any

import re
import time
import logging

flatten = lambda l: [item for sublist in l for item in sublist]
is_gfe = lambda x: True if re.search("\d+-\d+-\d+", x) else False
lc = lambda x: x.lower() if not re.search("UTR", x) else x.lower().replace("utr", "UTR")
is_classI = lambda x: True if re.search("HLA-\Dw", x) else False

# Dict of seqann objects

# ** May want to eventually create a GFE service

# TODO: have option for not returning features

# seqann = BioSeqAnn(dbversion="3310")
# gfe = pyGFE(seann={"3310": seqann})
# gfe = pyGFE(seqann=seqann})
# gfe = pyGFE(seqann=[seqann1, seqann2]})


class pyGFE(object):
    '''
    classdocs
    '''
    def __init__(self, url="http://feature.nmdp-bioinformatics.org",
                 loci=['HLA-A', 'HLA-B', 'HLA-C',
                       'HLA-DRB1', 'HLA-DQB1', 'HLA-DRB4',
                       'HLA-DRB5', 'HLA-DPB1', 'HLA-DPA1',
                       'HLA-DQA1', 'HLA-DRB3'],
                 graph: Graph=None,
                 seqann: Any={},
                 features: Dict=None,
                 verbose: bool=False,
                 kir: bool=False,
                 pid: str="NA",
                 gfe2hla: Dict=None,
                 gfe_feats: DataFrame=None,
                 seq2hla: DataFrame=None,
                 load_gfe2hla: bool=False,
                 load_seq2hla: bool=False,
                 load_gfe2feat: bool=False,
                 verbosity=1):
        '''
        Constructor
        '''
        # TODO: Add catch if seqann or graph aren't defined
        self.kir = kir
        self.graph = graph
        self.logger = logging.getLogger("Logger." + __name__)
        if pid:
            self.logname = "ID {:<10} - ".format(str(pid))
        else:
            self.logname = ''

        if not isinstance(seqann, Dict) and seqann:
            if isinstance(seqann, BioSeqAnn):
                self.seqann = {seqann.refdata.dbversion: seqann}
            elif(isinstance(seqann, List)):
                self.seqann = {}
                for ann in seqann:
                    self.seqann.update({ann.refdata.dbversion: ann})
            else:
                raise SeqAnnException(inputtype=type(seqann),
                                      reason="Can't initalize seqann")
        else:
            self.seqann = seqann

        self.features = features
        self.gfe2hla = gfe2hla
        self.seq2hla = seq2hla
        self.gfe_feats = gfe_feats
        self.verbose = verbose
        self.structures = get_structures()

        # ISSUE: gfe_feats & seq2hla need to be loaded together
        #
        if load_gfe2feat:
            self.gfe_feats = self.graph.run(all_gfe2feats()).to_data_frame()
            self.gfe_feats['DBV'] = self.gfe_feats['DB'].apply(lambda db: "".join(db.split(".")))
            self.gfe_feats['DB'] = self.gfe_feats['DBV']
            self.gfe_feats = self.gfe_feats.drop(['DBV'], axis=1)

        if load_seq2hla:
            self.seq2hla = self.graph.run(all_seq2hla()).to_data_frame()
            self.seq2hla['DBV'] = self.seq2hla['DB'].apply(lambda db: "".join(db.split(".")))
            self.seq2hla['DB'] = self.seq2hla['DBV']
            self.seq2hla = self.seq2hla.drop(['DBV'], axis=1)

        if load_gfe2hla:
            tmp_gfe = {}
            gfehla_df = self.graph.run(all_gfe2hla()).to_data_frame()
            for loc in gfehla_df['LOC'].unique().tolist():
                if re.search("HLA-\D$", loc):
                    loc_df = gfehla_df.loc[gfehla_df['LOC'] == loc]
                    loc1 = self.structures[loc]['exon-2']
                    loc2 = self.structures[loc]['exon-3']
                    loc_df['EXON23'] = loc_df['GFE'].apply(lambda gfe: "-".join([gfe.split("-")[loc1],gfe.split("-")[loc2]]))
                    tmp_gfe.update({loc: loc_df})

                if is_classII(loc):
                    loc_df = gfehla_df.loc[gfehla_df['LOC'] == loc]
                    loc1 = self.structures[loc]['exon-2']
                    loc_df['EXON2'] = loc_df['GFE'].apply(lambda gfe:gfe.split("-")[loc1])
                    tmp_gfe.update({loc: loc_df})
            self.gfe2hla = tmp_gfe

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
                feat_str = ":".join([imgtdb_version,
                                     locus, str(rank),
                                     feat.term.upper(), feat.sequence])
                if not feat_str in self.features:
                    unique.append(feat)
        else:
            for feat in features:
                feat_q = search_feature(feat.term.upper(),
                                        feat.rank,
                                        feat.sequence)

                seq_features = self.graph.run(feat_q).to_data_frame()
                if seq_features.empty:
                    unique.append(feat)
        return unique

    # TODO: have parameter for filling in
    #       seqann parameters if it's a new
    #       dbversion
    def type_from_seq(self, locus: str=None,
                      sequence: str=None,
                      imgtdb_version: str="3.31.0",
                      nseqs: int=20,
                      alignseqs: int=10,
                      skip: List=[]):
        """
        creates GFE from HLA sequence and locus

        :param locus: string containing HLA locus.
        :param sequence: string containing sequence data.

        :return: GFEobject.
        """

        # TODO: Add full gene accession
        # TODO: reformt dbversion if missing .
        ac_object = Typing()
        ac_object.imgtdb_version = "".join(imgtdb_version.split("."))
        ac_object.pygfe_version = pygfe.__version__
        ac_object.seqann_version = seqann.__version__
        ac_object.gfedb_version = '0.0.2'

        # If sequence is now a biopython
        # sequence record convert it to one
        if isinstance(sequence, Seq):
            sequence = str(sequence)
        elif(isinstance(sequence, SeqRecord)):
            sequence = str(sequence.seq)

        if not ac_object.imgtdb_version in self.seqann:
            self.seqann.update({ac_object.imgtdb_version:
                                BioSeqAnn(dbversion=ac_object.imgtdb_version,
                                          #store_features=self.store_features,
                                          load_features=self.load_features,
                                          cached_features=self.cached_features)})

        # If sequence contains any characters
        # other than ATCG then the GFE notation
        # can not be created
        valid_seq = checkseq(sequence)

        if self.verbose and not valid_seq:
            self.logger.warning(self.logname + " Sequence alphabet "
                                + "contains non DNA")
            self.logger.warning(self.logname
                                + " No GFE string will be generated")
            raise Exception("Input sequence was not valid! {}".format(sequence))

        # Check it the locus exists
        if not locus:
            if self.verbose:
                self.logger.info(self.logname + " No locus provided! ")

            # Guessing locus with blastn
            locus = get_locus(sequence,
                              kir=self.kir,
                              refdata=self.seqann[ac_object.imgtdb_version].refdata)

            if locus and self.verbose:
                self.logger.info(self.logname + " Locus prediction = " + locus)

            if not locus:
                if self.verbose:
                    self.logger.error(self.logname
                                      + " Locus could not be determined!")
                # TODO: Raise exception
                raise Exception("Locus could not be determined! {}".format(sequence))

        sequence = sequence.upper()
        sequence_typing = self.sequence_lookup(locus,
                                               sequence,
                                               ac_object.imgtdb_version)
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
            gfe_o = self.gfe_create(locus, sequence, ac_object.imgtdb_version)
            if not 'annotation' in gfe_o:
                self.logger.error(self.logname + "Failed to create annotation!!")
                error = Error("Failed to create annotation!!", ac_object.pygfe_version, ac_object.gfedb_version, imgtdb_version)
                return error

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
            novel_features = self.unique_features(ac_object.features,
                                                  locus,
                                                  ac_object.imgtdb_version)
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
                                                imgtdb_version)
            if similar_results:
                ac_object.hla = similar_results[0]
                ac_object.closest_gfe = similar_results[1]
                if self.seqann[ac_object.imgtdb_version].align:
                    if self.verbose:
                        self.logger.info(self.logname + " finding sequence differences")
                    ac_object.seqdiff = self.diff_seq(similar_results[0],
                                                      annotation,
                                                      imgtdb_version)
                    ac_object.differences = len(ac_object.seqdiff)
            else:
                ac_object.hla = "NA"
                ac_object.closest_gfe = "NA"
                if self.verbose:
                    self.logger.warn(self.logname + " No allele call made!")
            return ac_object

    def sequence_lookup(self, locus, sequence, imgtdb_version):
        """
        Looks up sequence from

        :param locus: string containing HLA locus.
        :param sequence: string containing sequence data.

        :return: GFEobject.
        """
        # TODO: just initalize seq2hla as DataFrame
        if isinstance(self.seq2hla, DataFrame) and not self.seq2hla.empty:
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
                    seq_features = self.graph.run(get_features(gfe)).to_data_frame()
                    for i in range(0, len(seq_features['term'])):
                        feature = Feature(accession=seq_features['accession'][i],
                                          rank=seq_features['rank'][i],
                                          sequence=seq_features['sequence'][i],
                                          term=lc(seq_features['term'][i]))
                        features.append(feature)
                return [hla, gfe, features]
        else:
            lookup_query = sequence_search(locus, sequence)
            sequence_data = self.graph.run(lookup_query).to_data_frame()
            if not sequence_data.empty:
                features = list()
                gfe = list(set([x for x in sequence_data["GFE"]]))
                hla = list(set([x for x in sequence_data["HLA"]]))
                seq_features = self.graph.run(get_features(gfe[0])).to_data_frame()
                for i in range(0, len(seq_features['term'])):
                    feature = Feature(accession=seq_features['accession'][i],
                                      rank=seq_features['rank'][i],
                                      sequence=seq_features['sequence'][i],
                                      term=lc(seq_features['term'][i]))
                    features.append(feature)
                return [hla[0], gfe[0], features]
            else:
                return

    def gfe_create(self, locus, sequence, imgtdb_version):
        """
        creates GFE from HLA sequence and locus

        :param locus: string containing HLA locus.
        :param sequence: string containing sequence data.

        :return: Dict.
        """
        seq_rec = SeqRecord(Seq(sequence, IUPAC.unambiguous_dna), id="GFE")
        annotation = self.seqann[imgtdb_version].annotate(seq_rec, locus)
        if not hasattr(annotation, 'annotation') or not hasattr(annotation, 'gfe'):
            return {}
        else:
            return {'gfe': annotation.gfe, 'structure': annotation.structure, 'annotation': annotation}

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
                similar_data = self.graph.run(cypher).to_data_frame()
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
                similar_data = self.graph.run(cypher).to_data_frame()
                return self.create_typing(similar_data, gfe, features)
        elif is_kir(gfe):
            return self.find_gfe_kir(gfe, features)
        else:
            return

    def rename_feat(self, f):
        if isutr(f):
            return f
        else:
            return "_".join(f.split("-"))

    def diff_seq(self, ref_allele, annotation, imgtdb_version):
        # TODO: get neo4j seq array back as python array
        #       and then do the comparison in python.
        locus = ref_allele.split("*")[0]
        dbv = "".join(imgtdb_version.split("."))
        feat_order = list(self.seqann[dbv].refdata.struct_order[locus].keys())
        feat_order.sort()
        aligned_seq = list("".join([annotation.aligned[self.rename_feat(self.seqann[dbv].refdata.struct_order[locus][i])] for i in feat_order]))
        seq_l = ",".join(["".join(["\"", s, "\""]) for s in aligned_seq])
        query = "MATCH (a:IMGT_HLA)-[r:HAS_ALIGNMENT]-(n:GEN_ALIGN) " \
                + "WHERE a.locus = \"" + locus + "\" " \
                + "AND a.name = \"" + ref_allele + "\" " \
                + "AND r.imgt_release = \"" + imgtdb_version + "\" " \
                + "WITH [" + seq_l + "] AS RSEQ,n " \
                + "WITH [x in range(0,size(n.seq)-1) | x] AS ind,RSEQ,n " \
                + "UNWIND ind AS number " \
                + "WITH DISTINCT number,RSEQ,n " \
                + "WHERE NOT(n.seq[number] = RSEQ[number]) "\
                + "RETURN number, n.seq[number] AS REF, RSEQ[number] AS SEQ"
        diff_data = self.graph.run(query).to_data_frame()
        diffs = []
        for i in range(0, len(diff_data)):
            index = diff_data['number'][i]
            diff_feat = self.seqann[dbv].refdata.align_coordinates[locus.split("-")[1]][index]
            ref_seq = diff_data['REF'][i]
            in_seq = diff_data['SEQ'][i]
            relative_location = index + self.seqann[dbv].refdata.location[locus]
            # TODO: Change rank!
            diff = Seqdiff(term=diff_feat, location=str(relative_location), ref=ref_seq, inseq=in_seq)
            diffs.append(diff)
        return diffs

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
        similar_data = self.graph.run(cypher).to_data_frame()
        return self.create_typing(similar_data, gfe, features)

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
                    if not feat == "FIVE_PRIME_UTR" and not feat == "THREE_PRIME_UTR":
                        feat_term, feat_rank = feat.split('-')
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
        for feature_rank in self.structures[locus]:
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

    def list_all_db_releases(self):
        cypher = all_db_imgt()
        response_data = self.graph.run(cypher).to_data_frame()
        return list(response_data.HLA_DB)


class SeqAnnException(Exception):

    def __init__(self, inputtype=None, reason=None):
        self.inputtype = inputtype
        self.reason = reason

    def __str__(self):
        """
        Custom error messages for exception
        """
        error_message = "Reason: ({0})\n"\
                        "Type: {1}\n".format(self.reason, self.inputtype)

        return error_message

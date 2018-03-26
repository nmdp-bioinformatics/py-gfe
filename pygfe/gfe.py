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

import os
import glob
import sys
import re
import logging
from pprint import pprint

from pygfe.feature_client.apis.features_api import FeaturesApi
from pygfe.feature_client.api_client import ApiClient
from pygfe.feature_client.rest import ApiException
from pygfe.feature_client.models.feature import Feature
from pygfe.feature_client.models.sequence import Sequence
from pygfe.feature_client.models.feature_request import FeatureRequest

is_kir = lambda x: True if re.search("KIR", x) else False
isutr = lambda f: True if re.search("UTR", f) else False

logging.basicConfig(format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
                    datefmt='%m/%d/%Y %I:%M:%S %p',
                    level=logging.INFO)


class GFE(object):
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
    def __init__(self, url="http://feature.nmdp-bioinformatics.org",
                 loci=['HLA-A', 'HLA-B', 'HLA-C',
                       'HLA-DRB1', 'HLA-DQB1', 'HLA-DRB4',
                       'HLA-DRB5', 'HLA-DPB1', 'HLA-DPA1',
                       'HLA-DQA1', 'HLA-DRB3'],
                 load_features=False, store_features=False,
                 verbose=False,
                 verbosity=1):

        self.loci = loci
        self.verbose = verbose
        self.verbosity = verbosity
        self.store_features = store_features

        client = ApiClient(host=url)
        api_instance = FeaturesApi(api_client=client)
        self.api = api_instance

        structures = {}
        struct_order = {}
        self.all_feats = {loc: {} for loc in loci}
        data_dir = os.path.dirname(__file__)
        struture_files = glob.glob(data_dir + '/data/*.structure')
        for inputfile in struture_files:
            file_path = inputfile.split("/")
            locus = file_path[len(file_path)-1].split(".")[0]
            # TODO: add try
            with open(inputfile, 'r') as f:
                features_order = {}
                features = {}
                n = 0
                for line in f:
                    line = line.rstrip()
                    [feature, rank] = line.split("\t")
                    feature_name = "_".join([feature, rank])
                    if feature == "three_prime_UTR" or feature == "five_prime_UTR":
                        feature_name = feature
                    n += 1
                    features.update({feature_name: n})
                    features_order.update({n: feature_name})
                    if is_kir(locus):
                        structures.update({locus: features})
                        struct_order.update({locus: features_order})
                    else:
                        structures.update({"HLA-" + locus: features})
                        struct_order.update({"HLA-" + locus: features_order})
            f.close()
        self.structures = structures
        self.struct_order = struct_order

        # Load all features from feature service
        if load_features:
            if verbose:
                logging.info("Loading features...")

            # Calling load_features() to load
            # features at each locus
            self.load_features()

    def load_features(self):
        """
        creates GFE from HLA sequence and locus

        :param locus: string containing HLA locus.
        :param sequence: string containing sequence data.

        :return: GFEobject.
        """
        # Loading all loci that
        # are in self.loci variable defined
        # when the pyGFE object is created
        for loc in self.loci:
            if self.verbose:
                logging.info("Loading features for " + loc)

            # Loading all features for loc from feature service
            self.all_feats.update({loc: self.locus_features(loc)})

            if self.verbose:
                logging.info("Finished loading features for " + loc)

        if self.verbose:
            mem = "{:4.4f}".format(sys.getsizeof(self.all_feats) / 1000000)
            logging.info("Finished loading all features * all_feats = " + mem + " MB *")

    def locus_features(self, locus):
        """
        creates GFE from HLA sequence and locus

        :param locus: string containing HLA locus.
        :param sequence: string containing sequence data.

        :return: GFEobject.
        """
        features = self.api.list_features(locus=locus)
        feat_dict = {":".join([a.locus, str(a.rank), a.term, a.sequence]): a.accession for a in features}
        return feat_dict

    def get_gfe(self, annotation, locus):
        """
        creates GFE from HLA sequence and locus

        :param locus: string containing HLA locus.
        :param sequence: string containing sequence data.

        :return: GFEobject.
        """
        features = []
        accessions = {}
        for feat in annotation.annotation:
            seq = str(annotation.annotation[feat].seq)

            # TODO: Drop this if statement
            if isutr(feat):
                feat_str = ":".join([locus, str(1), feat, seq])

                # If the feature has been loaded or stored
                # then use that instead of making a feature request
                if self.verbose and self.verbosity > 2:
                    logging.info("Getting accession " + feat_str)

                if feat_str in self.all_feats[locus]:

                    if self.verbose and self.verbosity > 2:
                        logging.info("Feature found " + feat_str)

                    accession = self.all_feats[locus][feat_str]
                    feature = Feature(term=feat,
                                      rank=1,
                                      locus=locus,
                                      sequence=seq,
                                      accession=accession)
                    accessions.update({feat: accession})
                    features.append(feature)
                else:

                    if self.verbose and self.verbosity > 2:
                        logging.info("Making FeatureRequest " + feat_str)

                    # Create FeatureRequest object
                    request = FeatureRequest(locus=locus,
                                             term=feat,
                                             rank=1,
                                             sequence=seq)

                    # Attempt to make feature request
                    try:
                        feature = self.api.create_feature(body=request)
                        accessions.update({feat: feature.accession})
                        features.append(feature)
                    except ApiException as e:
                        logging.warn("Exception when calling DefaultApi->create_feature" + e)
                        blank_feat = Feature(term=feat, rank=1, locus=locus,
                                             sequence=seq)
                        accessions.update({feat: 0})
                        features.append(blank_feat)

                    # Store new features for quick retrieval if flag passed
                    if self.store_features:

                        # Adding new feature to all_feats
                        self.all_feats[locus].update({feat_str: feature.accession})

                        # Calculating memory size of all_feats
                        if self.verbose and self.verbosity > 1:
                            logging.info("Storing new feature " + feat_str)
                            mem = "{:4.4f}".format(sys.getsizeof(self.all_feats) / 1000000)
                            logging.info("Updated * all_feats " + mem + " MB *")

            else:
                term, rank = feat.split("_")
                feat_str = ":".join([locus, str(rank), term, seq])

                # If the feature has been loaded or stored
                # then use that instead of making a feature request
                if feat_str in self.all_feats[locus]:

                    if self.verbose and self.verbosity > 2:
                        logging.info("Feature found " + feat_str)

                    accession = self.all_feats[locus][feat_str]
                    feature = Feature(term=term,
                                      rank=rank,
                                      locus=locus,
                                      sequence=seq,
                                      accession=accession)
                    accessions.update({feat: accession})
                    features.append(feature)
                else:

                    if self.verbose and self.verbosity > 2:
                        logging.info("Making FeatureRequest " + feat_str)

                    # Create FeatureRequest object
                    request = FeatureRequest(locus=locus,
                                             term=term,
                                             rank=rank,
                                             sequence=seq)

                    # Attempt to make feature request
                    try:
                        feature = self.api.create_feature(body=request)
                        accessions.update({feat: feature.accession})
                        features.append(feature)
                    except ApiException as e:
                        print(request)
                        logging.warn("Exception when calling DefaultApi->create_feature %e" + e)
                        blank_feat = Feature(term=term, rank=rank, locus=locus,
                                             sequence=seq)
                        accessions.update({feat: 0})
                        features.append(blank_feat)

                    # Store new features for quick retrieval if flag passed
                    if self.store_features:

                        # Adding new feature to all_feats
                        self.all_feats[locus].update({feat_str: feature.accession})

                        # Calculating memory size of all_feats
                        if self.verbose and self.verbosity > 1:
                            logging.info("Storing new feature " + feat_str)
                            mem = "{:4.4f}".format(sys.getsizeof(self.all_feats) / 1000000)
                            logging.info("Updated * all_feats " + mem + " MB *")

        # Creating GFE
        gfe = self._make_gfe(accessions, locus)

        if self.verbose:
            logging.info("GFE = " + gfe)

        return features, gfe

    def get_sequence(self, gfe):
        """
        creates GFE from HLA sequence and locus

        :param locus: string containing HLA locus.
        :param sequence: string containing sequence data.

        :return: GFEobject.
        """
        feats = []
        seqs = []
        loc, accessions = gfe.split("w")
        features = self._breakup_gfe(gfe)
        for f in features:
            if int(features[f]) != 0:
                if isutr(f):
                    feat = self._seq(loc, f, 1,  features[f])
                    seqs.append(feat.sequence)
                    feats.append(feat)
                else:

                    feat = self._seq(loc, f.split("_")[0], f.split("_")[1],
                                     features[f])
                    seqs.append(feat.sequence)
                    feats.append(feat)
        seq = "".join(seqs)
        sequence_o = Sequence(sequence=seq, structure=feats)
        return sequence_o

    def _seq(self, locus, term, rank, accession):
        """
        creates GFE from HLA sequence and locus

        :param locus: string containing HLA locus.
        :param sequence: string containing sequence data.

        :return: GFEobject.
        """
        try:
            feature = self.api.get_feature_by_path(locus,
                                                   term,
                                                   rank,
                                                   accession)
            return feature
        except ApiException as e:
            print("Exception when calling DefaultApi->get_feature_by_path: %s\n" % e)
            return ''

    def _breakup_gfe(self, gfe):
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
            features.update({feature_rank: accession})
            i += 1

        return(features)

    def _make_gfe(self, features, locus):
        """
        creates GFE from HLA sequence and locus

        :param locus: string containing HLA locus.
        :param sequence: string containing sequence data.

        :return: GFEobject.
        """
        gfe_list = []
        for feat in sorted(self.structures[locus],
                           key=lambda k: self.structures[locus][k]):

            acc = str(0)
            if feat in features:
                acc = str(features[feat])
            gfe_list.append(acc)

        gfea = '-'.join(gfe_list)
        return locus + "w" + gfea


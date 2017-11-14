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
import re


from pygfe.feature_client.apis.features_api import FeaturesApi
from pygfe.feature_client.api_client import ApiClient
from pygfe.feature_client.rest import ApiException
from pprint import pprint
from pygfe.feature_client.models.feature import Feature
from pygfe.feature_client.models.feature_request import FeatureRequest

is_kir = lambda x: True if re.search("KIR", x) else False
isutr = lambda f: True if re.search("UTR", f) else False


class pyGFE(object):
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
    def __init__(self, url="http://feature.nmdp-bioinformatics.org", verbose=False):
        self.verbose = verbose
        client = ApiClient(host=url)
        api_instance = FeaturesApi(api_client=client)
        self.api = api_instance

        structures = {}
        struct_order = {}
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

    def get_gfe(self, annotation, locus):

        features = []
        accessions = {}
        for feat in annotation.annotation:
            if isutr(feat):
                seq = str(annotation.annotation[feat].seq)
                request = FeatureRequest(locus=locus,
                                         term=feat,
                                         rank=1,
                                         sequence=seq)

                feature = self.api.create_feature(body=request)
                accessions.update({feat: feature.accession})
                features.append(feature)
            else:
                term, rank = feat.split("_")
                seq = str(annotation.annotation[feat].seq)
                request = FeatureRequest(locus=locus,
                                         term=term,
                                         rank=rank,
                                         sequence=seq)
                try:
                    feature = self.api.create_feature(body=request)
                    accessions.update({feat: feature.accession})
                    features.append(feature)
                except ApiException as e:
                    print("Exception when calling DefaultApi->create_feature: %s\n" % e)
                    blank_feat = Feature(term=term, rank=rank, locus=locus,
                                         sequence=seq)
                    accessions.update({feat: 0})
                    features.append(blank_feat)

        gfe = self._make_gfe(accessions, locus)
        return features, gfe

    def get_sequence(self, gfe):

        loc, accessions = gfe.split("w")
        features = self._breakup_gfe(gfe)
        seqs = []
        for f in features:
            if isutr(f):
                seqs.append(self._seq(loc, f, 1,  features[f]))
            else:
                seqs.append(self._seq(loc, f.split("_")[0], f.split("_")[1],
                            features[f]))
        seq = "".join(seqs)
        return seq

    def _seq(self, locus, term, rank, accession):
        try:
            feature = self.api.get_feature_by_path(locus,
                                                   term,
                                                   rank,
                                                   accession)
            return feature.sequence
        except ApiException as e:
            print("Exception when calling DefaultApi->get_feature_by_path: %s\n" % e)
            return ''

    def _breakup_gfe(self, gfe):
        [locus, feature_accessions] = gfe.split("w")
        accessions = feature_accessions.split("-")

        if locus == "HLA-DQB1":
            if len(accessions) < len(self.structures[locus]):
                i = 0
                features = {}
                old_dq = ["intro_1","exon_2","intron_2","exon_3","intro_3"]
                for feature_rank in old_dq:
                    accession = accessions[i]
                    features.update({feature_rank: accession})
                    i += 1

                return(features)
            else:
                i = 0
                features = {}
                for feature_rank in self.structures[locus]:
                    accession = accessions[i]
                    features.update({feature_rank: accession})
                    i += 1

                return(features)
        else:
            i = 0
            features = {}
            for feature_rank in self.structures[locus]:
                accession = accessions[i]
                features.update({feature_rank: accession})
                i += 1

            return(features)

    def _make_gfe(self, features, locus):

        gfe_list = []
        for feat in sorted(self.structures[locus],
                           key=lambda k: self.structures[locus][k]):

            acc = str(0)
            if feat in features:
                acc = str(features[feat])
            gfe_list.append(acc)

        gfea = '-'.join(gfe_list)
        return locus + "w" + gfea


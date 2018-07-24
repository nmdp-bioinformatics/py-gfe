'''
Created on Feb 8, 2017

@author: mhalagan
'''

from typing import List, Dict


def hla_seqdiff(locus: str, db: str, typing1: str, typing2: str):
    query = "MATCH (a:IMGT_HLA)-[r:HAS_ALIGNMENT]-(n:GEN_ALIGN)," \
            + "(a2:IMGT_HLA)-[r2:HAS_ALIGNMENT]-(n2:GEN_ALIGN)" \
            + " WHERE a.locus = \"" + locus + "\"" \
            + " AND a2.locus = \"" + locus + "\"" \
            + " AND a2.name = \"" + typing2 + "\"" \
            + " AND a.name = \"" + typing1 + "\"" \
            + " AND r2.imgt_release = \"" + db + "\"" \
            + " AND r.imgt_release = \"" + db + "\"" \
            + " WITH [x in range(0,size(n.seq)-1) | x] AS ind,n,n2" \
            + " UNWIND ind AS number" \
            + " WITH DISTINCT number,n,n2" \
            + " WHERE NOT(n.seq[number] = n2.seq[number])" \
            + " RETURN number AS POS, n.seq[number] AS VAR1, n2.seq[number] AS VAR2"
    print(query)
    return query

# def hla_seqdiff(typing1: str, typing2: str) -> str:


# def gfe_seqdiff(typing1: str, typing2: str) -> str:


# def get_variants(locus: str, positions: List) -> List:


# def get_sequence(typing: str, position: int, seq_type: str) -> str:


# def seq_seqdiff(sequence: str) -> str:
# MATCH (a:GFE)-[r:HAS_ALIGNMENT]-(n:NUC_ALIGN)
# WHERE a.locus = "HLA-A"
# AND r.imgt_release = "3.31.0"
# WITH ["-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","G","C","T","C","C","C","A","C","T","C","C","A","T","G","A","G","G","T","A","T","T","T","C","T","T","C","A","C","A","T","C","C","G","T","G","T","C","C","C","-","-","-","-","-","G","G","C","C","C","G","G","C","C","G","C","G","G","G","G","A","G","C","C","C","C","G","C","T","T","C","A","T","C","G","C","C","G","T","G","G","G","C","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","T","A","C","G","T","G","G","A","C","G","A","C","A","C","G","C","A","G","T","T","C","G","T","G","C","G","G","T","T","C","G","A","C","A","G","C","G","A","C","G","C","C","G","C","G","A","G","C","C","A","G","A","-","-","A","G","A","T","G","G","A","G","C","C","G","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","C","G","G","G","C","G","C","C","G","T","G","G","A","T","A","G","A","G","C","A","G","G","A","G","G","G","G","C","C","G","G","A","G","T","A","T","T","G","G","G","A","C","C","A","G","G","A","G","A","C","A","C","G","G","A","A","T","A","T","G","A","A","G","G","C","C","C","-","A","C","T","C","A","C","A","G","A","C","T","G","A","C","C","G","A","G","C","G","A","A","C","C","T","G","G","G","G","A","C","C","C","T","G","C","G","C","G","G","C","T","A","C","T","A","C","A","A","C","C","A","G","A","G","C","G","A","G","G","A","C","G","G","T","T","C","T","C","A","C","A","C","C","A","T","C","C","A","G","A","T","A","A","T","G","T","A","T","G","G","C","T","G","C","G","A","C","G","T","G","G","G","G","C","C","G","G","A","C","G","G","G","C","G","C","T","T","C","C","T","C","C","G","C","G","G","G","T","A","C","C","G","G","C","A","G","G","A","C","G","C","C","T","A","C","G","A","C","G","G","C","A","A","G","G","A","T","T","A","C","A","T","C","G","C","C","C","T","G","A","A","C","G","A","G","G","A","C","C","T","G","C","G","C","T","C","T","T","G","G","A","C","C","G","C","G","G","C","G","G","A","C","A","T","G","G","C","A","G","C","T","C","A","G","A","T","C","A","C","C","A","A","G","C","G","C","A","A","G","T","G","G","G","A","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","G","-","-","G","C","G","G","T","C","C","A","T","G","C","G","G","C","-","G","G","A","G","C","A","G","C","G","G","A","G","A","G","T","C","T","A","C","C","T","G","G","A","G","G","G","C","C","G","G","T","G","C","G","T","G","G","A","C","G","G","G","C","T","C","C","G","C","A","G","A","T","A","C","T","T","G","G","A","G","A","A","C","G","G","G","A","A","G","G","A","G","A","C","G","C","T","G","C","A","G","C","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","G","C","A","C","G","G","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-"] AS RSEQ,n,a
# WITH [x in range(0,size(n.seq)-1) | x] AS ind,RSEQ,n,a
# UNWIND ind AS number
# WITH DISTINCT number,RSEQ,n,a
# WHERE NOT(n.seq[number] = RSEQ[number])
# RETURN a.name,number, n.seq[number], RSEQ[number]
# LIMIT 10

# MATCH (a:GFE)-[r:HAS_ALIGNMENT]-(n:NUC_ALIGN)
# WHERE a.locus = "HLA-A"
# AND r.imgt_release = "3.31.0"
# WITH ["-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","G","C","T","C","C","C","A","C","T","C","C","A","T","G","A","G","G","T","A","T","T","T","C","T","T","C","A","C","A","T","C","C","G","T","G","T","C","C","C","-","-","-","-","-","G","G","C","C","C","G","G","C","C","G","C","G","G","G","G","A","G","C","C","C","C","G","C","T","T","C","A","T","C","G","C","C","G","T","G","G","G","C","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","T","A","C","G","T","G","G","A","C","G","A","C","A","C","G","C","A","G","T","T","C","G","T","G","C","G","G","T","T","C","G","A","C","A","G","C","G","A","C","G","C","C","G","C","G","A","G","C","C","A","G","A","-","-","A","G","A","T","G","G","A","G","C","C","G","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","C","G","G","G","C","G","C","C","G","T","G","G","A","T","A","G","A","G","C","A","G","G","A","G","G","G","G","C","C","G","G","A","G","T","A","T","T","G","G","G","A","C","C","A","G","G","A","G","A","C","A","C","G","G","A","A","T","A","T","G","A","A","G","G","C","C","C","-","A","C","T","C","A","C","A","G","A","C","T","G","A","C","C","G","A","G","C","G","A","A","C","C","T","G","G","G","G","A","C","C","C","T","G","C","G","C","G","G","C","T","A","C","T","A","C","A","A","C","C","A","G","A","G","C","G","A","G","G","A","C","G","G","T","T","C","T","C","A","C","A","C","C","A","T","C","C","A","G","A","T","A","A","T","G","T","A","T","G","G","C","T","G","C","G","A","C","G","T","G","G","G","G","C","C","G","G","A","C","G","G","G","C","G","C","T","T","C","C","T","C","C","G","C","G","G","G","T","A","C","C","G","G","C","A","G","G","A","C","G","C","C","T","A","C","G","A","C","G","G","C","A","A","G","G","A","T","T","A","C","A","T","C","G","C","C","C","T","G","A","A","C","G","A","G","G","A","C","C","T","G","C","G","C","T","C","T","T","G","G","A","C","C","G","C","G","G","C","G","G","A","C","A","T","G","G","C","A","G","C","T","C","A","G","A","T","C","A","C","C","A","A","G","C","G","C","A","A","G","T","G","G","G","A","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","G","-","-","G","C","G","G","T","C","C","A","T","G","C","G","G","C","-","G","G","A","G","C","A","G","C","G","G","A","G","A","G","T","C","T","A","C","C","T","G","G","A","G","G","G","C","C","G","G","T","G","C","G","T","G","G","A","C","G","G","G","C","T","C","C","G","C","A","G","A","T","A","C","T","T","G","G","A","G","A","A","C","G","G","G","A","A","G","G","A","G","A","C","G","C","T","G","C","A","G","C","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","G","C","A","C","G","G","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-","-"] AS RSEQ,n,a
# WITH [x in range(0,size(n.seq)-1) | x] AS ind,RSEQ,n,a
# UNWIND ind AS number
# WITH DISTINCT number,RSEQ,n,a
# WHERE NOT(n.seq[number] = RSEQ[number])
# RETURN a.name,collect(distinct {number: number, seq: n.seq[number], rseq: RSEQ[number]})
# LIMIT 10

def all_gfe2feats() -> str:
    query = "MATCH(gfe:GFE)-[r:HAS_FEATURE]-(feat:FEATURE)" \
        + "RETURN r.imgt_release AS DB,gfe.name AS GFE," \
        + "collect(distinct {term: feat.name, rank: feat.rank, accession: r.accession, sequence: feat.sequence}) AS FEATS"
    return query


def all_seq2hla() -> str:
    query = "MATCH (hla:IMGT_HLA)-[f1:HAS_GFE]-(gfe:GFE)-[f2:HAS_SEQUENCE]-(seq:SEQUENCE) " \
        + "WHERE f1.imgt_release = f2.imgt_release " \
        + "RETURN DISTINCT f1.imgt_release AS DB,hla.locus AS LOC,hla.name AS HLA, gfe.name AS GFE, seq.sequence AS SEQ"
    return query


def all_gfe2hla() -> str:
    query = "MATCH(hla:IMGT_HLA)-[r1:HAS_GFE]-(gfe:GFE)" \
        + "RETURN DISTINCT r1.imgt_release AS DB,hla.locus AS LOC," \
        + "hla.name AS HLA, gfe.name AS GFE"
    return query


def all_feats() -> str:
    query = "MATCH(hla:IMGT_HLA)-[r1:HAS_FEATURE]-(feat:FEATURE)" \
        + "RETURN DISTINCT hla.locus AS LOC,r1.imgt_release AS DB," \
        + "feat.rank AS RANK,feat.name AS TERM,r1.accession AS ACCESSION,feat.sequence AS SEQ"
    return query


def search_hla_features(locus: str, gfe_feats: List) -> str:

    # TODO: Fix error observed with DQB1
    match = "MATCH(hla:IMGT_HLA)-[:HAS_GFE]-(gfe:GFE)-[f1:HAS_FEATURE]-(feat1:FEATURE)"
    if(len(gfe_feats)) > 1:
        for i in range(1, len(gfe_feats)):
            j = i + 1
            match = match + ",(hla:IMGT_HLA)-[:HAS_GFE]-(gfe:GFE)-[f" + str(j) + ":HAS_FEATURE]-(feat" + str(j) + ":FEATURE)"

    i = 1
    feat_q = "WHERE hla.locus = \"" + locus + "\""

    for feat in gfe_feats:
        [term, rank] = feat.split("-")
        feat_q = feat_q + " AND feat" + str(i) + ".name = \"" + term.upper() + "\""
        feat_q = feat_q + " AND feat" + str(i) + ".rank = \"" + rank + "\""
        acc_q = " AND( f" + str(i) + ".accession = \"" + gfe_feats[feat][0] + "\""
        if(len(gfe_feats[feat]) == 1):
            acc_q = acc_q + ")"
        else:
            for j in range(1, len(gfe_feats[feat])):
                acc_q = acc_q + " OR f" + str(i) + ".accession = \"" + gfe_feats[feat][j] + "\""
            acc_q = acc_q + ")"
        feat_q = feat_q + acc_q
        i += 1

    return_q = " RETURN DISTINCT hla.name AS HLA, gfe.name AS GFE"
    return match + feat_q + return_q


def fullseqid(seq: str) -> str:
    q1 = "MATCH(seq:SEQUENCE)"
    q2 = " WHERE seq.sequence = \"" + seq + "\""
    q3 = " RETURN seq.sequenceId AS ID"
    return q1 + q2 + q3


def seqid(seq: str, seqtype: str, rank: int) -> str:
    q1 = "MATCH(feat:FEATURE)"
    q2 = " WHERE feat.sequence = \"" + seq + "\""
    q3 = " AND feat.rank = \"" + str(rank) + "\""
    q4 = " AND feat.name = \"" + seqtype.upper() + "\""
    q5 = "RETURN feat.sequenceId AS ID"
    return q1 + q2 + q3 + q4 + q5


def hla_alleleid(hla: str) -> str:
    q1 = "MATCH(hla:IMGT_HLA) "
    q2 = "WHERE hla.name = \"" + hla + "\" "
    q3 = "RETURN hla.alleleId AS ID"
    return q1 + q2 + q3


def gfe_alleleid(gfe: str) -> str:
    q1 = "MATCH(gfe:GFE) "
    q2 = "WHERE gfe.name = \"" + gfe + "\" "
    q3 = "RETURN gfe.alleleId AS ID"
    return q1 + q2 + q3


def get_features(gfe: str) -> str:

    q1 = "MATCH(gfe:GFE)-[f1:HAS_FEATURE]-(f:FEATURE)"
    q2 = "WHERE gfe.name = \"" + gfe + "\""
    q4 = "RETURN f.name AS term,f.rank AS rank,f1.accession AS accession,f.sequence AS sequence"
    return q1 + q2 + q4


def sequence_search(locus: str, sequence: str) -> str:
    seq_query = "MATCH (hla:IMGT_HLA)-[f1:HAS_GFE]-(gfe:GFE)-[f2:HAS_SEQUENCE]-(seq:SEQUENCE)"
    seq_query1 = " WHERE seq.sequence = \"" + sequence + "\""
    seq_query2 = " AND hla.locus = \"" + locus + "\""
    seq_query5 = " RETURN hla.name AS HLA, gfe.name AS GFE"
    query = seq_query + seq_query1 + seq_query2 + seq_query5
    return(query)


def gfe_search(gfe: str) -> str:
    seq_query = "MATCH (hla:IMGT_HLA)-[f:HAS_GFE]-(gfe:GFE)"
    seq_query1 = " WHERE gfe.name = \"" + gfe + "\""
    seq_query2 = " AND NOT(f.status = \"persisted\")"
    seq_query3 = " RETURN hla.name AS HLA"
    query = seq_query + seq_query1 + seq_query2 + seq_query3
    return(query)


def similar_gfe_classII(gfe: str, exon2: int, db: str) -> str:
    [locus, feature_accessions] = gfe.split("w")
    typing_query = "MATCH (hla:IMGT_HLA)-[r:HAS_GFE]-(gfe1:GFE)-[f1:HAS_FEATURE]-(feat1:FEATURE)" \
            + " WHERE gfe1.locus = \"" + locus + "\"" \
            + " AND f1.accession = \"" + exon2 + "\"" \
            + " AND feat1.name = \"EXON\"" \
            + " AND r.imgt_release = \"" + str(db) + "\"" \
            + " AND feat1.rank = \"2\"" \
            + " WITH COLLECT( distinct {hla: hla.name, gfe: gfe1.name}) AS TYPING" \
            + " UNWIND TYPING AS HLA_GFE" \
            + "  RETURN HLA_GFE.hla AS HLA,HLA_GFE.gfe AS GFE"
    return(typing_query)


def similar_gfe_classI(gfe: str, exon2: int, exon3: int, db: str) -> str:

    [locus, feature_accessions] = gfe.split("w")
    typing_query = "MATCH (hla:IMGT_HLA)-[r:HAS_GFE]-(gfe1:GFE)-[f1:HAS_FEATURE]-(feat1:FEATURE)" \
            + " MATCH (hla:IMGT_HLA)-[r2:HAS_GFE]-(gfe1:GFE)-[f2:HAS_FEATURE]-(feat2:FEATURE)" \
            + " WHERE gfe1.locus = \"" + locus + "\"" \
            + " AND f1.accession = \"" + exon2 + "\"" \
            + " AND f2.accession = \"" + exon3 + "\"" \
            + " AND r.imgt_release = \"" + str(db) + "\"" \
            + " AND r2.imgt_release = \"" + str(db) + "\"" \
            + " AND feat1.name = \"EXON\"" \
            + " AND feat2.name = \"EXON\"" \
            + " AND feat1.rank = \"2\"" \
            + " AND feat2.rank = \"3\"" \
            + " WITH COLLECT( distinct {hla: hla.name, gfe: gfe1.name}) AS TYPING" \
            + " UNWIND TYPING AS HLA_GFE" \
            + "  RETURN HLA_GFE.hla AS HLA,HLA_GFE.gfe AS GFE"
    return(typing_query)


def similar_kir(locus: str) -> str:
    query = "MATCH(gfe:GFE)" \
        + "WHERE gfe.locus = \"" + locus + "\"" \
        + "RETURN DISTINCT gfe.name AS GFE"
    return query


def hla_Ggroups(hla: str) -> str:
    matchq = "MATCH (G:G_GROUP)-[:IN_GROUP]-(hla:IMGT_HLA) WHERE hla.name = \"" + hla + "\""
    returnq = " RETURN DISTINCT G.name as G_GROUP"
    cyper_query = matchq + returnq
    return(cyper_query)


def gfe_Ggroups(gfe: str) -> str:
    matchq = "MATCH (G:G_GROUP)-[:IN_GROUP]-(hla:IMGT_HLA)-[:HAS_GFE]-(gfe1:GFE) WHERE gfe1.name = \"" + gfe + "\""
    returnq = " RETURN DISTINCT G.name as G_GROUP"
    cyper_query = matchq + returnq
    return(cyper_query)


def gfe_hla(gfe: str) -> str:
    matchq = "MATCH (hla:IMGT_HLA)-[h:HAS_GFE]-(gfe1:GFE) WHERE gfe1.name = \"" + gfe + "\""
    withq = " WITH collect(DISTINCT hla.name) as HLA"
    returnq = " RETURN HLA"
    cyper_query = matchq + withq + returnq
    return(cyper_query)


def groups_classI(locus: str, group: str, exon2: int, exon3: int) -> str:

    match = "MATCH (group:" + group + ")-[:IN_GROUP]-(hla:IMGT_HLA)-[:HAS_GFE]-(gfe1:GFE)-[f1:HAS_FEATURE]->(feat1:FEATURE)"
    match2 = ",(gfe1:GFE)-[f2:HAS_FEATURE]->(feat2:FEATURE) "
    where = " WHERE feat1.rank = '2' AND feat2.name = \"EXON\" AND feat1.name = \"EXON\" AND hla.locus = \"" + locus + "\""
    where2 = " AND feat2.rank = '3' AND f2.accession = \"" + exon3 + "\" AND f1.accession = \"" + exon2 + "\""
    withst = " WITH collect(DISTINCT hla.name) as HLA,collect(DISTINCT gfe1.name) as GFE,collect(DISTINCT group.name) as ARS"
    returnc = " RETURN HLA,GFE,ARS"
    cypher = match + match2 + where + where2 + withst + returnc
    return(cypher)


def groups_classII(locus: str, group: str, exon3: int) -> str:

    match = "MATCH (group:" + group + ")-[:IN_GROUP]-(hla:IMGT_HLA)-[:HAS_GFE]-(gfe1:GFE)-[f1:HAS_FEATURE]->(feat1:FEATURE)"
    match2 = ",(gfe1:GFE)-[f2:HAS_FEATURE]->(feat2:FEATURE) "
    where = " WHERE feat1.rank = '3' AND feat1.name = \"EXON\" AND hla.locus = \"" + locus + "\""
    where2 = " AND f1.accession = \"" + exon3 + "\""
    withst = " WITH collect(DISTINCT hla.name) as HLA,collect(DISTINCT gfe1.name) as GFE"
    returnc = " RETURN HLA,GFE"
    cypher = match + match2 + where + where2 + withst + returnc
    return(cypher)


def search_feature(term: str, rank: int, sequence: str) -> str:
    q1 = "MATCH(gfe:GFE)-[h:HAS_FEATURE]-(feat:FEATURE) WHERE feat.name = \"" + term.upper() + "\" "
    q2 = "AND feat.rank = \"" + str(rank) + "\" "
    q3 = "AND feat.sequence = \"" + sequence + "\" "
    q4 = "RETURN feat.name"
    return q1 + q2 + q3 + q4


def hla_gfe(hla: str) -> str:
    seq_query = "MATCH (hla:IMGT_HLA)-[h:HAS_GFE]-(gfe:GFE)"
    seq_query1 = " WHERE hla.name = \"" + hla + "\""
    unq = " AND NOT(h.status = \"persisted\") "
    seq_query2 = " RETURN DISTINCT gfe.name AS GFE"
    query = seq_query + seq_query1 + unq + seq_query2
    return(query)


def hla_ars(group: str, hla: str) -> str:
    matchq = "MATCH (group:" + group + ")-[:IN_GROUP]-(hla:IMGT_HLA)-[h:HAS_GFE]-(g:GFE) WHERE hla.name = \"" + hla + "\""
    unq = " AND NOT(h.status = \"persisted\")"
    returnq = " RETURN group.name as ARS, hla.name as HLA,g.name AS GFE"
    cyper_query = matchq + unq + returnq
    return(cyper_query)


def gfe_ars(group: str, gfe: str) -> str:
    matchq = "MATCH (group:" + group + ")-[:IN_GROUP]-(hla:IMGT_HLA)-[:HAS_GFE]-(g:GFE) WHERE g.name = \"" + gfe + "\""
    returnq = " RETURN group.name as ARS, hla.name as HLA,g.name AS GFE"
    cyper_query = matchq + returnq
    return(cyper_query)


def get_sequence(seqtype: str, allele: str) -> str:
    seq_query = "MATCH (allele:" + seqtype + ")-[h:HAS_FEATURE]-(seq:SEQUENCE)"
    seq_query1 = " WHERE allele.name = \"" + allele + "\""
    unq = " AND NOT(h.status = \"persisted\")"
    seq_query2 = " RETURN DISTINCT seq.sequence AS SEQ"
    query = seq_query + seq_query1 + unq + seq_query2
    return(query)


def persisted_query() -> str:
    query = "MATCH(hla:IMGT_HLA)-[g:HAS_GFE]-(gfe:GFE)-[f:HAS_FEATURE]-(feat:FEATURE) "
    q2 = "WHERE g.status = \"persisted\" "
    q3 = "AND g.status = \"persisted\" "
    q4 = "RETURN hla.name AS HLA,gfe.name AS GFE,feat.name AS TERM,feat.rank AS RANK,f.accession AS ACCESSION,feat.sequence AS SEQUENCE"
    return query + q2 + q3 + q4


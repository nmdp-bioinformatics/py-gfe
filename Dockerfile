FROM python:3.7
MAINTAINER NMDP Bioinformatics

RUN apt-get update -q \
    && apt-get install clustalo -y \
	  && apt-get install ncbi-blast+ -y \
    && apt-get autoremove \
    && apt-get clean

RUN pip install --no-cache-dir seq-ann==1.1.0
RUN pip install --no-cache-dir py-gfe==1.1.3

RUN touch blank.fasta \
    && seq2gfe -f blank.fasta -l HLA-A \
	  && seq2gfe -f blank.fasta -l KIR3DL2 -k

ENV BIOSQLHOST=imgt_biosqldb 
ENV BIOSQLPORT=3306

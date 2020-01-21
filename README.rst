===============================
py-gfe
===============================


.. image:: https://img.shields.io/pypi/v/pygfe.svg
        :target: https://pypi.python.org/pypi/pygfe

.. image:: https://img.shields.io/travis/mhalagan-nmdp/pygfe.svg
        :target: https://travis-ci.org/mhalagan-nmdp/pygfe

.. image:: https://readthedocs.org/projects/pygfe/badge/?version=latest
        :target: https://pygfe.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status

.. image:: https://pyup.io/repos/github/mhalagan-nmdp/pygfe/shield.svg
     :target: https://pyup.io/repos/github/mhalagan-nmdp/pygfe/
     :alt: Updates


Python Boilerplate contains all the boilerplate you need to create a Python package.


* Free software: LGPL 3.0
* Documentation: https://pygfe.readthedocs.io.

Docker
--------

.. code-block:: shell

  docker pull nmdpbioinformatics/py-gfe

.. code-block:: 

	docker run -it --rm -v $PWD:/opt nmdpbioinformatics/py-gfe seq2gfe \
		-f /opt/your_fastafile.fasta -l HLA-A


Example
--------

.. code-block:: python3

    >>> from Bio import SeqIO
    >>> from BioSQL import BioSeqDatabase
    >>> from seqann.sequence_annotation import BioSeqAnn
    >>> import pygfe
    >>> seq_file = 'test_dq.fasta'
    >>> gfe = pygfe.pyGFE()
    >>> server = BioSeqDatabase.open_database(driver="pymysql", user="root",
    ...                                       passwd="", host="localhost",
    ...                                      db="bioseqdb")
    >>> seqann = BioSeqAnn(server=server)
    >>> seq_rec = list(SeqIO.parse(seq_file, 'fasta'))[0]
    >>> annotation = seqann.annotate(seq_rec, "HLA-DQB1")
    >>> gfe = gfe.get_gfe(annotation, "HLA-DQB1")
    >>> print(gfe)
    HLA-DQB1w0-4-0-141-0-12-0-4-0-0-0-0-0

Credits
---------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage


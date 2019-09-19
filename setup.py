#!/usr/bin/env python
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

from setuptools import setup

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = [
    'atomicwrites==1.1.5',
    'attrs==18.1.0',
    'biopython==1.70',
    'bson==0.5.5',
    'certifi==2017.11.5',
    'chardet==3.0.4',
    'click==6.7',
    'colorama==0.3.9',
    'idna==2.6',
    'more-itertools==4.2.0',
    'neo4j-driver==1.6.1',
    'neotime==1.0.0',
    'numpy==1.17.2',
    'pandas==0.25.1',
    'pkginfo==1.4.1',
    'pluggy==0.6.0',
    'prompt-toolkit==1.0.15',
    'py==1.5.4',
    'py2neo==4.0.0',
    'Pygments==2.2.0',
    'PyMySQL==0.7.11',
    'pytest==3.6.3',
    'python-dateutil==2.6.1',
    'pytz==2017.3',
    'requests==2.18.4',
    'requests-toolbelt==0.8.0',
    'seqann',
    'six==1.11.0',
    'tqdm==4.19.4',
    'twine==1.9.1',
    'urllib3==1.22',
    'wcwidth==0.1.7'
]

test_requirements = [
    # TODO: put package test requirements here
    'bson',
    'certifi',
    'pytz',
    'six',
    'urllib3'
]

setup(
    name='pygfe',
    version='0.0.25',
    description="Python package for converting sequence annotations to gene feature enumerations (GFE).",
    long_description=readme + '\n\n' + history,
    author="Mike Halagan",
    author_email='mhalagan@nmdp.org',
    url='https://github.com/mhalagan-nmdp/pygfe',
    packages=[
        'pygfe',
        'pygfe.models',
        'pygfe.feature_client',
        'pygfe.feature_client.apis',
        'pygfe.feature_client.models'
    ],
    package_dir={'pygfe':
                 'pygfe'},
    install_requires=requirements,
    license="LGPL 3.0",
    zip_safe=False,
    keywords='pygfe',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)',
        'Natural Language :: English',
        "Programming Language :: Python :: 2",
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
    ],
    scripts=['scripts/seq2gfe'],
    test_suite='tests',
    tests_require=test_requirements
)

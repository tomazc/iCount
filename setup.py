#!/usr/bin/env python
"""Computational pipeline for analysis of iCLIP data.

See:
https://github.com/tomazc/iCount
"""
from os import path

# Use codecs' open for a consistent encoding
from codecs import open
from setuptools import setup, find_packages

# base_dir = path.abspath(path.dirname(__file__))
base_dir = path.dirname(path.realpath(__file__))

# Get the long description from the README file
# with open(path.join(base_dir, 'README.md'), encoding='utf-8') as f:
#    long_description = f.read()

try:
    import pypandoc
    LONG_DESC = pypandoc.convert_file('README.md', 'rst')
    LONG_DESC = LONG_DESC.replace('\r', '')
except(IOError, ImportError):
    with open(path.join(base_dir, 'README.md'), encoding='utf-8') as f:
        LONG_DESC = f.read()
    # Skip header with badges
    LONG_DESC = LONG_DESC.split('\n')[8:]
    LONG_DESC = '\n'.join(LONG_DESC)

# Get package metadata from 'iCount/__about__.py' file
about = {}
with open(path.join(base_dir, 'iCount', '__about__.py'), encoding='utf-8') as f:
    exec(f.read(), about)

setup(
    name=about['__title__'],

    version=about['__version__'],

    description=about['__summary__'],
    long_description=LONG_DESC,

    url=about['__url__'],

    author=about['__author__'],
    author_email=about['__email__'],

    license=about['__license__'],

    # exclude tests from built/installed package
    packages=find_packages(exclude=['*.tests', '*.tests.*', 'docs/presentations',
                                    'docs/presentations/*']),
    package_data={
        'iCount': [
            'examples/*.sh',
        ]
    },
    install_requires={
        'numpy',
        'pandas',
        'cutadapt>=1.10',
        'pysam',
        'pybedtools',
        'numpydoc',
        'sphinx>=1.4',
        'matplotlib',
    },
    extras_require={
        'docs': [
            'docutils',
            'releases',
            'sphinx_rtd_theme',
        ],
        'package': [
            'pypandoc'
            'twine',
            'wheel',
        ],
        'test': [
            'check-manifest',
            'pylint>=1.6.4',
            'pycodestyle>=2.1.0',
            'pydocstyle>=1.0.0',
            'pytest-cov',
            'readme_renderer',
            'coverage>=4.2',
        ],
    },

    entry_points={
        'console_scripts': [
            'iCount = iCount.cli:main',
        ],
    },

    classifiers=[
        'Development Status :: 4 - Beta',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Operating System :: POSIX',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
    ],

    keywords='iCLIP protein-RNA',
)

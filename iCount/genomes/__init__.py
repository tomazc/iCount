"""
Genomes
=======

This module provides access to `Ensembl`_ genome sequence and annotation. Segmentation into genes
and into segments of same type (exon, intron, UTR, ...) is supported.

.. automodule:: iCount.genomes.ensembl
   :members:

.. automodule:: iCount.genomes.segment
   :members:

.. _Ensembl:
    http://www.ensembl.org/index.html

"""

from . import ensembl
from . import segment
from . import genes

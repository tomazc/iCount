# pylint: disable=invalid-name
""".. Line to protect from pydocstyle D205, D400.

iCount: protein-RNA interaction analysis
========================================

iCount is a Python module and associated command-line interface (CLI), which provides all the
commands needed to process protein-RNA `iCLIP`_ interaction data and to identify and quantify
sites of protein-RNA interactions on RNA.

iCount's main input are FASTQ files with `iCLIP`_ sequencing data, its main output are BED files
with identified and quantified cross-linked sites.

A number of analyses are included in iCount that provide insights into the properties of
protein-RNA interaction.

.. _iCLIP: http://www.ulelab.info

"""
import os

from iCount.__about__ import __author__, __copyright__, __email__, __license__, __summary__, \
    __title__, __url__, __version__
from . import analysis
from . import externals
from . import files
from . import demultiplex
from . import mapping
from . import genomes
from . import examples
from .logger import log_to_stdout, log_to_file, log_inputs, _log_progress
from .metrics import Metrics

#: Default output folder for iCount functions/commands. It points to the value
#: of eviroment variable ``ICOUNT_OUTPUT_ROOT`` if set. Otherwise, current
#: working directory is used.
OUTPUT_ROOT = os.environ.get('ICOUNT_OUTPUT_ROOT', '.')

#: Default temporary folder for iCount functions/commands. It is used to store
#: temporary files, used in intermediate steps of analyses. It points to the value
#: of eviroment variable ``ICOUNT_TMP_ROOT`` if set. Otherwise, ``/tmp/iCount`` is
#: used.
TMP_ROOT = os.environ.get('ICOUNT_TMP_ROOT', '/tmp/iCount')

# Create folders if needed:
if not os.path.exists(OUTPUT_ROOT):
    print("OUTPUT_ROOT folder does not exist. Creating it at: {}".format(OUTPUT_ROOT))
    os.makedirs(OUTPUT_ROOT)
if not os.path.exists(TMP_ROOT):
    print("TMP_ROOT folder does not exist. Creating it at: {}".format(TMP_ROOT))
    os.makedirs(TMP_ROOT)

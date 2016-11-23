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

.. _iCLIP:
    http://www.ulelab.info

"""

import os

from iCount.__about__ import (
    __author__, __copyright__, __email__, __license__, __summary__, __title__,
    __url__, __version__,
)

from . import analysis
from . import externals
from . import files
from . import demultiplex
from . import mapping
from . import genomes

from . import examples
from .logger import log_to_stdout, log_to_file, log_inputs, _log_progress
from .metrics import Metrics

# CONFIG
#: Output path points to the current working folder by default.
OUTPUT_ROOT = '.'

#: Default output folder can be set via the environment variable specified here.
#: Note, output path may be overridden by optional parameters of individual commands.
ICOUNT_OUTPUT_ROOT_VAR = 'ICOUNT_OUTPUT_ROOT'
OUTPUT_ROOT = os.environ.get(ICOUNT_OUTPUT_ROOT_VAR, OUTPUT_ROOT)


#: Path to temporary folder.
TMP_ROOT = '/tmp/iCount'

#: Default temporary folder can be set via the environment variable specified here.
#: Note, temporary path may be overridden by optional parameters of individual commands.
ICOUNT_TMP_ROOT_VAR = 'ICOUNT_TMP_ROOT'
TMP_ROOT = os.environ.get(ICOUNT_TMP_ROOT_VAR, TMP_ROOT)


# create folders if needed
if not os.path.exists(OUTPUT_ROOT):
    print("Output root folder does not exist. Will create it at: "
          "{0:s}".format(OUTPUT_ROOT))
    os.makedirs(OUTPUT_ROOT)

if not os.path.exists(TMP_ROOT):
    print("tmp root folder does not exist. Will create it at: "
          "{0:s}".format(TMP_ROOT))
    os.makedirs(TMP_ROOT)

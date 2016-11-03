"""
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

from . import analysis
from . import externals
from . import files
from . import demultiplex
from . import mapping
from . import genomes
from ._version import __version__

from . import examples
from .logger import log_to_stdout, log_to_file, log_inputs, _log_progress
from .metrics import Metrics

# CONFIG
#: Output path points to the current working folder by default.
output_root = '.'

#: Default output folder can be set via the environment variable specified here.
#: Note, output path may be overridden by optional parameters of individual commands.
ICOUNT_OUTPUT_ROOT_VAR = 'ICOUNT_OUTPUT_ROOT'
output_root = os.environ.get(ICOUNT_OUTPUT_ROOT_VAR, output_root)


#: Path to temporary folder.
tmp_root = '/tmp/iCount'

#: Default temporary folder can be set via the environment variable specified here.
#: Note, temporary path may be overridden by optional parameters of individual commands.
ICOUNT_TMP_ROOT_VAR = 'ICOUNT_TMP_ROOT'
tmp_root = os.environ.get(ICOUNT_TMP_ROOT_VAR, tmp_root)


# create folders if needed
if not os.path.exists(output_root):
    print("Output root folder does not exist. Will create it at: "
          "{0:s}".format(output_root))
    os.makedirs(output_root)

if not os.path.exists(tmp_root):
    print("tmp root folder does not exist. Will create it at: "
          "{0:s}".format(tmp_root))
    os.makedirs(tmp_root)

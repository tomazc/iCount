"""iCount pipeline.

iCount is a Python library for processing iCLIP and other NGS data.

=======
Modules
=======

.. automodule:: iCount.analysis

.. automodule:: iCount.externals

.. automodule:: iCount.files

.. automodule:: iCount.genomes

.. automodule:: iCount.mapping

.. automodule:: iCount.demultiplex

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

# CONFIG
# Because iCount is used in command-line, all paths point to the current
# working folder by default.
genomes_root = '.'
output_root = '.'
tmp_root = '/tmp/iCount'

# Paths can be specified with environment variables.
ICOUNT_GENOMES_ROOT_VAR = 'ICOUNT_GENOMES_ROOT'
ICOUNT_OUTPUT_ROOT_VAR = 'ICOUNT_OUTPUT_ROOT'
genomes_root = os.environ.get(ICOUNT_GENOMES_ROOT_VAR, genomes_root)
output_root = os.environ.get(ICOUNT_OUTPUT_ROOT_VAR, output_root)

# Paths are overridden by optional parameters of individual commands.


# create folders if not present yet
if not os.path.exists(genomes_root):
    print("Genomes root folder does not exist. Will create it at: "
          "{0:s}".format(genomes_root))
    os.makedirs(genomes_root)

if not os.path.exists(output_root):
    print("Output root folder does not exist. Will create it at: "
          "{0:s}".format(output_root))
    os.makedirs(output_root)

if not os.path.exists(tmp_root):
    print("tmp root folder does not exist. Will create it at: "
          "{0:s}".format(tmp_root))
    os.makedirs(tmp_root)

"""
iCount is a Python library for processing iCLIP and other NGS data.
"""

import os
import inspect

from . import analysis
from . import externals
from . import files
from . import demultiplex
from . import mapping
from . import genomes
from ._version import __version__

from . import examples
from .logger import log_to_stdout, log_to_file, log_inputs

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


class Result:

    def __init__(self, context=None):
        # If context is not given, determine it from calling function.
        if context is None:
            previous_frame = inspect.getouterframes(inspect.currentframe())[1]
            module = inspect.getmodulename(previous_frame[1])
            context = module + '.' + previous_frame[3]
        self.context = context

    def __repr__(self):
        return '{:s}\n{:s}'.format(self.context, self.__dict__.__repr__())

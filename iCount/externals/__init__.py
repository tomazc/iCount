""".. Line to protect from pydocstyle D205, D400.

Externals
=========

A number of external software is needed for iCount to work.

.. automodule:: iCount.externals.cutadapt
   :members:

.. automodule:: iCount.externals.star
   :members:

"""

from . import cutadapt
from . import star

EXPECTED_CUTADAPT_VERSION = '1.10'

"""
GTF
---

Reading `GTF`_ files.

"""

import pybedtools

def load(fn):
    return pybedtools.BedTool(fn)

import os
import gzip

from . import bam
from . import bed
from . import fastq
from . import gtf
from . import temp


def gz_open(fname, omode):
    if fname.endswith(".gz"):
        base_fname = os.path.basename(fname)
        return gzip.GzipFile(base_fname, fileobj=open(fname, omode),
                             mode=omode)
    return open(fname, omode)

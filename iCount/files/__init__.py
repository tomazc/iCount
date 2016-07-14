from . import bam
from . import bed
from . import fastq
from . import gtf


def gz_open(fname, omode):
    if fname.endswith(".gz"):
        base_fname = os.path.basename(fname)
        return gzip.GzipFile(base_fname, fileobj=open(fname, omode), mode=omode)
        #return gzip.open(fname, omode)
    return open(fname, omode)

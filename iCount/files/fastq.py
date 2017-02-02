# pylint: skip-file
""".. Line to protect from pydocstyle D205, D400.

FASTQ
-----

Reading and writting `FASTQ`_ files.
"""
import gzip
import math
import os


# precompute, so mapping is faster
_dX2L = [x - 64 for x in range(59, 105)]
_dX2L = [10.0 * math.log(math.pow(10.0, x / 10.0) + 1.0, 10.0) for x in _dX2L]
_dX2L = [chr(math.floor(x) + 33) for x in _dX2L]

_Dx2l = {chr(ij): v for ij, v in zip(range(59, 105), _dX2L)}


def _x2l(s):
    """TODO."""
    # Q_solexa = [ord(x) - 64 for x in s]
    # Q_PHRED = [10.0 * math.log(math.pow(10.0, x / 10.0) + 1.0, 10.0) for x in Q_solexa]
    # L = [chr(math.floor(x) + 33) for x in Q_PHRED]
    # L = [_dX2L[ord(x) - 59] for x in x]
    # return ''.join(L)
    return ''.join((_Dx2l[ij] for ij in s))


# -31 = -64 + 33
_Dij2l = {chr(ij): chr(ij - 31) for ij in range(64, 104 + 1 + 1)}


def _ij2l(s):
    """TODO."""
    return ''.join((_Dij2l[ij] for ij in s))


# 0 = -33 + 33
_Dsl2l = {chr(ij): chr(ij - 0) for ij in range(64, 104 + 1 + 1)}


def _sl2l(s):
    """TODO."""
    # return ''.join((_Dsl2l[ij] for ij in s))
    return s


def _determine_format(fname):
    """
    Read first few records and determine quality encoding in FASTQ file.

    Set function for conversion of quality scores to Phred scores
    in format L. See format description
    `http://en.wikipedia.org/wiki/FASTQ_format`_.

    S - Sanger        Phred+33,  raw reads typically (0, 40)   [33..73]
    X - Solexa        Solexa+64, raw reads typically (-5, 40)  [59..104]
    I - Illumina 1.3+ Phred+64,  raw reads typically (0, 40)   [64..104]
    J - Illumina 1.5+ Phred+64,  raw reads typically (3, 40)   [66..104]
    with 0=unused, 1=unused, 2=Read Segment Quality Control Indicator (bold)
    (Note: See discussion above).
    L - Illumina 1.8+ Phred+33,  raw reads typically (0, 41)   [33..74]

    """
    # Read a few reads to determine format of quality scores.

    r2q = {
        (33, 73): 'S',
        (59, 104): 'X',
        (64, 104): 'I',
        (66, 104): 'J',
        (33, 74): 'L',
    }

    # During scan assume that index is in header line.
    format = 'L'
    quals = set()
    for rec_cn, _, r in _reader(fname, _sl2l):
        quals.update(set(r.r_qual))

        if rec_cn % 10000 == 0:
            minq = ord(min(quals))
            maxq = ord(max(quals))

            format = r2q.get((minq, maxq), None)
            if format is not None:
                return format
    return '?'


def reader(fname, check=False):
    """TODO."""
    if check:
        format = _determine_format(fname)
    else:
        format = 'L'
        qual_xform = _sl2l

    assert(format in ['L', 'S', 'X', 'I', 'J'])
    if format == 'L':
        qual_xform = _sl2l
    elif format == 'S':
        qual_xform = _sl2l
    elif format == 'X':
        qual_xform = _x2l
    elif format == 'I':
        qual_xform = _ij2l
    elif format == 'J':
        qual_xform = _ij2l

    return _reader(fname, qual_xform)


def is_gzip(fname):
    """TODO."""
    try:
        f = gzip.open(fname, 'r')
        f.peek(10)
        return True
    except (OSError, AttributeError):
        return False


class Read:
    """TODO."""

    def __init__(self, r_id, r_seq, r_plus, r_qual):
        """TODO."""
        self.r_id = r_id
        self.r_seq = r_seq
        self.r_plus = r_plus
        self.r_qual = r_qual

    def __str__(self):
        """TODO."""
        return '\n'.join(map(str, [self.r_id, self.r_seq, self.r_plus, self.r_qual])) + '\n'


def read_one(fin):
    """TODO."""
    # r_id = fin.readline().decode('ascii').rstrip('\n')
    # r_seq = fin.readline().decode('ascii').rstrip('\n')
    # r_plus = fin.readline().decode('ascii').rstrip('\n')
    # r_qual = fin.readline().decode('ascii').rstrip('\n')
    r_id = str(fin.readline()).rstrip('\n')
    r_seq = str(fin.readline()).rstrip('\n')
    r_plus = str(fin.readline()).rstrip('\n')
    r_qual = str(fin.readline()).rstrip('\n')
    print("AAAA", r_id, r_seq, r_plus, r_qual)
    return Read(r_id, r_seq, r_plus, r_qual)


def _reader(fname, qual_xform):
    """TODO."""
    if is_gzip(fname):
        f = gzip.open(fname, 'rt')
    else:
        f = open(fname, 'rt')

    # get file size
    f_size = os.path.getsize(fname)
    if f_size == 0:
        f_size = 1

    r_id = f.readline()
    rec_cn = 1
    per = 0

    while r_id:
        r_id = r_id.rstrip('\n')
        r_seq = f.readline().rstrip('\n')
        r_plus = f.readline().rstrip('\n')
        r_qual = f.readline().rstrip('\n')
        r_qualL = qual_xform(r_qual)

        if rec_cn % 1000000 == 0:
            f_cur_pos = f.tell()
            per = min(100, round(100 * f_cur_pos / f_size))

        yield rec_cn, per, Read(r_id, r_seq, r_plus, r_qualL)
        r_id = f.readline()
        rec_cn += 1


class Writer:
    """TODO."""

    def __init__(self, fname):
        """TODO."""
        if not fname:
            return

        if fname.endswith('.gz'):
            self.f = gzip.open(fname, 'wt')
        else:
            self.f = open(fname, 'wt')

    def __del__(self):
        """TODO."""
        self.close()

    def write(self, r_id, r_seq, r_plus, r_qual):
        """TODO."""
        self.f.write('\n'.join(map(str, [r_id, r_seq, r_plus, r_qual])) + '\n')

    def close(self):
        """TODO."""
        if self.f:
            self.f.close()
            self.f = None

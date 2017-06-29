""".. Line to protect from pydocstyle D205, D400.

FASTQ
-----

Reading and writting `FASTQ`_ files.
"""
import os
import math

import iCount


# From here: https://en.wikipedia.org/wiki/FASTQ_format#Quality
# one can deduce the following equation, connecting Q_solexa with Q_sanger:
# 10**(Q_solexa/10) = 10**(Q_sanger/10) - 1
# and therefore:
# Q_sanger = 10 * log_10(10**(Q_solexa/10) + 1)
# For all possible Q_solexa values (-5 to 40), here are corresponding Q_sanger values
X2L_ = [10.0 * math.log(math.pow(10.0, x / 10.0) + 1.0, 10.0) for x in range(-5, 41)]
# Mapping from Solexa quality encoding to Illumina 1.8+
X2L = {chr(i): chr(math.floor(v) + 33) for i, v in zip(range(59, 105), X2L_)}

# Mapping from Illumina 1.3+/1.5+ quality encoding to Illumina 1.8+:
I2L = {chr(i): chr(i - 31) for i in range(64, 105)}


def _x2l(quals):
    """Change Solexa quality encoding to Illumina 1.8+."""
    return ''.join([X2L[qual] for qual in quals])


def _i2l(quals):
    """Change Illumina 1.3+/1.5+ quality encoding to Illumina 1.8+."""
    return ''.join([I2L[qual] for qual in quals])


def _l2l(quals):
    """Change Sanger quality encoding to Illumina 1.8+."""
    return quals


def _transform_encoding(encoding):
    """Transform any quality encoding to Illumina 1.8+."""
    assert(encoding in ['S', 'X', 'I', 'J', 'L'])

    if encoding == 'S':
        return _l2l
    elif encoding == 'X':
        return _x2l
    elif encoding == 'I':
        return _i2l
    elif encoding == 'J':
        return _i2l
    elif encoding == 'L':
        return _l2l


def _get_qual_encoding(fname):
    """
    Read first few records and determine quality encoding in FASTQ file.

    See format description
    `http://en.wikipedia.org/wiki/FASTQ_format`_.

    S - Sanger        Phred+33,  raw reads typically (0, 40)   [33..73]
    X - Solexa        Solexa+64, raw reads typically (-5, 40)  [59..104]
    I - Illumina 1.3+ Phred+64,  raw reads typically (0, 40)   [64..104]
    J - Illumina 1.5+ Phred+64,  raw reads typically (3, 40)   [66..104]
    L - Illumina 1.8+ Phred+33,  raw reads typically (0, 41)   [33..74]
    """
    def get_encoding(quals, count=0, check_count=False):
        """Determine encoding from quality scores raneg."""
        minq, maxq = ord(min(quals)), ord(max(quals))
        if minq < 59:
            if maxq > 73:
                return 'L'
            else:
                return 'S'
        elif 59 <= minq < 64 and maxq > 74:
            return 'X'
        elif minq >= 64 and maxq > 74:
            if (check_count and count > 10000) or not check_count:
                if minq < 66:
                    return 'I'
                else:
                    return 'J'

    quals = set()
    for count, read in enumerate(FastqFile(fname).read(encoding='L')):
        quals.update(set(read.qual))
        if count % 10000 == 0:
            encoding = get_encoding(quals, count, check_count=True)
            if encoding:
                return encoding

    encoding = get_encoding(quals, check_count=False)
    return encoding


class FastqEntry:
    """Single FASTQ entry."""

    def __init__(self, id, seq, plus, qual):  # pylint: disable=redefined-builtin
        """Initialize attributes."""
        self.id = str(id).split(' ')[0]  # pylint: disable=invalid-name
        self.seq = str(seq)
        self.plus = str(plus)
        self.qual = str(qual)

    def __repr__(self):
        """Represent object."""
        return self.id


class FastqFile:
    """Write FASTQ files."""

    def __init__(self, fname, mode='rt'):
        """Open file handle in desired mode."""
        self.fname = fname
        if 'r' in mode and not os.path.isfile(fname):
            self.file = None
            raise FileNotFoundError('File not found.')

        self.file = iCount.files.gz_open(fname, mode)

    def __del__(self):
        """Close file."""
        self.close()

    def read(self, encoding=None):
        """Read FASTQ file."""
        if not encoding:
            encoding = _get_qual_encoding(self.fname)
        # Function used for quality transform:
        qual_xform = _transform_encoding(encoding)

        for read_id in self.file:
            read_seq = next(self.file).rstrip('\n')
            read_plus = next(self.file).rstrip('\n')
            read_qual = qual_xform(next(self.file).rstrip('\n'))
            yield FastqEntry(read_id.rstrip('\n'), read_seq, read_plus, read_qual)

    def write(self, fq_entry):
        """Write single FASTQ entry."""
        content = [fq_entry.id, fq_entry.seq, fq_entry.plus, fq_entry.qual]
        self.file.write('\n'.join(map(str, content)) + '\n')

    def close(self):
        """Close file if it is stil open."""
        if self.file and not self.file.closed:
            self.file.close()

""".. Line to protect from pydocstyle D205, D400.

FASTQ
-----

Reading and writting `FASTQ`_ files.
"""
import os

import iCount

ENCODING_TO_OFFSET = {
    'S': 33,
    'X': 64,
    'I': 64,
    'J': 64,
    'L': 33,
}


def get_qual_encoding(fname):
    """
    Read first few records and determine quality encoding in FASTQ file.

    See format description:
    `http://en.wikipedia.org/wiki/FASTQ_format`

    S - Sanger        Phred+33,  raw reads typically (0, 40)   [33..73]
    X - Solexa        Solexa+64, raw reads typically (-5, 40)  [59..104]
    I - Illumina 1.3+ Phred+64,  raw reads typically (0, 40)   [64..104]
    J - Illumina 1.5+ Phred+64,  raw reads typically (3, 40)   [66..104]
    L - Illumina 1.8+ Phred+33,  raw reads typically (0, 41)   [33..74]
    """
    def get_encoding(quals, count=0, check_count=False):
        """Determine encoding from quality scores range."""
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
    for count, read in enumerate(FastqFile(fname).read()):
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

    def read(self):
        """Read FASTQ file."""
        for read_id in self.file:
            read_seq = next(self.file).rstrip('\n')
            read_plus = next(self.file).rstrip('\n')
            read_qual = next(self.file).rstrip('\n')
            yield FastqEntry(read_id.rstrip('\n'), read_seq, read_plus, read_qual)

    def write(self, fq_entry):
        """Write single FASTQ entry."""
        content = [fq_entry.id, fq_entry.seq, fq_entry.plus, fq_entry.qual]
        self.file.write('\n'.join(map(str, content)) + '\n')

    def close(self):
        """Close file if it is stil open."""
        if self.file and not self.file.closed:
            self.file.close()

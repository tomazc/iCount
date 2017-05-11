"""Utility functions for testing."""
# pylint: disable=protected-access
import os
import tempfile

import pysam

from numpy import random

import pybedtools
from pybedtools import create_interval_from_list

import iCount

BASES = ['A', 'C', 'G', 'T']


def get_temp_file_name(tmp_dir=None, extension=''):
    """Return an availiable name for temporary file."""
    tmp_name = next(tempfile._get_candidate_names())
    if not tmp_dir:
        tmp_dir = tempfile._get_default_tempdir()
    if extension is not None:
        tmp_name = tmp_name + '.' + extension
    return os.path.join(tmp_dir, tmp_name)


def get_temp_dir(tmp_dir=None):
    """Return a temporary directory."""
    return tempfile.mkdtemp()


def make_file_from_list(data, bedtool=True, extension=''):
    """Return path to file with the content from `data` (list of lists)."""
    tfile = get_temp_file_name(extension=extension)
    if bedtool:
        tmp_bt = pybedtools.BedTool(
            pybedtools.create_interval_from_list(line) for line in data).saveas()
        tmp_bt.saveas(tfile)
    else:
        with open(tfile, 'wt') as file_:
            for list_ in data:
                file_.write('\t'.join(map(str, list_)) + '\n')
    return os.path.abspath(tfile)


def make_list_from_file(fname, fields_separator=None):
    """Read file to list of lists."""
    data = []
    with iCount.files.gz_open(fname, 'rt') as file_:
        for line in file_:
            data.append(line.strip().split(fields_separator))
    return data


def list_to_intervals(data):
    """Transform list of lists to list of pybedtools.Intervals."""
    return [create_interval_from_list(list_) for list_ in data]


def intervals_to_list(data):
    """Transform list of pybedtools.Intervals to list of lists."""
    return [interval.fields for interval in data]


def reverse_strand(data):
    """
    Reverse the strand of every element in data.

    Elements can be pybedtools.Intervals or lists with same content as
    interval.fields.
    """
    rstrands = ['-' if i[6] == '+' else '+' for i in data]
    if isinstance(data[0], pybedtools.Interval):
        return [create_interval_from_list(
            data[i][:6] + [rstrands[i]] + data[i][7:]) for i in range(len(data))]
    else:
        return [
            data[i][:6] + [rstrands[i]] + data[i][7:] for i in range(len(data))]


def make_sequence(size, include_n=False):
    """Make random DNA segment of length `size`."""
    if include_n:
        bases = ['A', 'C', 'G', 'T', 'N']
    else:
        bases = ['A', 'C', 'G', 'T']
    return ''.join(random.choice(bases, size))  # pylint: disable=no-member


def make_quality_scores(size, min_chr=33, max_chr=127):
    """Make random DNA segment of length `size`."""
    scores = [chr(i) for i in range(min_chr, max_chr + 1)]
    return ''.join(random.choice(scores, size))  # pylint: disable=no-member


def make_aligned_segment(data):
    """
    Return pysam.AlignedSegment() object with the data in `data`.

    Each pysam.AlignedSegment element has 11 mandatory tab-separated
    fields. Additionaly there is TAGS field for additional info::

        QNAME FLAG RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN SEQ QUAL TAGS

    Since only some fields are required for our purpuses, the input
    parameter `data` should be a tuple with the following content:

    (qname, flag, refname, pos, mapq, cigar, tags)

        * qname - querry name
        * flag - bitwise flag, detailed description in [1]
        * refname - Index of the reference_sequence in header
        * pos - starting position of read in reference sequence
        * mapq - mapping quality
        * cigar - CIGAR values, detailed description in [1]
        * tags - additional information

    Flag is a bitwise value. For detailed explanation, see [1], but for our
    needs, we only need values 4 (is_unmapped) and 16 (reverse_strand).

    Columns RNEXT and RNAME are left undefined, since we are only interested
    in sigle-end reads. TLEN is also set to undefined (value 0). SEQ and
    QUAL are randomly generated, their size is determined by cigar value.

    Example of input parameter `data`:

        data = ('name123', 20, 0, 30, 20, [(0, 90), (2, 5), (0, 5)], {'NH': 5})

    Read more about SAM file specifications here:
    [1] https://samtools.github.io/hts-specs/SAMv1.pdf

    Parameters
    ----------
    data : tuple
        Input data for AlignedSegment.

    Returns
    -------
    pysam.AlignedSegment
        AlignedSegment, filled with given data.

    """
    # pylint: disable=no-member
    segment = pysam.AlignedSegment()
    segment.query_name = data[0]
    segment.flag = data[1]
    segment.reference_id = data[2]
    segment.reference_start = data[3]
    segment.mapping_quality = data[4]
    segment.cigar = data[5]

    segment.next_reference_id = 0
    segment.next_reference_start = 0
    segment.template_length = 0

    length = sum([n2 for (n1, n2) in segment.cigar])
    segment.query_sequence = make_sequence(size=length, include_n=True)
    segment.query_qualities = pysam.qualitystring_to_array(
        make_quality_scores(size=length))

    segment.tags = data[6].items()

    return segment


def make_bam_file(data):
    """
    Make a BAM file and fill it with content of `data`.

    Each BAM file consists of header and allignment entries. Main content of the
    header are the reference sequence definitions (chromosomes). Allignment
    entries (segments) are records of how each read is mapped to a reference.

    The structure of input parameter `data` should be the following:

        data = {
            # define all the chromosome names and their lengths:
            'chromosomes': [('chr1', 3000), ('chr2', 2000)],
            'segments': [
                # Each segment should be a 7-tuple - for detailed description
                # see `make_aligned_segment` function docstring
                # (qname, flag, refname, pos, mapq, cigar, tags)
                ('name1', 4, 0, 100, 20, [(0, 100)], {'NH': 2}),
                ...
                ('namex', 4, 0, 100, 20, [(0, 100)], {'NH': 2})]}

    Links for further read:
    https://samtools.github.io/hts-specs/SAMv1.pdf
    http://pysam.readthedocs.io/en/stable/usage.html#creating-bam-cram-sam-files-from-scratch

    Parameters
    ----------
    data : dict
        Input data for BAM file

    Returns
    -------
    str
        Absoulte path to bamfile.

    """
    fname = get_temp_file_name(extension='.bam')

    # Make header:
    chromosomes = [dict(
        [('SN', chrom), ('LN', size)]) for chrom, size in data['chromosomes']]
    header = {'HD': {'VN': '1.0'}, 'SQ': chromosomes}

    with pysam.AlignmentFile(fname, "wb", header=header) as outf:  # pylint: disable=no-member
        for segment_data in data['segments']:
            outf.write(make_aligned_segment(segment_data))

    return os.path.abspath(fname)


def make_fasta_file(headers=None, out_file=None, num_sequences=10, seq_len=80):
    """Make artificial FASTA file."""
    if headers is None:
        headers = ['{}'.format(i + 1) for i in range(num_sequences)]
    else:
        num_sequences = len(headers)

    if out_file is None:
        out_file = get_temp_file_name()
    with open(out_file, 'wt') as ofile:
        for header in headers:
            ofile.write('>' + header + '\n')
            ofile.write(make_sequence(seq_len) + '\n')

    return os.path.abspath(out_file)


def make_fastq_file(genome=None, barcodes=None, adapter='', out_file=None,
                    num_sequences=10, seq_len=80):
    """
    Make artificial FASTQ file.

    TODO: refactor and write more descriptive doscrting.
    """
    if barcodes is None:
        barcodes = ['NNN']

    if genome is not None:
        genome_data = iCount.files.fasta.read_fasta(genome)
        num_sequences = len(genome_data)
        seq_len = len(genome_data[0])

    seq_len_reduced = seq_len - (len(barcodes[0]) if len(barcodes) != 0 else 0) - len(adapter)

    num_barcodes = len(barcodes)

    def make_fastq_entry(ofile):  # pylint: disable=missing-docstring
        description = 'artificial header {}'.format(random.random())

        # Select random barcode and transform 'N'-s to nucleotids:
        pre_barcode = barcodes[random.randint(0, num_barcodes)]  # pylint: disable=no-member
        barcode = ''
        for letter in pre_barcode:
            if letter == 'N':
                # Randomly select one of bases:
                barcode += BASES[random.randint(0, 4)]  # pylint: disable=no-member
            else:
                barcode += letter

        if genome:
            int1 = random.randint(0, high=len(genome_data))  # pylint: disable=no-member
            int2 = random.randint(0, high=len(genome_data[int1][1]))  # pylint: disable=no-member
            random_piece = genome_data[int1][1][int2:int2 + seq_len_reduced]
            seq = barcode + random_piece + adapter
        else:
            seq = barcode + make_sequence(seq_len_reduced) + adapter
        quality_scores = make_quality_scores(len(seq))

        ofile.write('@' + description + '\n')
        ofile.write(seq + '\n')
        ofile.write('+' + '\n')
        ofile.write(quality_scores + '\n')

    if out_file is None:
        out_file = get_temp_file_name()
    with open(out_file, 'wt') as ofile:
        for _ in range(num_sequences):
            make_fastq_entry(ofile)

    return os.path.abspath(out_file)

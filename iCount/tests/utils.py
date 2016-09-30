import os
import tempfile

import pysam
import pybedtools

from numpy import random
from pybedtools import create_interval_from_list

from iCount.files.fasta import read_fasta


BASES = ['A', 'C', 'G', 'T']


def get_temp_file_name(tmp_dir=None, extension=''):
    """
    Returns an availiable name for temporary file.
    """
    tmp_name = next(tempfile._get_candidate_names())
    if not tmp_dir:
        tmp_dir = tempfile._get_default_tempdir()
    if extension is not None:
        tmp_name = tmp_name + '.' + extension
    return os.path.join(tmp_dir, tmp_name)


def get_temp_dir(tmp_dir=None):
    """
    Returns a temporary directory.
    """
    return tempfile.mkdtemp()


def make_file_from_list(data, bedtool=True):
    """Return path to file with the content from `data` (list of lists)"""
    tfile = get_temp_file_name()
    if bedtool:
        pybedtools.BedTool(pybedtools.create_interval_from_list(list_)
                           for list_ in data).saveas(tfile)
    else:
        with open(tfile, 'wt') as file_:
            for list_ in data:
                file_.write('\t'.join(map(str, list_)) + '\n')
    return os.path.abspath(tfile)


def make_list_from_file(fname, fields_separator=None):
    """Read file to list of lists"""
    data = []
    with open(fname) as file_:
        for line in file_:
            data.append(line.strip().split(fields_separator))
    return data


def list_to_intervals(data):
    """Transform list of lists to list of pybedtools.Intervals."""
    return [create_interval_from_list(list_) for list_ in data]


def intervals_to_list(data):
    """Transform list of pybedtools.Intervals to list of lists"""
    return [interval.fields for interval in data]


def reverse_strand(data):
    """Reverses the strand of every element in data. Elements can be
    pybedtools.Intervals or lists with same content as interval.fields."""
    rstrands = ['-' if i[6] == '+' else '+' for i in data]
    if isinstance(data[0], pybedtools.Interval):
        return [create_interval_from_list(
            data[i][:6] + [rstrands[i]] + data[i][7:]) for i in range(len(data))]
    else:
        return [
            data[i][:6] + [rstrands[i]] + data[i][7:] for i in range(len(data))]


def make_sequence(size, include_N=False):
    """
    Makes random DNA segment of length `size`
    """
    if include_N:
        bases = ['A', 'C', 'G', 'T', 'N']
    else:
        bases = ['A', 'C', 'G', 'T']
    return ''.join(random.choice(bases, size))


def make_quality_scores(size, min_chr=33, max_chr=127):
    """
    Makes random DNA segment of length `size`
    """
    scores = [chr(i) for i in range(min_chr, max_chr + 1)]
    return ''.join(random.choice(scores, size))


def make_aligned_segment(data):
    """
    Return pysam.AlignedSegment() object with the data in `data`

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
    a = pysam.AlignedSegment()
    a.query_name = data[0]
    a.flag = data[1]
    a.reference_id = data[2]
    a.reference_start = data[3]
    a.mapping_quality = data[4]
    a.cigar = data[5]

    a.next_reference_id = 0
    a.next_reference_start = 0
    a.template_length = 0

    length = sum([n2 for (n1, n2) in a.cigar])
    a.query_sequence = make_sequence(size=length, include_N=True)
    a.query_qualities = pysam.qualitystring_to_array(
        make_quality_scores(size=length))

    a.tags = data[6].items()

    return a


def make_bam_file(data):
    """
    Make a BAM file and fill it with content of `data`

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
    fname = get_temp_file_name()

    # Make header:
    chromosomes = [dict(
        [('SN', chrom), ('LN', size)]) for chrom, size in data['chromosomes']]
    header = {'HD': {'VN': '1.0'}, 'SQ': chromosomes}

    with pysam.AlignmentFile(fname, "wb", header=header) as outf:
        for segment_data in data['segments']:
            outf.write(make_aligned_segment(segment_data))

    return os.path.abspath(fname)


def make_fasta_file(headers=None, out_file=None, num_sequences=10, seq_len=80):
    """
    Make artificial FASTA file.

    Headers should be a list with length equal num_sequences.
    """
    if headers is None:
        headers = ['{}'.format(i + 1) for i in range(num_sequences)]

    def make_fasta_entry(ofile, headers):
        ofile.write('>' + headers.pop(0) + '\n')
        ofile.write(make_sequence(seq_len) + '\n')

    if out_file is None:
        out_file = get_temp_file_name()
    with open(out_file, 'wt') as ofile:
        for i in range(num_sequences):
            make_fasta_entry(ofile, headers)

    return os.path.abspath(out_file)


def make_fastq_file(genome=None, barcodes=None, adapter='', out_file=None,
                    num_sequences=10, seq_len=80):
    """
    Make artificial FASTQ file.

    File can be customized to with given barcodes and
    """
    if barcodes is None:
        barcodes = ['NNN']

    if genome is not None:
        genome_data = read_fasta(genome)
        num_sequences = len(genome_data)
        seq_len = len(genome_data[0])

    seq_len_reduced = seq_len - (len(barcodes[0]) if len(barcodes) != 0 else 0) - len(adapter)

    qualities = [chr(i) for i in range(33, 127)]
    num_barcodes = len(barcodes)

    def make_fastq_entry(ofile):
        description = 'artificial header {}'.format(random.random())

        # Select random barcode and transform 'N'-s to nucleotids:
        pre_barcode = barcodes[random.randint(0, num_barcodes)]
        barcode = ''
        for letter in pre_barcode:
            if letter == 'N':
                # Randomly select one of bases:
                barcode += BASES[random.randint(0, 4)]
            else:
                barcode += letter

        if genome:
            r1 = random.randint(0, high=len(genome_data))
            r2 = random.randint(0, high=len(genome_data[r1][1]))
            random_piece = genome_data[r1][1][r2:r2 + seq_len_reduced]
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
        for i in range(num_sequences):
            make_fastq_entry(ofile)

    return os.path.abspath(out_file)

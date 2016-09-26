"""
Demultiplexing
--------------

Split FASTQ file into separate FASTQ files, one for each sample barcode.
Saved FASTQ files contain sequences where sample barcode, random
barcode, and adapter sequences were removed. Random barcode is moved into
the header line, because needed in later steps, when removing PCR duplicates
and counting number of cross-link events.

"""


import os
import logging

import iCount


LOGGER = logging.getLogger(__name__)


def run(fastq_fn, adapter, mismatches, barcodes,
        minimum_length=15, prefix='demux', outdir='.'):
    """
    Demultiplex FASTQ file.

    Split input FASTQ file into separate files, one for each barcode, and
    additional file for non-matching barcodes.

    Parameters
    ----------
    fastq_fn : str
        Path to reads from a sequencing library.
    adapter : str
        Adapter sequence to remove from ends of reads.
    mismatches : int
        Number of tolerated mismatches when comparing barcodes.
    barcodes : list_str
        List of barcodes used for library.
    minimum_length : int
        Minimum length of trimmed sequence to keep.
    prefix : str
        Prefix of generated FASTQ files.
    outdir : str
        Output folder. Use local folder if none given.

    Returns
    -------
    str
        List of filenames where separated reads are stored.

    """
    iCount.log_inputs(LOGGER, level=logging.INFO)

    if not os.path.isdir(outdir):
        raise FileNotFoundError('Output directory does not exist. Make sure it does.')
    out_fn_prefix = []
    for bc in ['nomatch'] + barcodes:
        out_fn_prefix.append(
            os.path.join(outdir, '{:s}_{:s}'.format(prefix, bc))
        )

    if adapter:
        # need to remove adapter
        out_fns = ['{:s}_raw.fastq.gz'.format(fn) for fn in out_fn_prefix]
    else:
        out_fns = ['{:s}.fastq.gz'.format(fn) for fn in out_fn_prefix]

    # demultiplex
    LOGGER.info('Allowing max %d mismatches in barcodes.', mismatches)
    LOGGER.info('Demultiplexing file: %s', fastq_fn)
    demultiplex(fastq_fn, out_fns[1:], out_fns[0], barcodes, mismatches,
                minimum_length)
    LOGGER.info('Saving results to:')
    for fn in out_fns:
        LOGGER.info('    %s', fn)

    LOGGER.info('Trimming adapters (discarding shorter than %d)...', minimum_length)
    # remove adapter, if requested
    if adapter:
        out_fns_intermediate = out_fns
        out_fns = ['{:s}.fastq.gz'.format(fn) for fn in out_fn_prefix]
        for fn_in, fn_out in zip(out_fns_intermediate, out_fns):
            remove_adapter(fn_in, fn_out, adapter, minimum_length=minimum_length)
            os.remove(fn_in)

    return out_fns


def demultiplex(in_fastq_fname, out_fastq_fnames, not_matching_fastq_fname,
                barcodes, mismatches=1, minimum_length=15):
    """Extract reads and save to individual FASTQ files.

    All non-matching reads are stored in FASTQ file not_matching_fastq_fname.

    Parameters
    ----------
    in_fastq_fname : str
        Filename to read,  must be FASTQ in order to avoid duplicate read records.
    out_fastq_fnames : list_str
        List of filenames to store reads as determined by barcodes.
    not_matching_fastq_fname : str
        Fastq filename where to store non matching reads.
    barcodes : list_str
        Experiment and randomer barcode definition.
    mismatches : int
        Number of allowed mismatches in sample barcode.
    minimum_length : int
        Discard reads with length less than this.

    Returns
    -------
    str
        List of filenames where reads are stored (should be same as
        out_fastq_fnames).

    """

    out_fastq = [iCount.files.fastq.Writer(fn) for fn in
                 out_fastq_fnames + [not_matching_fastq_fname]]
    reader = iCount.files.fastq.Reader(in_fastq_fname)
    for r_id, exp_id, r_randomer, r_seq, r_plus, r_qual in \
            extract(reader, barcodes, mismatches=mismatches,
                    minimum_length=minimum_length):
        if r_randomer:
            if r_id[-2] == '/':
                r_pair = r_id[-2:]
                r_id = '%s:rbc:%s%s' % (r_id[:-2], r_randomer, r_pair)
            else:
                r_id = '%s:rbc:%s' % (r_id, r_randomer)

        out_fastq[exp_id].write(r_id, r_seq, '+', r_qual)
    while out_fastq:
        f = out_fastq.pop()
        f.close()
    return out_fastq_fnames


def extract(seqs, barcodes, mismatches=1, minimum_length=15):
    """Return iterator that returns experiment, randomer, and remaining
    sequence.

    :param seqs: input sequence of tuples (r_id, r_seq, r_qualL)
    :param barcodes: dictionary of experiment annotation
    :param mismatches: number of allowed mismatches in sample barcode
    :return: iterator
    """
    valid_nucs = {'A', 'T', 'C', 'G'}

    all_barcodes = len(barcodes)
    # at each barcode position, map nucleotide to corresponding exp_id
    p2n2i = {}
    i2start_pos = {}
    for bi, bar5 in enumerate(barcodes):
        for p, n in enumerate(bar5):
            if n == 'N':
                ns = ['A', 'T', 'C', 'G']
            else:
                ns = [n]
            assert set(ns) & valid_nucs == set(ns)
            for n2 in ns:
                p2n2i.setdefault(p, {}).setdefault(n2, set()).add(bi)
        assert i2start_pos.setdefault(bi, len(bar5)) == len(bar5)

    # determine the randomer positions - those that cannot be distinguished
    # based on barcodes
    random_pos = []
    p2i = []
    max_votes = 0
    for p, n2i in sorted(p2n2i.items()):
        if all(len(i) == all_barcodes for i in n2i.values()):
            random_pos.append(p)
        else:
            p2i.append((p, n2i))
            max_votes += 1

    # process sequences
    for r_cn, r_per, r in seqs:
        votes = [0]*all_barcodes
        for p, n2i in p2i:
            bis = n2i.get(r.r_seq[p], [])
            for bi in bis:
                votes[bi] += 1
        mvotes = max(votes)
        if mvotes < max_votes - mismatches:  # allowed mismatches
            # not recognized
            exp_id = -1

            r_randomer = ''
            r_seq = r.r_seq
            r_qual = r.r_qual
        else:
            # recognized as valid barcode
            exp_id = votes.index(mvotes)

            r_randomer = ''.join(r.r_seq[p] for p in random_pos)
            ps = i2start_pos[exp_id]
            r_seq = r.r_seq[ps:]
            r_qual = r.r_qual[ps:]

        if len(r_seq) >= minimum_length:
            yield r.r_id, exp_id, r_randomer, r_seq, r.r_plus, r_qual


def remove_adapter(in_fastq, out_fastq, adapter, qual_trim=10,
                   minimum_length=15):
    """

    :param in_fastq: FASTQ file with sequences
    :param out_fastq: FASTQ file where to save adapter-trimmed sequences
    :param adapter: sequence to be found at 3'
    :param qual_trim: quality cutoff threshold (see cutadapt manual)
    :return: exit status of cutadapt command line
    """

    return iCount.externals.cutadapt.run(in_fastq, out_fastq, adapter,
                                         qual_trim=qual_trim,
                                         minimum_length=minimum_length)

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
import iCount.files
import iCount.externals

# description and parameters needed for the analysis
analysis_name = 'demultiplex'
analysis_description_short = 'demultiplex FASTQ file'
analysis_description = 'Split input FASTQ file into separate files, one for ' \
                       'each barcode, and additional file for non-matching ' \
                       'barcodes.'

params_opt = [
    (
        'barcodes', 'str_list', [], True,
        'List of barcodes used for library.'
    ),
    (
        'adapter', 'string', 'AGATCGGAAGAGCGGTTCAG', False,
        'Adapter sequence to remove from ends of reads.'
    ),
    (
        'mismatches', 'int_range', (1, 0, 5), False,
        'Number of tolerated mismatches when comparing barcodes.'
    ),
    (
        'minimum_length', 'int_range', (15, 0, 100), False,
        'Minimum length of trimmed sequence to keep.'
    ),
    (
        'prefix', 'string', 'exp', False,
        'Prefix of generated FASTQ files.'
    ),
    (
        'outdir', 'string', '.', False,
        'Output folder. Use local folder if none given.'
    ),
]

params_pos = [
    (
        'sequences', 'FASTQ', 'in', 1,
        '(input) sequences from a sequencing library.'
    ),
]


def run(fastq_fn, barcodes, adapter, mismatches,
        minimum_length=15, prefix='demux', outdir='.'):
    assert os.path.isdir(outdir)
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
    demultiplex(fastq_fn, out_fns[1:], out_fns[0], barcodes, mismatches,
                minimum_length)

    # remove adapter, if requested
    if adapter:
        out_fns_intermediate = out_fns
        out_fns = ['{:s}.fastq.gz'.format(fn) for fn in out_fn_prefix]
        for fn_in, fn_out in zip(out_fns_intermediate, out_fns):
            iCount.externals.cutadapt.run(fn_in, fn_out, adapter,
                                          minimum_length=minimum_length)
            os.remove(fn_in)

    return out_fns


def demultiplex(in_fastq_fname, out_fastq_fnames, not_matching_fastq_fname,
                barcodes, mismatches=1, minimum_length=15):
    """Extract reads and save to individual FASTQ files.

    All non-matching reads are stored in FASTQ file not_matching_fastq_fname.

    :param in_fastq_fname: filename to read,  must be FASTQ in order to
    avoid duplicate read records.
    :param out_fastq_fnames: list of filenames to store reads as determined
    by barcodes
    : not_matching_fastq_fname: fastq filename where to store reads not
    matching
    :param not_matching_fastq_fname: FASTQ to store reads non matchning reads
    :param barcodes: experiment and randomer barcode definition
    :param mismatches: number of allowed mismatches in sample barcode
    """

    out_fastq = [iCount.files.fastq.Writer(fn) for fn in
                 out_fastq_fnames + [not_matching_fastq_fname]]
    reader = iCount.files.fastq.Reader(in_fastq_fname)
    for r_id, exp_id, r_randomer, r_seq, r_plus, r_qual in \
            extract(reader, barcodes, mismatches=mismatches,
                    minimum_length=minimum_length):
        if r_randomer:
            r_id = "%s:%s" % (r_id, r_randomer)

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
        assert(i2start_pos.setdefault(bi, len(bar5)) == len(bar5))

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

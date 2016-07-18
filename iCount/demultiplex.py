"""
Demultiplexing
--------------

Remove sample barcode, random barcode, and adapter sequences.

"""



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
        'outdir', 'string', '.', False,
        'Output folder. Use local folder if none given.'
    ),
    (
        'prefix', 'string', 'exp', False,
        'Prefix of generated FASTQ files.'
    ),
    (
        'adapter', 'string', 'AGATCGGAAGAGCGGTTCAG', False,
        'Adapter sequence to remove from ends of reads.'
    ),
    (
        'mismatches', 'int_range', (1, 0, 5), False,
        'Number of tolerated mismatches when comparing barcodes.'
    )
]

params_pos = [
    (
        'sequences', 'FASTQ', 'in',
        '(input) sequences from a sequencing library.'
    ),
]

def demultiplex(in_fastq_fname, out_fastq_fnames, not_matching_fastq_fname,
                barcodes):
    """Extract reads and save to individual FASTQ files.

    All non-matching reads are stored in FASTQ file not_matching_fastq_fname.

    :param in_fastq_fname: filename to read,  must be FASTQ in order to
    avoid duplicate read records.
    :param out_fastq_fnames: list of filenames to store reads as determined
    by barcodes
    : not_matching_fastq_fname: fastq filename where to store reads not
    matching
    :param barcodes: experiment and randomer barcode definition
    """

    out_fastq = [iCount.files.fastq.Writer(fn) for fn in
                 out_fastq_fnames + [not_matching_fastq_fname]]
    reader = iCount.files.fastq.Reader(in_fastq_fname)
    for r_id, exp_id, r_randomer, r_seq, r_plus, r_qual in extract(reader,
                                                                 barcodes):
        if r_randomer:
            r_id = "%s:%s" % (r_id, r_randomer)

        out_fastq[exp_id].write(r_id, r_seq, '+', r_qual)
    while out_fastq:
        f = out_fastq.pop()
        f.close()
    return out_fastq_fnames

def extract(seqs, barcodes):
    """Return iterator that returns experiment, randomer, and remaining
    sequence.

    :param seqs: input sequence of tuples (r_id, r_seq, r_qualL)
    :param barcodes: dictionary of experiment annotation
    :return: iterator
    """
    valid_nucs = set(['A', 'T', 'C', 'G'])

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
            for n in ns:
                p2n2i.setdefault(p, {}).setdefault(n, set()).add(bi)
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
        #r_id, r_seq, r_plus, r_qualL
        votes = [0]*all_barcodes
        for p, n2i in p2i:
            bis = n2i.get(r.r_seq[p], [])
            for bi in bis:
                votes[bi] = votes[bi] + 1
        mvotes = max(votes)
        if mvotes < max_votes - 1: # allow one mismatch
            # not recognized
            exp_id = -1

            r_randomer = ''
            r_seq = r.r_seq
            r_qual = r.r_qual
        else:
            exp_id = votes.index(mvotes)

            r_randomer = ''.join(r.r_seq[p] for p in random_pos)
            ps = i2start_pos[exp_id]
            r_seq = r.r_seq[ps:]
            r_qual = r.r_qual[ps:]

        yield r.r_id, exp_id, r_randomer, r_seq, r.r_plus, r_qual


def remove_adapter(in_fastq, out_fastq, adapter):
    """

    :param adapter: sequence to be found at 3'
    :return:
    """

    iCount.externals.cutadapt.run(in_fastq, out_fastq, adapter, qual_trim=10)
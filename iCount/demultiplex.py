""".. Line to protect from pydocstyle D205, D400.

Demultiplexing
==============

Split FASTQ file into separate files, one for each sample barcode.

Saved FASTQ files contain sequences where sample barcode, random
barcode, and adapter sequences are removed. Random barcode is moved into
the header line, since it is needed in later steps (removing PCR duplicates
and counting number of cross-link events).

.. autofunction:: iCount.demultiplex.run

"""

import os
import logging
import shutil

import iCount

from iCount.externals.cutadapt import run as remove_adapter


LOGGER = logging.getLogger(__name__)

VALID_NUCS = {'A', 'T', 'C', 'G'}


def run(reads, adapter, barcodes, mismatches=1, minimum_length=15, prefix='demux', out_dir='.'):
    """
    Demultiplex FASTQ file.

    Split input FASTQ file into separate files, one for each barcode, and
    additional file for non-matching barcodes.

    Parameters
    ----------
    reads : str
        Path to reads from a sequencing library.
    adapter : str
        Adapter sequence to remove from ends of reads.
    barcodes : list_str
        List of barcodes used for library.
    mismatches : int
        Number of tolerated mismatches when comparing barcodes.
    minimum_length : int
        Minimum length of trimmed sequence to keep.
    prefix : str
        Prefix of generated FASTQ files.
    out_dir : str
        Output folder. Use local folder if none given.

    Returns
    -------
    str
        List of filenames where separated reads are stored.

    """
    iCount.log_inputs(LOGGER, level=logging.INFO)
    metrics = iCount.Metrics()
    metrics.reads_ok = 0
    metrics.reads_fail = 0

    if not os.path.isdir(out_dir):
        raise FileNotFoundError('Output directory does not exist. Make sure it does.')

    out_fnames = [os.path.abspath(os.path.join(out_dir, '{}_{}_tmp.fastq.gz'.format(
        prefix, barcode))) for barcode in barcodes + ['nomatch']]

    LOGGER.info('Demultiplexing...')
    # Make list of file handles - one for each output filename:
    out_fastqs = [iCount.files.fastq.FastqFile(fname, 'wt') for fname in out_fnames]

    # Determine experiment ID and random barcode for each fastq entry:
    kwargs = {'mismatches': mismatches, 'minimum_length': minimum_length}
    for fq_entry, exp_id, randomer in _extract(reads, barcodes, **kwargs):
        if randomer:
            metrics.reads_ok += 1
            if fq_entry.id[-2] == '/':
                # For early versions of Illumina, keep mate info at the end:
                r_pair = fq_entry.id[-2:]
                fq_entry.id = '{}:rbc:{}{}'.format(fq_entry.id[:-2], randomer, r_pair)
            else:
                fq_entry.id = '{}:rbc:{}'.format(fq_entry.id, randomer)
        else:
            metrics.reads_fail += 1
        out_fastqs[exp_id].write(fq_entry)

    for out_fastq in out_fastqs:
        out_fastq.close()

    # Finally, remove adapters (if requested) or just rename to final names:
    out_fnames_final = ['{}.fastq.gz'.format(fname[:-13]) for fname in out_fnames]
    for fname_in, fname_out in zip(out_fnames, out_fnames_final):
        if adapter and 'nomatch' not in fname_in:
            remove_adapter(
                fname_in, adapter, reads_trimmed=fname_out, minimum_length=minimum_length)
            os.remove(fname_in)
        else:
            shutil.move(fname_in, fname_out)

    return metrics


def _make_p2n2i(barcodes):
    """
    Create a mapping p2n2i (Position-to-Nucleotide-to-experiment_Id).

    * For each position tell which nucleotides are possible,
    * For each of possible nucletides, tell which experiments are possible
    * However, do not include positions where all nucleotides are N's
    """
    p2n2i = {}
    for exp_id, barcode in enumerate(barcodes):
        for pos, nuc in enumerate(barcode):
            if nuc != 'N':
                assert nuc in VALID_NUCS
                p2n2i.setdefault(pos, {}).setdefault(nuc, set()).add(exp_id)
    return p2n2i


def _extract(reads, barcodes, mismatches=1, minimum_length=15):
    """
    Get experiment ID, randomer and remaining sequence for each read in FASTQ file.

    Parameters
    ----------
    reads : str
        Fastq file name.
    barcodes : list_str
        Barcodes: N's are random, others are experiment positions.
    mismatches : int
        Number of allowed mismatches in experiment barcode.
    minimum_length : int
        Discard reads with length less than this.

    Yields
    ------
    tuple
        FastqEntry, experiment ID, randomer

    """
    p2n2i = _make_p2n2i(barcodes)
    all_barcodes = len(barcodes)
    # Precompute barcode length for each barcode
    barcode_len = {i: len(brc) for i, brc in enumerate(barcodes)}
    longest_barcode = max(barcode_len.values())
    # Precompute max possible votes for each barcode. Equals to the number of valid nucleotides
    max_votes = {i: len(brc.replace('N', '')) for i, brc in enumerate(barcodes)}
    # Precompute positions of randomer nucleotides for each barcode
    randomer_pos = {
        i: [j for j, nuc in enumerate(brc) if nuc == 'N'] for i, brc in enumerate(barcodes)}

    # Determine experiment ID and extract random barcode for each FASTQ entry.
    # Experiment ID si determined by "voting": at each non-randomer position
    # `pos` there is nucleotide X. Check which experiments (=barcodes) have X
    # at position `pos`. Increment votes for them.
    for read in iCount.files.fastq.FastqFile(reads).read():

        if len(read.seq) < longest_barcode:
            continue

        votes = [0] * all_barcodes
        for pos, n2i in p2n2i.items():
            for exp_id in n2i.get(read.seq[pos], []):
                votes[exp_id] += 1

        # Count mismatches for each barcode. Experiment with least mismatches is the winner.
        xmatches = [max_votes[i] - votes[i] for i, barcode in enumerate(barcodes)]
        min_mismatches = min(xmatches)
        if min_mismatches > mismatches:  # not recognized as valid barcode
            experiment_id = -1
            randomer = ''
            seq = read.seq
            qual = read.qual
        else:  # recognized as valid barcode
            experiment_id = xmatches.index(min_mismatches)
            randomer = ''.join(read.seq[p] for p in randomer_pos[experiment_id])
            # Keep only the remainder of the sequence / quality scores
            seq = read.seq[barcode_len[experiment_id]:]
            qual = read.qual[barcode_len[experiment_id]:]

        if len(seq) >= minimum_length:
            new_fq_entry = iCount.files.fastq.FastqEntry(read.id, seq, '+', qual)
            yield (new_fq_entry, experiment_id, randomer)

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
            if fq_entry.id[-2] == '/':
                # For early versions of Illumina, keep mate info at the end:
                r_pair = fq_entry.id[-2:]
                fq_entry.id = '{}:rbc:{}{}'.format(fq_entry.id[:-2], randomer, r_pair)
            else:
                fq_entry.id = '{}:rbc:{}'.format(fq_entry.id, randomer)
        out_fastqs[exp_id].write(fq_entry)

    for out_fastq in out_fastqs:
        out_fastq.close()

    # Finally, remove adapters (if requested) or just rename to final names:
    out_fnames_final = ['{}.fastq.gz'.format(fname[:-13]) for fname in out_fnames]
    for fname_in, fname_out in zip(out_fnames, out_fnames_final):
        if adapter and 'nomatch' not in fname_in:
            remove_adapter(fname_in, fname_out, adapter, minimum_length=minimum_length)
            os.remove(fname_in)
        else:
            shutil.move(fname_in, fname_out)

    return out_fnames_final


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
    all_barcodes = len(barcodes)
    barcode_len = len(barcodes[0])
    p2n2i = _make_p2n2i(barcodes)
    max_votes = len(p2n2i)
    randomer_pos = [i for i in range(barcode_len) if i not in p2n2i]

    # Determine experiment ID and extract random barcode for each FASTQ entry.
    # Experiment ID si determined by "voting": at each non-randomer position
    # `pos` there is nucleotide X. Check which experiments (=barcodes) have X
    # at position `pos`. Increment votes for them. Experiment with max votes is
    # the winner, but only if votes equal/exceed max_votes - mismatches.
    for read in iCount.files.fastq.FastqFile(reads).read():
        votes = [0] * all_barcodes
        for pos, n2i in p2n2i.items():
            for exp_id in n2i.get(read.seq[pos], []):
                votes[exp_id] += 1
        mvotes = max(votes)
        if mvotes < max_votes - mismatches:  # not recognized as valid barcode
            experiment_id = -1
            randomer = ''
            seq = read.seq
            qual = read.qual
        else:  # recognized as valid barcode
            experiment_id = votes.index(mvotes)
            randomer = ''.join(read.seq[p] for p in randomer_pos)
            # Keep only the remainder of the sequence / quality scores
            seq = read.seq[barcode_len:]
            qual = read.qual[barcode_len:]

        if len(seq) >= minimum_length:
            new_fq_entry = iCount.files.fastq.FastqEntry(read.id, seq, '+', qual)
            yield (new_fq_entry, experiment_id, randomer)

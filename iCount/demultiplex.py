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

import re
import os
import logging
from collections import OrderedDict

import iCount

from iCount.externals.cutadapt import run as remove_adapter


LOGGER = logging.getLogger(__name__)

VALID_NUCS = {'A', 'T', 'C', 'G'}


def get_relative5_positions(brc5):
    """Get relative positions of barcode nucleotides from 5'end."""
    return list(range(0, len(brc5)))


def get_relative3_positions(brc3):
    """Get relative positions of barcode nucleotides from 3'end."""
    return list(range(-len(brc3), 0))


def prepare_barcodes(barcodes5, barcodes3):
    """Prepare barcodes for further processing."""
    if not barcodes3:
        barcodes3 = ['.'] * len(barcodes5)
        if len(barcodes5) != len(set(barcodes5)):
            raise ValueError("Barcodes should not repeat.")
    else:
        combined = list(zip(barcodes5, barcodes3))
        if len(combined) != len(set(combined)):
            raise ValueError("Comnbination of 3' and 5' barcodes should be unique.")

    if len(barcodes3) != len(barcodes5):
        raise ValueError('Parameter barcodes3 should contain same amount of barcodes as barcodes5 parameter.')

    barcodes = OrderedDict()
    # In first pass fill data for barcodes5
    for brc5 in barcodes5:
        rnd_pos5 = [pos for pos in get_relative5_positions(brc5) if brc5[pos] == 'N']
        barcodes[brc5] = {
            'rel_pos': get_relative5_positions(brc5),
            'rnd_pos': rnd_pos5,
            'rnd_len': len(rnd_pos5),
            'len': len(brc5),
            'barcodes3': OrderedDict(),
        }

    # In second pass fill data for barcodes3
    for brc5, brc3 in zip(barcodes5, barcodes3):
        if brc3 == '.' or brc3 == '':
            continue

        rel_pos3 = get_relative3_positions(brc3)
        rnd_pos3 = [pos for pos in rel_pos3 if brc3[pos] == 'N']
        brc3_data = {
            'rel_pos': get_relative3_positions(brc3),
            'rnd_pos': rnd_pos3,
            'rnd_len': len(rnd_pos3),
            'len': len(brc3),
        }
        barcodes[brc5]['barcodes3'][brc3] = brc3_data

    return barcodes


def create_index(barcodes):
    """
    Create a mapping: Position-to-Nucleotide-to-barcode.

    * For each position tell which nucleotides are possible,
    * For each of possible nucletides, tell which barcodes are possible
    * However, do not include positions where all nucleotides are N's
    """
    index = {}
    for brc, brc_data in barcodes.items():
        for pos in brc_data['rel_pos']:
            nuc = brc[pos]
            if nuc != 'N':
                assert nuc in VALID_NUCS
                index.setdefault(pos, {}).setdefault(nuc, set()).add(brc)
    return index


def _extract(reads, barcodes, **kwargs):
    """
    Get experiment barcode, randomer and remaining sequence for each read in FASTQ file.

    Parameters
    ----------
    reads : str
        Fastq file name.
    barcodes : dict
        Barcodes.
    kwargs : dict
        Kwargs with analysis parameters.

    Yields
    ------
    tuple
        FastqEntry, experiment ID, randomer

    """
    index = create_index(barcodes)

    # Precompute longest barcode.
    longest_barcode = max([brc_data['len'] for brc_data in barcodes.values()])

    # Precompute max possible votes for each barcode. Equals to the number of valid nucleotides.
    max_votes = {brc: brc_data['len'] - brc_data['rnd_len'] for brc, brc_data in barcodes.items()}

    # Precompute positions of randomer nucleotides for each barcode.
    randomer_pos = {brc: brc_data['rnd_pos'] for brc, brc_data in barcodes.items()}

    # Determine experiment barcode and extract random barcode for each FASTQ entry.
    # Experiment barcode si determined by "voting": at each non-randomer position
    # `pos` there is nucleotide X. Check which experiments/barcodes have X
    # at position `pos`. Increment votes for them.
    for read in iCount.files.fastq.FastqFile(reads).read():

        if len(read.seq) < longest_barcode:
            continue

        votes = {brc: 0 for brc in barcodes}
        for pos, nuc_to_brc in index.items():
            for brc in nuc_to_brc.get(read.seq[pos], []):
                votes[brc] = votes.get(brc, 0) + 1

        # Count mismatches for each barcode. Experiment with least mismatches is the winner.
        xmatches = {brc: max_votes[brc] - votes[brc] for brc in barcodes}
        min_mismatches = min(xmatches.values())
        if min_mismatches > kwargs['mismatches']:  # not recognized as valid barcode
            winner_barcode = 'nomatch'
            randomer = ''
            seq = read.seq
            qual = read.qual
        else:  # recognized as valid barcode
            winner_barcode = next(brc for brc, xmatch in xmatches.items() if xmatch == min_mismatches)
            randomer = ''.join(read.seq[p] for p in randomer_pos[winner_barcode])
            # Keep only the remainder of the sequence / quality scores
            if barcodes[winner_barcode]['rel_pos'][0] < 0:
                # This are 3' end barcodes
                brc_len = barcodes[winner_barcode]['len']
                seq = read.seq[:-brc_len]
                qual = read.qual[:-brc_len]
            else:
                # This are 5' end barcodes
                brc_len = barcodes[winner_barcode]['len']
                seq = read.seq[brc_len:]
                qual = read.qual[brc_len:]

        if len(seq) >= kwargs['minimum_length']:
            new_fq_entry = iCount.files.fastq.FastqEntry(read.id, seq, '+', qual)
            yield (new_fq_entry, winner_barcode, randomer)


def add_randomer_to_header(randomer, fq_entry):
    """Add randomer info to FASTQ header."""
    match = re.match(r'.*(:rbc:)([ACGTN]+).*', fq_entry.id)
    if match:
        # Handle the case where randomer is already in the header.
        rbc = match.group(2)
        randomer = randomer + rbc

        fq_entry.id = re.sub(r':rbc:[ACGTN]+', '', fq_entry.id)

    if fq_entry.id[-2:] in ['/1', '/2']:
        # For early versions of Illumina, keep mate info at the end:
        fq_entry.id = '{}:rbc:{}{}'.format(fq_entry.id[:-2], randomer, fq_entry.id[-2:])
    else:
        fq_entry.id = '{}:rbc:{}'.format(fq_entry.id, randomer)


def demultiplex(reads, barcodes, **kwargs):
    """Demultiplex."""
    out_fastqs = {}
    for brc in list(barcodes.keys()) + ['nomatch']:
        basename = '{}_{}.fastq.gz'.format(kwargs['prefix'], brc)
        fname = os.path.abspath(os.path.join(kwargs['out_dir'], basename))
        out_fastqs[brc] = iCount.files.fastq.FastqFile(fname, 'wt')

    # Determine experiment ID and random barcode for each fastq entry:
    for fq_entry, wbrc, randomer in _extract(reads, barcodes, **kwargs):
        if randomer:
            kwargs['metrics'].reads_ok += 1
            add_randomer_to_header(randomer, fq_entry)
        else:
            kwargs['metrics'].reads_fail += 1
        out_fastqs[wbrc].write(fq_entry)

    for out_fastq in out_fastqs.values():
        out_fastq.close()


def run(reads, adapter, barcodes5, barcodes3=None, mismatches=1, minimum_length=15, min_adapter_overlap=7,
        prefix='demux', out_dir='.'):
    """
    Demultiplex FASTQ file.

    Split input FASTQ file into separate files, one for each barcode, and
    additional file for non-matching barcodes. Write random barcode of a read
    into it's FASTQ header row.

    Parameters
    ----------
    reads : str
        Sequencing reads.
    adapter : str
        Adapter sequence to remove from 3-prime end of reads.
    barcodes5 : list_str
        List of 5-prime end barcodes.
    barcodes3 : list_str
        List of 3-prime end barcodes.
    mismatches : int
        Number of tolerated mismatches when comparing barcodes.
    minimum_length : int
        Minimum length of trimmed sequence to keep.
    min_adapter_overlap : int
        Minimum length of adapter on 3' end if demultiplexing also on
        3' barcodes.
    prefix : str
        Prefix of generated FASTQ files.
    out_dir : str
        Output folder. Use current folder if none is given.

    Returns
    -------
    iCount.Metrics
        Metrics object, storing analysis metadata.

    """
    iCount.log_inputs(LOGGER, level=logging.INFO)
    metrics = iCount.Metrics()
    metrics.reads_ok = 0
    metrics.reads_fail = 0

    kwargs = {
        'mismatches': mismatches,
        'minimum_length': minimum_length,
        'prefix': prefix,
        'out_dir': out_dir,
        'metrics': metrics,
    }

    if not os.path.isdir(out_dir):
        raise FileNotFoundError('Output directory does not exist. Make sure it does.')

    # This should be a dict, where each 5' barcode has al list of it's 3' barcods
    barcodes = prepare_barcodes(barcodes5, barcodes3)

    LOGGER.info("Demultiplexing based on 5' barcodes...")
    # Demultiplex by 5' barcodes only. Store in files. Just like before.
    demultiplex(reads=reads, barcodes=barcodes, **kwargs)

    os.rename(
        os.path.join(out_dir, 'demux_nomatch.fastq.gz'),
        os.path.join(out_dir, 'demux_nomatch5.fastq.gz'),
    )

    LOGGER.info("Demultiplexing based on 3' barcodes...")
    for barcode5 in barcodes:
        reads5 = os.path.join(out_dir, 'demux_{}.fastq.gz'.format(barcode5))
        barcodes3 = barcodes[barcode5]['barcodes3']

        if not barcodes3:
            # This barcode has no 3' counterparts. Just remove the adapter and continue
            # TODO: polish the parameters for adapter removal in this case...
            remove_adapter(reads5, adapter, overwrite=True)
            continue

        # One must be sure that there actually are 3' barcodes on the
        # 3' end. In cca. 20-30% of cases read is so short that it does
        # not reach the 3' barcode. In such cases, demultiplexing by 3'
        # end would be done by random chance which is not acceptable not
        # good. To be sure that 3' barcode is reached, read needs to
        # contain at least ``adapter_overlap`` bp of the adapter.

        no_adapters = os.path.join(out_dir, "no_adapter_found_{}.fastq.gz".format(barcode5))
        remove_adapter(reads5, adapter, overwrite=True, overlap=min_adapter_overlap, untrimmed_output=no_adapters)

        # Fix the prefix, to include 5' barcode info:
        kwargs['prefix'] = '{}_{}'.format(prefix, barcode5)
        kwargs['mismatches'] = 0

        # Now, demutiplex based on 3' adapter
        demultiplex(reads=reads5, barcodes=barcodes3, **kwargs)

        # File that is demultiplexed only by 5'end is not needed anymore
        os.remove(reads5)

    # TODO: merge nomatch files 3' and 5' end... This may be many!

    return metrics

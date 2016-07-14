"""Pipeline for processing iCLIP data.



"""
import iCount


def bam2beds(f_in, f_out, max_multi_count):
    # single, unique hits to genome
    # multiple hits, up to 20 places, spread to all positions hit
    pass


def process_lib(in_fastq_fname, barcodes, adapter, map_to,
                max_multi_hits=20, out_folder='.'):
    """Take all steps to generate bedGraphs for each experiment


    """

    out_files = {}

    # demultiplex, make sure to convert to Illumina 1.8+ Phred+33
    demux_fastqs = []
    for b5, b3 in barcodes:
        demux_fastqs.append("%s/demux_%s.fastq" % (out_folder, b5))

    out_nomatch_fastq = "%s/demux_nomatch.fastq" % (out_folder)
    iCount.reads.demultiplex(in_fastq_fname, demux_fastqs, out_nomatch_fastq,
                             barcodes)

    # remove adapter sequences and trim low quality
    trimmed_fastqs = []
    for (b5, b3), demux_fastq in zip(barcodes, demux_fastqs):
        trimmed_fastq = "%s/trimmed_%s.fastq" % (out_folder, b5)
        trimmed_fastqs.append(trimmed_fastq)
        iCount.reads.remove_adapter(demux_fastq, trimmed_fastq, adapter)

    return out_files

    # map
    to_filter = {}
    for mt, fin_name in zip(map_to, demux_fastqs):
        fout_name = u"{0:s}_m.bam".format(fin_name[:-6])
        iCount.externals.tophat(fin_name, fout_name, mt, max_multi_hits)
        to_filter.setdefault(mt, []).append(fout_name)
        out_files.append(fout_name)

    # filter
    for mt, fin_names in to_filter.items():
        fout_names = [u"{0:s}_f.bam".format(fn[:-4]) for fn in fin_names]
        iCount.analysis.hits.filter(fin_names, fout_names)
        out_files.extend(fout_names)

        # bedgraphs
        for fn in fout_names:
            # unique, single hits
            fout_u = u"{0:s}_U.bed".format(fn[:-4])
            bam2beds(fn, fout_u, max_multi_count=1)
            out_files.append(fout_u)

            # all, multiple hits
            fout_a = u"{0:s}_A.bed".format(fn[:-4])
            bam2beds(fn, fout_a, max_multi_count=max_multi_hits)
            out_files.append(fout_a)

    return out_files

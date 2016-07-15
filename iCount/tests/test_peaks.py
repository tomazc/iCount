import unittest

import os
import shutil
from . import test_data

import pybedtools
import iCount

files_to_get = [
    (
        'hnRNPC_5_S.bed.gz',

        'http://icount.fri.uni-lj.si/results/ensembl59/20101116_LUjh03'
        '/bedGraph_cDNA'
        '/20101116_LUjh03_5_hg19_ensembl59_S_iCLIP_hnRNPC_HeLa_bedGraph-cDNA'
        '.bed.gz'
    ),
    (
        'hnRNPC_5_S_peaks_scores.tab',

        'http://icount.fri.uni-lj.si/results/analyses/ensembl59/peaks/51764'
        '/peaks_id51764_rnd100_flank15_fdr0'
        '.05_20101116_LUjh03_5_hg19_ensembl59_S_iCLIP_hnRNPC_HeLa_bedGraph'
        '-cDNA.bed.gz_scores.tab.gz'
    ),
    (
        'hnRNPC_5_S_peask.bed.gz',

        'http://icount.fri.uni-lj.si/results/analyses/ensembl59/peaks/51764'
        '/peaks_id51764_rnd100_flank15_fdr0'
        '.05_20101116_LUjh03_5_hg19_ensembl59_S_iCLIP_hnRNPC_HeLa_bedGraph'
        '-cDNA.bed.gz_lowFDR.bed.gz',
    ),
    (
        'hg19.gtf.gz',

        'ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode'
        '.v19.annotation.gtf.gz'
    ),
    (
        'hg19.gff3.gz',

        'ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode'
        '.v19.annotation.gff3.gz'
    )
]


def setUpModule():
    test_data.get_update(files_to_get)


def tearDownModule():
    pass


test_output_folder = os.path.join(test_data.test_folder, 'runs', 'peaks')
if not os.path.exists(test_output_folder):
    print("Test folder does not exist. Will create it at: "
          "{0:s}".format(test_output_folder))
    os.makedirs(test_output_folder)


class Testpeaks(unittest.TestCase):

    def setUp(self):
        if os.path.isdir(test_output_folder):
            print('Removing folder: {:s}'.format(test_output_folder))
            shutil.rmtree(test_output_folder)
            print('Creating folder: {:s}'.format(test_output_folder))
            os.makedirs(test_output_folder)

    def test_peaks(self):
        return
        fin_annotation = os.path.join(test_data.test_data_folder,
                                      'hg19.gtf.gz')
        fin_bedGraph = os.path.join(test_data.test_data_folder,
                                    'hnRNPC_5_S.bed.gz')
        fout_bedGraph = os.path.join(test_data.test_data_folder,
                                     'peaks_out.bed.gz')
        fout_scores = os.path.join(test_data.test_data_folder,
                                   'peaks_out_scores.tab.gz')

        fin_converted_bedGraph = os.path.join(test_output_folder,
                                    'hnRNPC_5_S.conv.bed.gz')

        def convert_legacy_bed_format(feature):
            # use BED6 format, see:
            # http://bedtools.readthedocs.io/en/latest/content/general-usage.html
            chrom = feature.chrom
            start = feature.start
            end = feature.stop
            name = '.'
            if feature.name[0] == '-' or feature.name[0] == '+':
                score = feature.name[1:]
                strand = feature.name[0]
            else:
                score = feature.name
                strand = '+'
            return pybedtools.create_interval_from_list(
                [chrom, start, end, name, score, strand]
            )

        sites = pybedtools.BedTool(fin_bedGraph).sort().saveas()
        sites.each(convert_legacy_bed_format).saveas(fin_converted_bedGraph)

        iCount.analysis.peaks.run(fin_annotation, fin_converted_bedGraph,
                                  fout_bedGraph, fout_scores=fout_scores)


    def test_sum_within_window(self):
        self.assertEqual(iCount.analysis.peaks.sum_within_window([]), [])
        sites = [
            (10, 1), (11, 1), (12, 1), (13, 1), (14, 1), (15, 1), (16, 1),
            (20, 2)
        ]
        summed_sites_w1 = [
            (10, 2), (11, 3), (12, 3), (13, 3), (14, 3), (15, 3), (16, 2),
            (20, 2)
        ]
        summed_sites_w3 = [
            (10, 4), (11, 5), (12, 6), (13, 7), (14, 6), (15, 5), (16, 4),
            (20, 2)
        ]
        self.assertEqual(
            iCount.analysis.peaks.sum_within_window(sites, w=1),
            summed_sites_w1
        )
        self.assertEqual(
            iCount.analysis.peaks.sum_within_window(sites, w=3),
            summed_sites_w3
        )

if __name__ == '__main__':
    unittest.main()

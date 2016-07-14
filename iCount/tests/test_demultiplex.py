import unittest

import os
import shutil
from . import test_data
import iCount

files_to_get = [
    (
        'LUjh03.fq.gz',
        'http://icount.fri.uni-lj.si/data/20101116_LUjh03/SLX-2605'
        '.CRIRUN_501.s_4.sequence.txt.gz'
    ),
]


def setUpModule():
    test_data.get_update(files_to_get)


def tearDownModule():
    pass


test_output_folder = os.path.join(test_data.test_folder, 'runs', 'iCLIP')
if not os.path.exists(test_output_folder):
    print("Test folder does not exist. Will create it at: "
          "{0:s}".format(test_output_folder))
    os.makedirs(test_output_folder)


class TestDemultiplex(unittest.TestCase):

    def setUp(self):
        if os.path.isdir(test_output_folder):
            print('Removing folder: {:s}'.format(test_output_folder))
            shutil.rmtree(test_output_folder)
            print('Creating folder: {:s}'.format(test_output_folder))
            os.makedirs(test_output_folder)

    def test_demultiplex(self):
        in_fastq_fname = os.path.join(test_data.test_data_folder,
                                      'LUjh03_reduced.fq.gz')

        # extract
        # list all experiment barcodes
        barcodes = [
            'NNNGGTTNN',
            'NNNTTGTNN',
            'NNNCAATNN',
            'NNNACCTNN',
            'NNNGGCGNN',
        ]
#        adapter = 'AGATCGGAAGAGCGGTTCAG'  # 'AGATCGGAAG_1,AGCGGTTCAG_2'

        out_fastq_fnames = []
        for bc in barcodes:
            out_fastq_fnames.append(
                os.path.join(test_output_folder, 'demulti.{}.fastq'.format(bc))
            )
        not_matching_fastq_fname = os.path.join(test_output_folder,
                                                'nomatch.fastq')

        iCount.demultiplex.demultiplex(in_fastq_fname, out_fastq_fnames,
                                       not_matching_fastq_fname, barcodes)



if __name__ == '__main__':
    unittest.main()

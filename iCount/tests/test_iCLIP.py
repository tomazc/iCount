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


class TestiCLIP(unittest.TestCase):

    def setUp(self):
        if os.path.isdir(test_output_folder):
            print('Removing folder: {:s}'.format(test_output_folder))
            shutil.rmtree(test_output_folder)
            print('Creating folder: {:s}'.format(test_output_folder))
            os.makedirs(test_output_folder)

    def test_iCLIP(self):
        in_fastq_fname = os.path.join(test_data.test_data_folder,
                                      'LUjh03_reduced.fq.gz')

        # extract
        # list all experiment barcodes
        barcodes = [
            ('NNNGGTTNN', ''),
            ('NNNTTGTNN', ''),
            ('NNNCAATNN', ''),
            ('NNNACCTNN', ''),
            ('NNNGGCGNN', ''),
        ]
        adapter = 'AGATCGGAAGAGCGGTTCAG'  # 'AGATCGGAAG_1,AGCGGTTCAG_2'

        map_to = [
            'hg19',
            'hg19',
            'hg19',
            'hg19',
            'hg19',
        ]

        iCount.analysis.iCLIP.process_lib(in_fastq_fname, barcodes, adapter,
                                          map_to,
                                          out_folder=test_output_folder)


if __name__ == '__main__':
    unittest.main()

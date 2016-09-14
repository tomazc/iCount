import os
import shutil
import unittest

from . import test_data
import iCount

files_to_get = [
    (
        'LUjh03.fq.gz',
        'http://icount.fri.uni-lj.si/data/20101116_LUjh03/SLX-2605'
        '.CRIRUN_501.s_4.sequence.reduced.txt.gz'
    ),
]


def setUpModule():
    test_data.get_update(files_to_get)


def tearDownModule():
    pass


test_output_folder = os.path.join(test_data.test_folder, 'runs', 'cutadapt')
if not os.path.exists(test_output_folder):
    print("Test folder does not exist. Will create it at: "
          "{0:s}".format(test_output_folder))
    os.makedirs(test_output_folder)


class TestCutadapt(unittest.TestCase):

    def setUp(self):
        if os.path.isdir(test_output_folder):
            print('Removing folder: {:s}'.format(test_output_folder))
            shutil.rmtree(test_output_folder)
            print('Creating folder: {:s}'.format(test_output_folder))
            os.makedirs(test_output_folder)

    def test_cutadapt(self):
        in_fastq_fname = os.path.join(
            test_data.test_data_folder,
            'LUjh03.fq.gz',
        )
        out_fastq_fname = os.path.join(
            test_output_folder,
            'out.fastq',
        )
        adapter = "AGATCGGAAGAGCGGTTCAG"

        rc = iCount.externals.cutadapt.run(in_fastq_fname, out_fastq_fname,
                                           adapter)
        self.assertEqual(rc, 0)

    def test_version(self):
        ver = iCount.externals.cutadapt.get_version()
        self.assertEqual(ver, iCount.externals.expected_cutadapt_version)


if __name__ == '__main__':
    unittest.main()

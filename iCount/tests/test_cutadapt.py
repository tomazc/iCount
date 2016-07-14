import unittest

import os
import iCount

class TestCutadapt(unittest.TestCase):
    def test_cutadapt(self):
        in_fastq_fname = os.path.join(
            iCount.tmp_root,
#           '20101116_LUjh03/SLX-2605.CRIRUN_501.s_4.sequence.txt.gz'
            'in_small.fastq'
        )
        out_fastq_fname = os.path.join(
            iCount.tmp_root,
            'out.fastq'
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

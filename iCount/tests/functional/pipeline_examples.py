"""
Test consistency of example scripts.

This script tests if results of iCount/examples/*.sh are consistent with the
reference results stored under
iCount/tests/pipeline_examples.[script_name].out.txt

If reference results are not present, they are created and saved under
iCount/tests/pipeline_examples.[script_name].out-cur.txt
"""
# pylint: disable=missing-docstring, protected-access

import os
import shutil
import unittest

from . import FUNCTIONAL_TEST_FOLDER


def _read_file(fn):
    with open(fn) as file_:
        # skip STAR status reports because time-stamped
        lines = [r.strip() for r in file_]
        return [r for r in lines if ' ..... ' not in r and ' ... ' not in r]


class TestPipeline(unittest.TestCase):

    def setUp(self):
        # prefix of files with reference results
        self.ref_pref = os.path.abspath(__file__[:-3])

        # folder where all examples will be run
        self.test_dir = os.path.join(FUNCTIONAL_TEST_FOLDER, 'pipeline')
        # remove it first
        if os.path.exists(self.test_dir):
            shutil.rmtree(self.test_dir)
        # create an empty folder
        if not os.path.exists(self.test_dir):
            os.makedirs(self.test_dir)
        cmd = 'cd "{:s}"; iCount examples'.format(self.test_dir)
        print('Executing: {:s}'.format(cmd))
        os.system(cmd)
        self.examples_dir = os.path.join(self.test_dir, 'examples')
        self.assertTrue(os.path.isdir(self.examples_dir))

        self.examples = ['hnRNPC_reduced.sh']
        self.files_to_remove = []

    def tearDown(self):
        for fn in self.files_to_remove:
            os.remove(fn)

    def test_examples(self):
        cur_results = {}
        for script in self.examples:
            cur_fn = '{:s}.{:s}.out-cur.txt'.format(self.ref_pref, script)
            cmd = 'cd "{:s}"; ./{:s} > {:s} 2>&1'.format(self.examples_dir,
                                                         script, cur_fn)
            print('Executing: {:s}'.format(cmd))
            os.system(cmd)
            # confirm that current results file was created
            self.assertTrue(os.path.isfile(cur_fn))
            cur_results[script] = cur_fn

        cur_lines = {}
        ref_lines = {}
        for script in self.examples:
            # confirm that reference file exists
            ref_fn = '{:s}.{:s}.out.txt'.format(self.ref_pref, script)
            self.assertTrue(os.path.isfile(ref_fn))

            cur_fn = cur_results[script]

            # remove temporary and time stamped entries
            ref = _read_file(ref_fn)
            ref = [r for r in ref if '/tmp/iCount/star' not in r]
            ref = [r for r in ref if 'URL:http://icount' not in r]
            ref_lines[script] = ref

            cur = _read_file(cur_fn)
            cur = [r for r in cur if '/tmp/iCount/star' not in r]
            cur = [r for r in cur if 'URL:http://icount' not in r]
            cur_lines[script] = cur

            if ref_lines[script] == cur_lines[script]:
                self.files_to_remove.append(cur_fn)

        for script in self.examples:
            self.assertIn(script, cur_lines)
            self.assertIn(script, ref_lines)

            self.assertEqual(ref_lines[script], cur_lines[script])


if __name__ == '__main__':
    unittest.main()

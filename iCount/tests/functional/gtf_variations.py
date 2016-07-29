"""
Tests how *.gtf file format changes during different releases.

The whole script can be run as a test and compare if reference
results agree with the current ones. If reference results are not
present, they are created.

If gtf file for given release and species is not found in `self.gtf_files_dir`
it is downloaded from ftp.ensembl.org by using function
iCount.genomes.ensembl.download_annotation. This may significantly extend
the execution time of this script.

By modifying the variables `self.releases_list` and `self.species_list` user can
determine which gtf files to put under test.
"""

import os
import re
import time
import sys
import pickle
import unittest

import pybedtools

from iCount.genomes import ensembl
from __init__ import functional_test_folder

SEPARATOR = '-' * 72


def produce_report(infile, outfile, report_columns=[1, 2, 8]):
    """
    Reports how column values in GTF file format change between releases.

    GTF file consits from lines with 9 columns. Function enables to examine
    desired set of colums, but the default is to invetigate the 2nd, 3rd
    and 9th. Detailed description and examples of of GTF file format can
    be found here:
    http://mblab.wustl.edu/GTF22.html

    :param str infile: path to file with input data (python pickle file)
    :param str outfile: desired path to report file
    :param list report_columns: list of columns for which to report diff's

    :return: absoulte path to report file
    :rtype: str
    """
    with open(infile, 'rb') as inputs:
        root_contents = sorted(pickle.load(inputs).items(), key=lambda x: x[0])

    previous_species = None
    with open(outfile, 'wt') as rfile:
        for i, key_report in enumerate(root_contents):
            key, report = key_report
            species, release = key.split('.')

            # Transofrm the 8th column from dict to list of dict keys
            if isinstance(report[8], dict):
                report[8] = list(report[8].keys())

            rfile.write('species: {}, release: {}\n'.format(species, release))

            # If first entry or species has changed:
            if i == 0 or previous_species != species:
                for column in report_columns:
                    rfile.write('column {}: {}\n'.format(column, report[column]))
            else:
                changes = 0
                for column in report_columns:
                    previous_report = root_contents[i - 1][1]
                    new_values = [e for e in report[column] if e not in previous_report[column]]
                    missing_values = [e for e in previous_report[column] if e not in report[column]]
                    if any([new_values, missing_values]):
                        changes += 1
                        rfile.write('column {}\n'.format(column))
                        if new_values:
                            rfile.write('new_values: {}\n'.format(new_values))
                        if missing_values:
                            rfile.write('missing_values: {}\n'.format(missing_values))
                if changes == 0:
                    rfile.write('No changes.\n')
            rfile.write(SEPARATOR + '\n')

            previous_species = species

        # Finally, also report which values are present in all releases (per column)
        rfile.write(SEPARATOR + '\n')
        rfile.write('Values present in all releases:\n')
        for column in report_columns:
            always_present = set(root_contents[0][1][column])
            for i in range(len(root_contents)):
                always_present = always_present & set(root_contents[i][1][column])
            rfile.write('column {}: {}\n'.format(column, always_present))

    return os.path.abspath(outfile)


class TestGTFVariation(unittest.TestCase):

    def setUp(self):
        self.reference_results_file = __file__.rstrip('.py') + '.reference.txt'
        self.pickle_file = __file__[:-3] + '.pickle'
        self.report_file = __file__[:-3] + '.report.txt'

        # Folder where all gtf files are stored:
        self.gtf_files_dir = os.path.join(functional_test_folder, 'gtf_files/')
        if not os.path.exists(self.gtf_files_dir):
            os.makedirs(self.gtf_files_dir)

        self.releases_list = ensembl.get_release_list()
        self.species_list = ['homo_sapiens', 'mus_musculus']

        self.current_results = []
        self.reference_results = []

        # The is the container with all filtered information. It has
        # 'species.release' keys and its values are extracted results
        self.root_container = {}

    def tearDown(self):
        if not os.path.isfile(self.reference_results_file):
            # Write new reference file:
            with open(self.reference_results_file, 'w') as ref_file:
                for result in self.current_results:
                    ref_file.write(result + '\n')

    def test_reference_variation(self):
        for species in self.species_list:
            for release in self.releases_list:
                t1 = time.time()

                gtf_file = os.path.join(self.gtf_files_dir, '{}.{}.gtf.gz'.format(species, release))
                # If gtf_in file is not preswent, download it:
                if not os.path.isfile(gtf_file):
                    print('Could not find file: {}'.format(gtf_file))
                    print('Donloading it to: {}'.format(self.gtf_files_dir))
                    _ = ensembl.download_annotation(
                        release, species, target_dir=self.gtf_files_dir)

                print("Reading file: {}".format(gtf_file))

                gs = pybedtools.BedTool(gtf_file)
                number_of_segments = gs.count()
                field_data = {name: [] for name in range(9)}

                print("Processing segments: ")
                for i, segment in enumerate(gs):
                    # Display progress in terminal:
                    sys.stdout.write("\r{0:.1f} %".format(i / number_of_segments * 100))
                    sys.stdout.flush()

                    for i, container in field_data.items():
                        # Case for all, but 9th column:
                        if ';' not in segment.fields[i]:
                            container.append(segment.fields[i])
                        # Case for 9th column:
                        else:
                            keys = list(segment.attrs.keys())
                            if not container:
                                container.append({})
                            # Count the number of times a key occures in GTF:
                            for key in keys:
                                if key in container[0]:
                                    container[0][key] += 1
                                else:
                                    container[0][key] = 1

                # Pack the results into report object:
                report = {}
                report['number_of_segments'] = number_of_segments

                for num, container in field_data.items():
                    # Case for all, but 9th column:
                    if isinstance(container[0], str):
                        set_ = sorted(list(set(container)))
                        if len(set_) > 100:
                            # This is probably a column that describes
                            # chromosome coordinates, gene_name etc. Just
                            # report on few example values:
                            report[num] = set_[:3]
                        else:
                            # Otherwise report the whole set of values:
                            report[num] = set_
                    # Case for 9th column
                    else:
                        report[num] = container[0]

                root_key = '{}.{}'.format(species, release)
                self.root_container[root_key] = report

                # Save current version of root_container into pickle file - this
                # is usefull if things break up for some reason
                with open(self.pickle_file, 'wb') as output:
                    pickle.dump(self.root_container, output, pickle.HIGHEST_PROTOCOL)

                print()
                print("This took {0:.2f} minutes to complete".format((time.time() - t1)/60))

        # Call function to generate report file:
        self.report_file = produce_report(self.pickle_file, self.report_file)

        # Confirm that reference file exists!
        self.assertTrue(os.path.isfile(self.reference_results_file))

        # Read the report and the reference file:
        with open(self.report_file) as rep, open(self.reference_results_file) as ref:
            for rep_line, ref_line in zip(rep, ref):
                self.current_results.append(rep_line.strip())
                self.reference_results.append(ref_line.strip())

        # Confirm that self.reference_results and self.current_results have
        # identical content:
        self.assertEqual(set(self.reference_results), set(self.current_results))

if __name__ == '__main__':
    unittest.main()

from .. files import bed

analysis_name = 'group'
analysis_description_short = 'Group multiple BED files'
analysis_description = 'Merge BED6 multiple files with cross link data into one file.'
params_opt = [
    ('outfile', 'string', None, True, 'Output filename.'),
]
params_pos = [('files', None, None, '+', 'List of BED6 files.')]


def run(files, outfile):
    bed.merge_bed(files, outfile)

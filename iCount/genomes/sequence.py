from . import ensembl

# description and parameters needed for the analysis
analysis_name = 'sequence'
analysis_description_short = 'ensembl sequence'
analysis_description = 'Download *.fasta sequence file for a specified ' \
                       'ensembl release and species.'

params_opt = [
    ('release', 'int_range', (ensembl.MAX_RELEASE_SUPPORTED,
                              ensembl.MIN_RELEASE_SUPPORTED,
                              ensembl.MAX_RELEASE_SUPPORTED),
     True, 'Ensembl release version.'),

    ('species', 'string', "homo_sapiens", True, 'Species name.'),

    ('target_dir', 'string', None, False, 'Target directory.'),

    ('target_fname', 'string', None, False, 'Target filename.'),
]

params_pos = []


def get(release, species, target_dir=None, target_fname=None):
    return ensembl.download_sequence(release, species, target_dir=target_dir, target_fname=target_fname)

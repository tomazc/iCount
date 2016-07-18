from . import ensembl

# description and parameters needed for the analysis
analysis_name = 'species'
analysis_description_short = 'ensembl species'
analysis_description = 'Retrive list of available species for a specified ' \
                       'ensembl release.'

# no parameters needed
params_opt = [
    (
        'release', 'int_range', (ensembl.MAX_RELEASE_SUPPORTED,
                                 ensembl.MIN_RELEASE_SUPPORTED,
                                 ensembl.MAX_RELEASE_SUPPORTED), True,
        'Ensembl release version.'
    ),
]
params_pos = []


def get(release):
    """Return list of species included in ensembl release.

    Spieces are ordered lexicographically.

    """
    return ensembl.get_species_list(release)

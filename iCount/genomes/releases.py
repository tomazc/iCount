from . import ensembl

# description and parameters needed for the analysis
analysis_name = 'releases'
analysis_description_short = 'ensembl releases'
analysis_description = 'Retrive list of available releases on ensembl FTP ' \
                       'site.'

# no parameters needed
params_opt = []
params_pos = []


def get():
    """Return list of ensembl releases available for download.

    Releases are ordered in decreasing order, latest release is listed first.

    """
    return ensembl.get_release_list()

import pybedtools

from . import ensembl


def get(release, species):
    return ensembl.download_annotation(release, species)

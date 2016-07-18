from . import ensembl

def get(release, species):
    return ensembl.download_sequence(release, species)

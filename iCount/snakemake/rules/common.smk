








##### Helper functions #####


def get_genome(wildcards):
    return ("{0}".format(samples.loc[wildcards.barcode, "mapto"]))


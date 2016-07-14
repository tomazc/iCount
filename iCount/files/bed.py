import pybedtools

def load(fname):
    return pybedtools.BedTool(fname)

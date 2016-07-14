import pybedtools

def load(fn):
    return pybedtools.BedTool(fn)

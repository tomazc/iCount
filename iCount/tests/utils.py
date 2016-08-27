import os
import tempfile

import pybedtools

from pybedtools import create_interval_from_list


def make_file_from_list(data, bedtool=True):
    """Return path to file with the content from `data` (list of lists)"""
    tfile = tempfile.NamedTemporaryFile(mode='w+', delete=False)
    tfile.close()
    if bedtool:
        pybedtools.BedTool(pybedtools.create_interval_from_list(list_)
                           for list_ in data).saveas(tfile.name)
    else:
        with open(tfile.name, 'wt') as file_:
            for list_ in data:
                file_.write('\t'.join(map(str, list_)) + '\n')
    return os.path.abspath(tfile.name)


def make_list_from_file(fname, fields_separator=None):
    """Read file to list of lists"""
    data = []
    with open(fname) as file_:
        for line in file_:
            data.append(line.strip().split(fields_separator))
    return data


def list_to_intervals(data):
    """Transform list of lists to list of pybedtools.Intervals."""
    return [create_interval_from_list(list_) for list_ in data]


def intervals_to_list(data):
    """Transform list of pybedtools.Intervals to list of lists"""
    return [interval.fields for interval in data]


def reverse_strand(data):
    """Reverses the strand of every element in data. Elements can be
    pybedtools.Intervals or lists with same content as interval.fields."""
    rstrands = ['-' if i[6] == '+' else '+' for i in data]
    if isinstance(data[0], pybedtools.Interval):
        return [create_interval_from_list(
            data[i][:6] + [rstrands[i]] + data[i][7:]) for i in range(len(data))]
    else:
        return [
            data[i][:6] + [rstrands[i]] + data[i][7:] for i in range(len(data))]

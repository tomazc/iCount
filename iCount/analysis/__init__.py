"""
Analysis
========

.. automodule:: iCount.analysis.annotate
   :members:

.. automodule:: iCount.analysis.clusters
   :members:

.. automodule:: iCount.analysis.group
   :members:

.. automodule:: iCount.analysis.kmers
   :members:

.. automodule:: iCount.analysis.peaks
   :members:

.. automodule:: iCount.analysis.rnamaps
   :members:

.. automodule:: iCount.analysis.summary
   :members:

"""

from . import annotate
from . import clusters
from . import group
from . import kmers
from . import peaks
from . import rnamaps
from . import summary


class IntRangeParser:
    def __init__(self, par_name, v_min, v_max):
        self.v_min = v_min
        self.v_max = v_max
        self.par_name = par_name

    def bound(self, v):
        if not self.v_min <= v <= self.v_max:
            new_v = max(self.v_min, min(self.v_max, v))
            print('WARNING: given parameter {:s} value ({:d}) '
                  'is out of range ({:d}..{:d}), '
                  'changed to {:d}.'.format(self.par_name, v, self.v_min,
                                            self.v_max, new_v)
                  )
            return new_v
        return v

    def __call__(self, string):
        return self.bound(int(string))


class FloatRangeParser:
    def __init__(self, par_name, v_min, v_max):
        self.v_min = v_min
        self.v_max = v_max
        self.par_name = par_name

    def bound(self, v):
        if not self.v_min <= v <= self.v_max:
            new_v = max(self.v_min, min(self.v_max, v))
            print('WARNING: given parameter {:s} value ({:f}) '
                  'is out of range ({:f}..{:f}), '
                  'changed to {:f}.'.format(self.par_name, v, self.v_min,
                                            self.v_max, new_v)
                  )
            return new_v
        return v

    def __call__(self, string):
        return self.bound(float(string))


def str_list(string):
    string = string.strip()
    if not string:
        return []
    return [x.strip() for x in string.split(',')]


def str_list_to_str(lst):
    if lst:
        return ",".join([str(x) for x in lst])
    else:
        return []


class IntervalsListParser:
    def __init__(self, par_name, v_min, v_max):
        self.v_min = v_min
        self.v_max = v_max
        self.par_name = par_name

    def bound(self, v):
        if not self.v_min <= v <= self.v_max:
            new_v = max(self.v_min, min(self.v_max, v))
            print('WARNING: given parameter {:s} value ({:d}) '
                  'is out of range ({:d}..{:d}), '
                  'changed to {:d}.'.format(self.par_name, v, self.v_min,
                                            self.v_max, new_v)
                  )
            return new_v
        return v

    def __call__(self, string):
        retl = [x.split(':') for x in string.split(',')]
        retl = [
            (self.bound(int(f)), self.bound(int(t)))
            for f, t in retl
            ]
        assert all(f <= t for f, t in retl)
        return retl


def intervals_to_str(lst):
    return ",".join(['{:d}:{:d}'.format(f, t) for f, t in lst])


def params_to_argparse(subparsers, a):
    # setup name and description of analysis
    import argparse
    parser = subparsers.add_parser(
        a.analysis_name,
        description=a.analysis_description,
        help=a.analysis_description_short,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # setup optional parameters
    for par_name, par_type, par_vals, par_req, par_help in a.params_opt:
        par_name = "-{:s}".format(par_name)
        if par_type == 'int_range':
            par_val, par_val_min, par_val_max = par_vals
            par_val = str(par_val)
            parser.add_argument(par_name, default=par_val,
                                type=IntRangeParser(par_name,
                                                    par_val_min,
                                                    par_val_max),
                                required=par_req,
                                help=par_help)
        elif par_type == 'float_range':
            par_val, par_val_min, par_val_max = par_vals
            par_val = str(par_val)
            parser.add_argument(par_name, default=par_val,
                                type=FloatRangeParser(par_name,
                                                      par_val_min,
                                                      par_val_max),
                                required=par_req,
                                help=par_help)
        elif par_type == 'choice-single':
            par_val = str(par_vals[0])
            par_choices = [str(x) for x in par_vals[1]]
            parser.add_argument(par_name, default=par_val, type=str,
                                choices=par_choices,
                                required=par_req,
                                help=par_help)
        elif par_type == 'file_in' or par_type == 'file_out':
            parser.add_argument(par_name,
                                required=par_req,
                                help=par_help)
        elif par_type == 'str_list':
            par_val = par_vals
            par_val = str_list_to_str(par_val)
            parser.add_argument(par_name, default=par_val, type=str_list,
                                required=par_req,
                                help=par_help)
        elif par_type == 'intervals':
            par_val, par_val_min, par_val_max = par_vals
            par_val = intervals_to_str(par_val)
            parser.add_argument(par_name, default=par_val,
                                type=IntervalsListParser(par_name,
                                                         par_val_min,
                                                         par_val_max),
                                required=par_req,
                                help=par_help)
        elif par_type == 'string':
            par_val = par_vals
            parser.add_argument(par_name, default=par_val,
                                type=str,
                                required=par_req,
                                help=par_help)
        else:
            print('Unhandled parameter type: {:s}'.format(par_type))
            assert False

    # setup positional parameters
    for par_name, par_type, par_mode, par_nargs, par_help in a.params_pos:
        if par_nargs != 1:
            parser.add_argument(par_name, nargs=par_nargs, help=par_help)
        else:
            parser.add_argument(par_name, help=par_help)

    return parser

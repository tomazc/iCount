#!/usr/bin/env python

"""
Command line interface for iCount Python package.
#################################################

This module offers CL interface to iCount Python package. To separate Python
package from CLI and to avoid code and docstring duplication, CL commands are
semi-automatically created from functions that already exist in the package.
For this automation to work, some rules need to respected when writing functions
that are later exposed to CLI.


.. autofunction:: iCount.scripts.make_parser_from_function

"""

import re
import sys
import logging
import inspect
import argparse
import traceback
import collections

import iCount
from iCount import logger

from sphinx.ext.napoleon.docstring import NumpyDocstring

LOGGER = logging.getLogger(__name__)


########################################
# Functions for atomated CLI generation:
########################################

def _list_str(string, separator=','):
    """Convert comma separated string to list"""
    return string.strip('{:s} \n\t'.format(separator))

VALID_TYPES = {
    'bool': bool,
    'str': str,
    'int': int,
    'float': float,
    'list_str': _list_str,
}


def _extract_parameter_data(function):
    """
    Extract name, type, description, etc. of each parameter in ``function``.

    The returned object should have the form::

        data = {
            'param1': {
                name = 'param1',
                type = 'str',
                help = 'The first parameter. Just for demo.',
            },
            'param2': {
                name = '--param2',
                type = 'int',
                help = 'Some optional input number...',
            },
        }

    Every parameter in returned object can have the following entries:

        * name - the name of parameter, preceeded by '--' if it is optional
        * default - the default value (only for optional parameters). Extracted
          from function signature.
        * type - type of parameter, extracted from function docstring. If not
          set, it is assumed to be str.
        * help - parameter description, extracted from function docstring. If
          not found, set to 'No description'

    """

    # Use OrderedDict to keep the order of parameters
    data = collections.OrderedDict()

    # Get function docstring and convert from Numpy to rst docstring style:
    function_docstring = str(NumpyDocstring(inspect.getdoc(function))) or ''

    # Determine weather argument is positional or optional: positional arguments
    # are the ones that do not have default values defined:
    params = inspect.getargspec(function).args
    default_values = list(inspect.getargspec(function).defaults or [])
    pos_opt_border = len(params) - len(default_values)
    optional = {par: default for par, default in
                zip(params[pos_opt_border:], default_values)}

    for param in params:
        if param in optional:
            data[param] = {'name': '--' + param, 'default': optional[param], 'metavar': ''}
        else:
            data[param] = {'name': param}

        # Extract help (parameter description) and type from function docstring:
        regex_help = r'.*:param {}: (.+?)\n.*'.format(param)
        regex_type = r'.*:type {}: (.+?)\n.*'.format(param)
        # Use DOTALL flag - this way '.' matches also newline characters:
        match_help = re.match(regex_help, function_docstring, flags=re.DOTALL)
        match_type = re.match(regex_type, function_docstring, flags=re.DOTALL)

        if match_help and match_type:
            param_type = match_type.group(1).strip()
            if param_type not in VALID_TYPES:
                raise ValueError('Invalid type {} in function {}'.format(
                    param_type, function.__name__))
            data[param]['type'] = VALID_TYPES[param_type]
            data[param]['help'] = match_help.group(1).strip().rstrip('.')

            if param_type == 'bool':
                data[param]['action'] = 'store_true'
                # If action == store_true, than `type` needs to be removed.
                data[param].pop('type')
                data[param].pop('metavar')
            if param_type == 'list_str':
                data[param]['nargs'] = '+'

        else:
            # TODO: raise ValueError or just make some arbitrary description?
            # raise ValueError("Bad docstring for parameter {}.".format(param))
            data[param]['help'] = 'No description'

    return data


def make_parser_from_function(function, subparsers, module=None, only_func=False):
    """
    Make a argparse parser from ``function`` and add it to ``subparsers``.

    A Python function is exposed in CLI by calling ``make_parser_from_function``.
    Example call for exposing function ``iCount.analysis.peaks.run`` that makes
    peaks analysis::

        make_parser_from_function(iCount.analysis.peaks.run, subparsers)

    What happens is such call?
        * CLI command ``iCount peaks`` is created, with description, positional
          and optional arguments exactly the same as the ones defined in
          function.
        * name of the command (peaks) equals to the name of the module where
          function is defined.
        * CLI help message is sourced from module's docstring
        * Positional and optional arguments (and default values) are sourced
          from function signature
        * Help text for each of the arguments is sourced from function
          docstring. For this to work functions needs to follow Nummpy docstrig
          formatting. All parameters should have meaningful description. All
          parameters should also have type equal to one of the keys in
          VALID_TYPES.
        * When function is executed, stdout logger with level INFO is registerd.
          This prints descriptive messages to CL. Function should therefore log
          its inputs (use iCount.logger.log_inputs function), outputs and most
          important steps. If there are no errors/execptions command exits with
          0 exit status code.
        * If error/exception occurs, is it caught and stack is printed to logger
          with ERROR level. Also, failed CLI command has status code 1.
        * CLI exposed function should return Metric object that is instance of
          iCount.metrics.Metric class. Function should set descriptive
          attributes to the Metric object for analysis results and analysis
          tatistics for inter-experimetal comparison.
        * Each command also gets these additional arguments:
            * ``--stdout_log`` - Threshold value (0-50) for logging to stdout. If 0
              logging to stdout if turned OFF.
            * ``--file_log`` - Threshold value (0-50) for logging to file. If 0
              logging to stdout if turned OFF.
            * ``--file_logpath`` - Path to log file.
            * ``--results_file`` - File into which to store Metrics (result object).
            * ``--help`` - Help message fot a command

            They control the level and location
            of log inputs as well as string the result objects.


    Exceptional cases:

        * In some cases function that performs the work is only *imported* in the
          correct module ("exposed module"), but the actual function definition is
          loctated somewhere else ("source module"). It that case, one can use
          ``module`` parameter with "source module" value. In this case, the
          command name and CLI docstring will be defined from "exposed module" but
          function, default values and parameter descriptions will be sourced from
          "source module".

        * In some cases, there will be more CLI exposed functions in the same
          module. In such case, use set ``only_func`` parameter to ``True``. This will
          use function name for CLI command name and use function docstring (form
          beggining until "Parameters" section) for CLI help message.

    """

    if module:
        # Command name is determined from module name:
        name = module.__name__.split('.')[-1]
        # Description is determined from module docstring:
        description = inspect.getdoc(module)
    elif only_func:
        # If only_func=True than name of the command equals function name
        name = function.__name__
        # Take only the docstring until the 'Parameters' section:
        desc = inspect.getdoc(function).split('\n')
        idx = [i for i, d in enumerate(desc) if d == 'Parameters' or d == 'Returns'][0]
        description = '\n'.join(desc[:idx])
    else:
        # Command name is determined from module name:
        name = inspect.getmodule(function).__name__.split('.')[-1]
        # Description is determined from module docstring:
        description = inspect.getdoc(inspect.getmodule(function))

    if not description:
        description = 'No description provided.'
    # Short description is determined from first line of module docstring:
    short_description = description.split('\n')[0]

    parser = subparsers.add_parser(
        name,
        description=description,
        help=short_description,
        formatter_class=argparse.RawTextHelpFormatter  # ArgumentDefaultsHelpFormatter
    )

    # Add each of the function parameter as CLI argument:
    params = _extract_parameter_data(function)
    for values in params.values():
        name = values.pop('name')
        parser.add_argument(name, **values)

    parser.add_argument('-S', '--stdout_log', default=logging.INFO, type=int, metavar='',
                        help='Threshold value (0-50) for logging to stdout. If 0,'
                        ' logging to stdout if turned OFF.')
    parser.add_argument('-F', '--file_log', default=0, type=int, metavar='',
                        help='Threshold value (0-50) for logging to file. If 0,'
                        ' logging to file if turned OFF.')
    parser.add_argument('-P', '--file_logpath', default=None, type=str, metavar='',
                        help="Path to log file.")
    parser.add_argument('-M', '--results_file', default=None, type=argparse.FileType('wt'),
                        metavar='', help="File into which to store Metrics.")

    # Connect parser with function to execute:
    parser.set_defaults(func=function)

# #########################################################################
# #########################################################################


def main():

    #####################
    # Define root parser:
    #####################

    root_parser = argparse.ArgumentParser(
        description=inspect.getdoc(iCount).split('\n..')[0],
        formatter_class=argparse.RawTextHelpFormatter  # ArgumentDefaultsHelpFormatter,
    )
    # Parser (general) arguments
    root_parser.add_argument('-v', '--version', action='version',
                             version='%%(prog)s %s' % iCount.__version__)

    subparsers = root_parser.add_subparsers(title='Commands', metavar=' ')
    ##################
    # Define commands:
    ##################

    # Genomes:
    make_parser_from_function(
        iCount.genomes.ensembl.get_release_list, subparsers, only_func=True)
    make_parser_from_function(
        iCount.genomes.ensembl.get_species_list, subparsers, only_func=True)
    make_parser_from_function(
        iCount.genomes.ensembl.download_annotation, subparsers, only_func=True)
    make_parser_from_function(
        iCount.genomes.ensembl.download_sequence, subparsers, only_func=True)
    make_parser_from_function(
        iCount.genomes.segment.get_regions, subparsers)
    make_parser_from_function(
        iCount.genomes.genes.get_genes, subparsers)

    # Demultiplex and mapping:
    make_parser_from_function(iCount.demultiplex.run, subparsers)
    make_parser_from_function(iCount.demultiplex.remove_adapter, subparsers,
                              module=iCount.externals.cutadapt)
    make_parser_from_function(iCount.mapping.mapindex.run, subparsers,
                              module=iCount.mapping.mapindex)
    make_parser_from_function(iCount.mapping.map.run, subparsers, module=iCount.mapping.map)
    make_parser_from_function(iCount.mapping.xlsites.run, subparsers)

    # # Analysis:
    make_parser_from_function(
        iCount.analysis.annotate.annotate_cross_links, subparsers)
    make_parser_from_function(
        iCount.analysis.clusters.run, subparsers)
    make_parser_from_function(
        iCount.analysis.group.run, subparsers, module=iCount.analysis.group)
    make_parser_from_function(
        iCount.analysis.peaks.run, subparsers)
    make_parser_from_function(
        iCount.analysis.summary.make_summary_report, subparsers)

    def verbose_help(mode):
        if mode == 'txt':
            print(root_parser.format_help() + '\n')
            for name, parser in root_parser._action_groups[2]._actions[2].choices.items():
                if name == 'man':
                    continue
                print(name)
                print('=' * len(name) + '\n')
                print(parser.format_help() + '\n')
        else:
            raise NotImplementedError('"{}" not supported.'.format(mode))

    # Add the man command:
    parser = subparsers.add_parser('man', help='Print help for all commands.')
    parser.add_argument('--mode', type=str, default='txt', metavar='')
    parser.set_defaults(func=verbose_help)

    #############################
    # Parse and execute commands:
    #############################

    # Parse agrs and run:
    parsed_args = root_parser.parse_args()

    if not vars(parsed_args):
        root_parser.print_help()
        root_parser.exit(1)
    try:
        args = vars(parsed_args)
        func = args.pop('func')

        # Stdout logger
        stdout_loglevel = args.pop('stdout_log', 20)
        is_on = not stdout_loglevel == 0
        logger.log_to_stdout(is_on=is_on, level=stdout_loglevel)

        # File logger
        file_loglevel = args.pop('file_log', 0)
        file_logpath = args.pop('file_logpath', None)
        is_on = all([file_loglevel, file_logpath])
        logger.log_to_file(is_on=is_on, level=file_loglevel, path=file_logpath)

        # Switch to output results to a file:
        results_file = args.pop('results_file', None)

        # Execute the wrapped function:
        command = 'iCount ' + ' '.join(sys.argv[1:])
        LOGGER.info("Executing the following command: %s" % command)
        result_object = func(**args)

        # Save results to results_file. If not given, don't do anything (don't print the result!)
        if results_file:
            print(result_object, file=results_file)
        sys.exit(0)

    except Exception as exception:
        exception_message = exception.args[0]
        exception_type = exception.__class__.__name__

        LOGGER.error('[%s] %s' % (exception_type, exception_message))
        for line in traceback.format_list(traceback.extract_tb(exception.__traceback__)):
            LOGGER.error(line)
        sys.exit(1)

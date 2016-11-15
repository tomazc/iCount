""".. Line to protect from pydocstyle D205, D400.

Logging
=======

A parent logger for all modules in iCount library is created.

Two utility functions (``log_to_stdout``, ``log_to_file``) ease the creation of
most frequently used handlers.

Also, ``log_inputs`` function is created. It can be called from a function and
logs all the function input params. This is especially usefull for analyses that
are logged to file (for reproduction reasons).

Additonally, sys excepthook mechanism is modified: All python exceptions are
handled by function, stored in ``sys.excepthook.`` By rewriting the default
implementation, we can modify it for our puruses - to log all uncaught
exceptions.

Note#1: Modified behaviour (logging of all uncaught exceptions) applies
only when runing in non-interactive mode.

Note#2: Any exception can be caught/uncaought and it can happen in
interactive/non-interactive mode. This makes 4 different scenarios.
The sys.excepthook modification takes care of uncaught exceptions in
non-interactive mode. In interactive mode, user is notified directly
if exception is raised. If exception is caught and not reraised, it
should  be logged somehow, since it can provide valuable information
for  developer when debugging. Therefore, we should use the following
convention for logging: "Exceptions are explicitly logged
only when they are caught and not re-raised."

"""
import os
import sys
import logging
import inspect
import iCount

LEVEL_MAP = {
    "DEBUG": logging.DEBUG,
    "INFO": logging.INFO,
    "WARNING": logging.WARNING,
    "ERROR": logging.ERROR,
    "EXCEPTION": logging.ERROR,
    "CRITICAL": logging.CRITICAL,
}

# Create root logger:
LOGGER_NAME = __name__.split('.')[0]
ROOT_LOGGER = logging.getLogger(LOGGER_NAME)
# Set level to lowest threshold possible, since specific handlers can
# increase the level, if neccessary:
ROOT_LOGGER.setLevel(logging.DEBUG)  # lowest possible value

# Create "do-nothing" handler. Why? If both (file and stdout handlers) are turned off,
# having NullHandler:
# (in py27) - prevents message: "No handlers could be found for logger XXXX".
# (in py3+) - prevents messages higher than warning to be printed to stderr.
ROOT_LOGGER.addHandler(logging.NullHandler())


def _configure_handler(handler, is_on=None, level=None):
    """
    Helper function for configuring handlers binded to ROOT_LOGGER.

    Parameters
    ----------
    is_on: bool
        Switch to turn handler ON or OFF.
    level: int
        Logging threshold level (should be integer between 1 and 50, or a key
        from LEVEL_MAP).

    Returns
    -------
    None
    """
    if is_on is not None:
        if isinstance(is_on, bool):
            if is_on:
                ROOT_LOGGER.addHandler(handler)
            else:
                ROOT_LOGGER.removeHandler(handler)
        else:
            raise ValueError("Wrong type of 'is_on' parameter: only True/False/None posssible.")

    if level is not None:
        if isinstance(level, str) and level.upper() in LEVEL_MAP.keys():
            handler.setLevel(LEVEL_MAP[level.upper()])
        elif isinstance(level, int) and level >= 0 and level <= 50:
            handler.setLevel(level)
        else:
            raise ValueError("Wrong value of 'level' parameter.")


def log_to_stdout(is_on=True, level=logging.INFO):
    """Configure logging to stdout.

    Parameters
    ----------
    is_on: bool
        Switch to turn handler ON or OFF.
    level: int
        Logging threshold level (should be integer between 1 and 50, or a key
        from LEVEL_MAP).

    Returns
    -------
    None
    """
    stdout_handler = logging.StreamHandler()
    stdout_handler.setFormatter(logging.Formatter(fmt='%(message)s'))

    _configure_handler(stdout_handler, is_on=is_on, level=level)


def log_to_file(is_on=False, level=logging.WARNING, path=None):
    """Configure logging to file.

    Each time this function is called, new file handler is created.

    Parameters
    ----------
    is_on: bool
        Switch to turn handler ON or OFF.
    level: int
        Logging threshold level (should be integer between 1 and 50, or a key
        from LEVEL_MAP).
    path: str
        Path to logfile.

    Returns
    -------
    None
    """
    if not path:
        path = os.path.join(iCount.OUTPUT_ROOT, 'iCount.log')
    file_handler = logging.FileHandler(path)
    formatter = logging.Formatter(
        fmt='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S')
    file_handler.setFormatter(formatter)

    _configure_handler(file_handler, is_on=is_on, level=level)


def log_inputs(logger, level=logging.INFO):
    """
    Log the calling function input params to `logger` with `level` severity.

    Parameters
    ----------
    logger: logging.Logger
        Loger instance to which to log.
    level: int
        Logging threshold level (should be integer between 1 and 50, or a key
        from LEVEL_MAP).

    Returns
    -------
    None
    """
    # Get frame of calling function and function object
    function_frame = inspect.currentframe().f_back
    function_object = function_frame.f_globals[function_frame.f_code.co_name]

    args = [(arg, function_frame.f_locals[arg]) for arg in
            inspect.signature(function_object).parameters]

    logger.log(level, "Input parameters for function '{}' in {}".format(
        function_object.__name__,  # function name
        function_object.__module__,  # file/module name
    ))
    for arg_name, arg_value in args:
        logger.log(level, "    {}: {}".format(arg_name, arg_value))


def _log_progress(new_ratio, old_ratio, logger, decimals=2):
    """
    If progress has insreased sufficiently, log it to ``logger``.

    If ``new_ratio``, rounded to ``decimals`` differs from ``old_ratio``, log to
    logger with INFO level and return rounded new_ratio. Else return unmodified
    ``old_ratio``.
    """
    new_ratio = round(new_ratio, decimals)
    if new_ratio != old_ratio:
        logger.info('%s', '{}%'.format(new_ratio * 100))
        return new_ratio
    else:
        return old_ratio


def _log_all_uncaught_exceptions(exc_type, exc_value, exc_traceback):
    """Log all uncaught exceptions in non-interactive mode.

    All python exceptions are handled by function, stored in
    ``sys.excepthook.`` By rewriting the default implementation, we
    can modify handling of all uncaught exceptions.

    Warning: modified behaviour (logging of all uncaught exceptions)
    applies only when runing in non-interactive mode.

    """
    # ignore KeyboardInterrupt
    if not issubclass(exc_type, KeyboardInterrupt):
        ROOT_LOGGER.error("", exc_info=(exc_type, exc_value, exc_traceback))

    sys.__excepthook__(exc_type, exc_value, exc_traceback)
    return


# Rewrite the default implementation os sys.excepthook to log all
# uncaught exceptions:
sys.excepthook = _log_all_uncaught_exceptions

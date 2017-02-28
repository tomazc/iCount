""".. Line to protect from pydocstyle D205, D400.

Examples
========

Provide a set of example bash scripts.

.. autofunction:: iCount.examples.run


"""

import os
import shutil
import logging


LOGGER = logging.getLogger(__name__)


def run(out_dir='.'):
    """
    Create an examples subfolder with example scripts.

    This will create an examples subfolder in current working directory and
    will copy bash scripts needed to run the iCount pipeline on a few
    examples (for now, only the hnRNP C data from Konig et al.).

    Parameters
    ----------
    out_dir : str
        Directory to which example scripts should be copied.

    Returns
    -------
    str
        Path of folder to where examples scripts were copied.

    """
    examples_folder = os.path.abspath(os.path.join(out_dir, 'examples'))
    if os.path.exists(examples_folder):
        raise FileExistsError('Examples folder already exists.')

    try:
        os.makedirs(examples_folder)
        LOGGER.info('Setting up folder with examples.')
    except OSError:
        raise OSError('Error creating examples folder.')

    cur_folder = os.path.dirname(os.path.abspath(__file__))
    for script in ['hnRNPC.sh', 'hnRNPC_reduced.sh']:
        src_fn = os.path.join(cur_folder, script)
        dst_fn = os.path.join(examples_folder, script)
        LOGGER.info('   copying example script %s', script)
        try:
            shutil.copy(src_fn, dst_fn)
        except OSError:
            message = 'Error copying example script from: {} to: {}'.format(src_fn, dst_fn)
            raise OSError(message)

    LOGGER.info('Done. Check example bash scripts in subfolder \'examples\'')
    return examples_folder

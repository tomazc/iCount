import os
import shutil

analysis_name = 'examples'
analysis_description_short = 'setup examples folder'
analysis_description = 'Create a subfolder with example scripts for running ' \
                       'iCount pipeline for iCLIP data.'

params_opt = []
params_pos = []


def run():
    """Creates an examples subfolder with example scripts.

    This will create an examples subfolder in current working directory and
    will copy bash scripts needed to run the iCount pipeline on a few
    examples (for now, only the hnRNP C data from Konig et al.).

    :return: None
    """

    examples_folder = os.path.abspath(os.path.join('.', 'examples'))
    if os.path.exists(examples_folder):
        print('Examples folder already exists: {:s}'.format(examples_folder))
        print('Scripts will NOT be copied into it.')
        print('Make sure subfolder \'examples\' does not exist before '
              'rerunning the command.')
        return

    try:
        os.makedirs(examples_folder)
    except OSError as e:
        print('Error creating examples folder {:s}'.format(e))
        return
    print('Setting up folder with examples.')

    cur_folder = os.path.dirname(os.path.abspath(__file__))
    for script in ['hnRNPC.sh', 'hnRNPC_reduced.sh']:
        src_fn = os.path.join(cur_folder, script)
        dst_fn = os.path.join(examples_folder, script)
        print('   copying example script {:s}'.format(script))
        try:
            shutil.copy(src_fn, dst_fn)
        except OSError as e:
            print('Error copying example script: ', e)
            print('   from: {:s}'.format(src_fn))
            print('   to: {:s}'.format(dst_fn))
            return

    print('Check example bash scripts in subfolder \'examples\'')
    return examples_folder

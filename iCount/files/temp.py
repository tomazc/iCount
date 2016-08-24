import os
import gzip
import shutil
import tempfile

import iCount


def decompress_to_tempfile(fname, context='misc'):
    if fname.endswith('.gz'):
        tmp_dir = os.path.join(iCount.tmp_root, context)
        if not os.path.exists(tmp_dir):
            os.makedirs(tmp_dir)

        suffix = '_{:s}'.format(os.path.basename(fname))
        fout = tempfile.NamedTemporaryFile(suffix=suffix, dir=tmp_dir,
                                           delete=False)
        fin = gzip.open(fname, 'r')
        shutil.copyfileobj(fin, fout)
        fin.close()
        fout.close()
        print(fout.name)
        return fout.name

    return fname

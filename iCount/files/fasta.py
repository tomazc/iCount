"""
FASTA
-----

Reading `FASTA`_ files.

"""

def read_fasta(fasta_file):
    """
    Read fasta file and return list

    The retuned list (named data) has the following structure::

        fasta = [
            [header1, sequence1],
            [header2, sequence2],
            [header2, sequence2],
            ...
        ]

    """
    data = []
    with open(fasta_file) as ffile:
        for line in ffile:
            line = line.strip()
            if line.startswith('>'):
                data.append([line])
            else:
                if len(data[-1]) == 1:
                    data[-1].append(line)
                else:
                    data[-1][1] = data[-1][1] + line

    return data

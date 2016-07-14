import os
import cProfile


# time and profile iCLIP flows
import iCount.analysis.iCLIP

in_fastq_fname = os.path.join(
        iCount.storage_root,
        '20101116_LUjh03/SLX-2605.CRIRUN_501.s_4.sequence.txt.gz'
)

# extract
# list all experiment barcodes
barcodes = [
    ('NNNGGTTNN', ''),
    ('NNNTTGTNN', ''),
    ('NNNCAATNN', ''),
    ('NNNACCTNN', ''),
    ('NNNGGCGNN', ''),
]
adapter = 'AGATCGGAAGAGCGGTTCAG'  # 'AGATCGGAAG_1,AGCGGTTCAG_2''

map_to = [
    'hg19',
    'hg19',
    'hg19',
    'hg19',
    'hg19',
]

cProfile.run('iCount.analysis.iCLIP.process_lib(in_fastq_fname, barcodes, '
             'adapter, map_to)', sort="cumtime")

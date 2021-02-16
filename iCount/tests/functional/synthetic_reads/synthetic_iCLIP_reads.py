#!/usr/bin/python2.7

import random
import os
import gzip

### Script to generate synthetic reads to test demultiplexing on iMaps
#
#   python synthetic_iCLIP_reads.py fastq_file_path.fq
#
#
#
#
#

'''
## Get fasta
# From the chr 20 I'll get synthetic reads from position chr20_3470487_4016368
# chr20_2652426_2664336.bed (11kb 2 genes both orientations lots of snora snords) 
bedtools getfasta -fi homo_sapiens_chr20.fa -bed chr20_2652426_2664336.bed > chr20_2652426_2664336.fasta 
/Users/mozosi/Programs/art_bin_MountRainier/art_illumina -ss HS20 --noALN -i "/Users/mozosi/Dropbox (UCL-MN Team)/GitHub/iCount/iCount/tests/functional/synthetic_reads/chr20_2652426_2664336.fasta" -l 100 -f 3 -o "/Users/mozosi/Desktop/synthetic_reads/chr20_2652426_2664336_1" --rndSeed 111
/Users/mozosi/Programs/art_bin_MountRainier/art_illumina -ss HS20 --noALN -i "/Users/mozosi/Dropbox (UCL-MN Team)/GitHub/iCount/iCount/tests/functional/synthetic_reads/chr20_2652426_2664336.fasta" -l 100 -f 3 -o "/Users/mozosi/Desktop/synthetic_reads/chr20_2652426_2664336_2" --rndSeed 222
/Users/mozosi/Programs/art_bin_MountRainier/art_illumina -ss HS20 --noALN -i "/Users/mozosi/Dropbox (UCL-MN Team)/GitHub/iCount/iCount/tests/functional/synthetic_reads/chr20_2652426_2664336.fasta" -l 100 -f 3 -o "/Users/mozosi/Desktop/synthetic_reads/chr20_2652426_2664336_3" --rndSeed 333

# Include some reads in snord86 chr20_2656093_2656195  (TTA)
bedtools getfasta -fi homo_sapiens_chr20.fa -bed chr20_2656093_2656195.bed > chr20_2656093_2656195.fasta 
/Users/mozosi/Programs/art_bin_MountRainier/art_illumina -ss HS20 --noALN -i "/Users/mozosi/Desktop/synthetic_reads/chr20_2656093_2656195.fasta" -l 100 -f 20 -o "/Users/mozosi/Desktop/synthetic_reads/chr20_2656093_2656195" --rndSeed 111

# intron-exon reads in NOP56 chr20_2656398_2656507 (AGG)
bedtools getfasta -fi homo_sapiens_chr20.fa -bed chr20_2656398_2656507.bed > chr20_2656398_2656507.fasta 
/Users/mozosi/Programs/art_bin_MountRainier/art_illumina -ss HS20 --noALN -i "/Users/mozosi/Desktop/synthetic_reads/chr20_2656398_2656507.fasta" -l 100 -f 100 -o "/Users/mozosi/Desktop/synthetic_reads/chr20_2656398_2656507" --rndSeed 111


## Merge
cat chr20_2652426_2664336_1_NNNN_GTAAC_NNN_NN_ATT.fq chr20_2652426_2664336_2_NNNN_GTAAC_NNN_NN_AGG.fq chr20_2652426_2664336_3_NNNN_GTAAC_NNN_NN_TTA.fq chr20_2656093_2656195_NNNN_GTAAC_NNN_NN_TTA.fq chr20_2656398_2656507_NNNN_GTAAC_NNN_NN_TTA.fq chr20_2656398_2656507_NNNN_GTAAC_NNN_NN_AGG.fq > merge_chr20_2652426_2664336_GTAAC.fq
cat merge_chr20_2652426_2664336_GTAAC.fq | wc -l
gzip merge_chr20_2652426_2664336_GTAAC.fq

HiSeq 2000 (100bp)



    ====================ART====================
             ART_Illumina (2008-2016)
          Q Version 2.5.8 (June 6, 2016)
     Contact: Weichun Huang <whduke@gmail.com>
    -------------------------------------------

Warning: your simulation will not output any ALN or SAM file with your parameter settings!
                  Single-end Simulation

Total CPU time used: 849.236

The random seed for the run: 666

Parameters used during run
	Read Length:	100
	Genome masking 'N' cutoff frequency: 	1 in 100
	Fold Coverage:            1X
	Profile Type:             Combined
	ID Tag:

Quality Profile(s)
	First Read:   HiSeq 2000 Length 100 R1 (built-in profile)

Output files

  FASTQ Sequence File:
	single_dat.fq


sed -n 1,10p single_dat.fq > synthetic_iCLIP_reads_10000.fq

NNNNGTAACNNN	NNATT
NNNNGTAACNNN	NNAGG
NNNNGTAACNNN	NNTTA
NNNNGTAACNNN	NNTGC
NNNNGTAACNNN	NNCTG
NNNNGTAACNNN	NNCGT
NNNNGTAACNNN	NNGTC
NNNNGTAACNNN	NNGGA
'''


### Imputs
# fastq file
fastq_file_path="chr20_2656398_2656507.fq"
fastq_file_name = os.path.basename(fastq_file_path)
fastq_file_name = fastq_file_name.replace(".fq", "")


# Synthetic 5' barcode
barcode5 = "GTAAC"
number_random_upstream = 4
number_random_downstream = 3

# Synthetic 5' barcode
barcode3 = "AGG"
number_random_barcode3 = 2

# Synthetic 3' Illumina adapter
adapter3 = "AGATCGGAAGAGCGGTTCAG"
adapter3_qual = "FFFFFFFFFFFFFFFFFFFF"


# Output name including the 5'barcode and UMI lenght
print (fastq_file_name)
number_N_upstream = 'N' * number_random_upstream
number_N_downstream = 'N' * number_random_downstream
number_N_barcode3 = 'N' * number_random_barcode3
modified_fastq_file_out="%s_%s_%s_%s_%s_%s.fq" % (fastq_file_name, number_N_upstream, barcode5, number_N_downstream, number_N_barcode3, barcode3)

print (modified_fastq_file_out)




# fastq parser class
class ParseFastQ(object):
    """Returns a read-by-read fastQ parser analogous to file.readline()"""

    def __init__(self, filePath, headerSymbols=['@', '+']):
        """Returns a read-by-read fastQ parser analogous to file.readline().
        Exmpl: parser.next()
        -OR-
        Its an iterator so you can do:
        for rec in parser:
            ... do something with rec ...

        rec is tuple: (seqHeader,seqStr,qualHeader,qualStr)
        """
        if filePath.endswith('.gz'):
            self._file = gzip.open(filePath)
        else:
            self._file = open(filePath, 'rU')
        self._currentLineNumber = 0
        self._hdSyms = headerSymbols

    def __iter__(self):
        return self

    def next(self):
        return self.__next__()

    def __next__(self):
        """Reads in next element, parses, and does minimal verification.
        Returns: tuple: (seqHeader,seqStr,qualHeader,qualStr)"""
        # ++++ Get Next Four Lines ++++
        elemList = []
        for i in range(4):
            line = self._file.readline()
            self._currentLineNumber += 1  ## increment file position
            if line:
                elemList.append(line.strip('\n'))
            else:
                elemList.append(None)

        # ++++ Check Lines For Expected Form ++++
        trues = [bool(x) for x in elemList].count(True)
        nones = elemList.count(None)
        # -- Check for acceptable end of file --
        if nones == 4:
            raise StopIteration
        # -- Make sure we got 4 full lines of data --
        assert trues == 4, \
            "** ERROR: It looks like I encountered a premature EOF or empty line.\n\
            Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (
            self._currentLineNumber)
        # -- Make sure we are in the correct "register" --
        assert elemList[0].startswith(self._hdSyms[0]), \
            "** ERROR: The 1st line in fastq element does not start with '%s'.\n\
            Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (
            self._hdSyms[0], self._currentLineNumber)
        assert elemList[2].startswith(self._hdSyms[1]), \
            "** ERROR: The 3rd line in fastq element does not start with '%s'.\n\
            Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (
            self._hdSyms[1], self._currentLineNumber)
        # -- Make sure the seq line and qual line have equal lengths --
        assert len(elemList[1]) == len(elemList[3]), "** ERROR: The length of Sequence data and Quality data of the last record aren't equal.\n\
               Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (
        self._currentLineNumber)

        # ++++ Return fatsQ data as tuple ++++
        return tuple(elemList)







def modify_read(fastq_file_path, barcode5):

    parser = ParseFastQ(fastq_file_path)
    counter_reads = 0

    for record in parser:

        # Initialise counter of total reads
        counter_reads += 1

        # Define each element on a fastq file
        header = record[0]
        seq = record[1]
        header2 = record[2]
        qual = record[3]


        random_upstream = ''.join(random.choice('ACTG') for _ in range(number_random_upstream))
        random_downstream = ''.join(random.choice('ACTG') for _ in range(number_random_downstream))

        len_barcode5 = len(barcode5) + number_random_upstream + number_random_downstream

        barcode5_qual = 'F' * len_barcode5  ## All the lenght including randomer

        len_barcode3 = len(barcode3) + number_random_barcode3
        random_barcode3 = ''.join(random.choice('ACTG') for _ in range(number_random_barcode3))
        barcode3_qual = 'F' * len_barcode3

        # print ("Original read")
        # print (seq)
        # print qual
        # print "\n"

        seq = random_upstream + barcode5 + random_downstream + seq + random_barcode3 + barcode3 + adapter3
        qual = barcode5_qual + qual + barcode3_qual + adapter3_qual

        # print "Modified read"
        # print seq
        # print qual
        # print "\n"

        f = open(modified_fastq_file_out, 'a+')
        f.write("%s\n%s\n%s\n%s\n" % (header, seq, header2, qual))
        f.close()


modify_read(fastq_file_path, barcode5)
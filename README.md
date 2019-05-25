[![Build Status](https://travis-ci.org/tomazc/iCount.svg?branch=master)](https://travis-ci.org/tomazc/iCount)
[![PyPI version](https://badge.fury.io/py/iCount.svg)](https://badge.fury.io/py/iCount)
[![codecov](https://codecov.io/gh/tomazc/iCount/branch/master/graph/badge.svg?token=JhUJ66rnJ3)](https://codecov.io/gh/tomazc/iCount)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/b77d104b59a74946bf8905f82dd381e4)](https://www.codacy.com/app/tomazc/iCount?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=tomazc/iCount&amp;utm_campaign=Badge_Grade)
[![CLA assistant](https://cla-assistant.io/readme/badge/tomazc/iCount)](https://cla-assistant.io/tomazc/iCount)
[![Docker Automated build](https://img.shields.io/docker/automated/jrottenberg/ffmpeg.svg)](https://hub.docker.com/r/tomazc/icount/)
[![Documentation Status](https://readthedocs.org/projects/icount/badge/?version=latest)](http://icount.readthedocs.io/en/latest/?badge=latest)

# iCount: protein-RNA interaction analysis

iCount is a Python module and associated command-line interface (CLI),
which provides all the commands needed to process iCLIP data on 
protein-RNA interactions and generate:
 
+ demultiplexed and adapter-trimmed FASTQ files,
+ BAM files with mapped iCLIP reads,
+ identified protein-RNA cross-linked sites, saved to BED files,
+ peaks composed of statistically significant cross-linked sites, saved to BED files,
+ clusters of significant cross-linked sites, saved to BED files,
+ grouping of individual replicate experiments,
+ RNAmap generation showing the positional distribution of cross-linked sites relative to genomic landmarks,
+ kmer enrichment analysis,
+ and other.

You may start with the [tutorial](http://icount.readthedocs.io/en/latest/tutorial.html) or dive into the 
[documentation](http://icount.readthedocs.io/en/latest/index.html).


## Authors

iCount is developed and supported by [Tomaž Curk](http://curk.info) from the 
[Bioinformatics Laboratory](http://biolab.si) at the [University of Ljubljana](http://www.uni-lj.si), 
[Faculty of Computer and Information Science](http://www.fri.uni-lj.si) and in collaboration with 
the laboratory of [Jernej Ule](http://ulelab.info).

The development started in late 2008 when Tomaž Curk and Gregor Rot
wrote a first prototype of iCount. A lot has happend since then. For details and full acknowledgments,
see the [how to cite](http://icount.readthedocs.io/en/latest/cite.html) section in the documentation.


## Contributing

Contributions (pull requests) are welcome! Please submit your contributions by following the
[guidelines](http://icount.readthedocs.io/en/latest/contributing.html).


## Reporting problems

Use the [issues page](https://github.com/tomazc/iCount/issues) to report issues and request enhancements. Before submitting, please check the list of already reported issues. Also, have a look at the [FAQ](http://icount.readthedocs.io/en/latest/faq.html) to see if the problem has been addressed already. 

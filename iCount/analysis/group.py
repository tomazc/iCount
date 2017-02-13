""".. Line to protect from pydocstyle D205, D400.

Group BED files
---------------

Merge multiple BED files with crosslinks into one.

First, concatenate files into one file. Then, merge crosslinks from different
files that are on same position and sum their scores.
"""
# pylint: disable=unused-import
from .. files.bed import merge_bed as run

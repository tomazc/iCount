**************************
Frequently asked questions
**************************

.. contents::
  :local:
  :depth: 1



When running ``peaks`` or other commands I get an "[OSError] 28" message.
===========================================================================
The error number indicates that you are running out of storage ("No space left on device").
Make sure that you have enough free storage and that ``iCount.TMP_ROOT`` (or bash variable ``ICOUNT_TMP_ROOT``) points to a folder with enough free storage.

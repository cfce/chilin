.. ChiLin documentation master file, created by
   sphinx-quickstart on Mon Jun 25 17:22:14 2012.
   You can adapt this file completely to your liking, but it should at least
   rep_control_output
   contain the root `toctree` directive.

ChIP-seq Data quality and Analysis Pipeline
===========================================

Introduction
============
* ChIP-seq experiment has been a mature and wide-spread technique for detecting the transcription factor, histone modification and chromatin regulator distribution on the genome scale. DNase-seq is a high-throughput technoque for finding the markers of open chromatin regions.

* Along with the popularity of the technique and the increasingly huge number of highthroughput datasets, it may be confusing for biologists to get a quick and easy access to understand their biological meaning, and the same important thing is the unbiased judgement of the data quality. So, it's necessary for us to establish a universal and user-friendly ChIP-seq data analysis pipeline for biologists.

* For bioinformaticians, this may be one of the most extensible and flexible ChIP-seq and DNase-seq pipeline implemented with python so far.

* ChiLin has been used in many Cistrome related projects which generated thousands of ChIP-seq and DNase-seq datasets.

**Features**

ChiLin_ provides a more flexible solution for understanding the ChIP-seq analysis workflow. ChiLin_ is intimately part of the Cistrome_ project. ChiLin has been fully tested on Ubuntu, CentOS, and Mac, and validate nearly 9000 datasets on centos 6.0 with slurm_ system, the metrics will be released after the publish of ChiLin paper. Cistrome related quality metrics will be released along with the paper.

.. _ChiLin: http://cistrome.org/chilin
.. _Cistrome: http://cistrome.org
.. _slurm: https://computing.llnl.gov/linux/slurm/

* :download:`QC report Summary Table instructions<../chilin2/modules/summary/instructions.pdf>`
* :download:`demo reports<FileS.tgz>`

.. * :download:`PDF version manual <_build/latex/ChiLin.pdf>`

* :ref:`Instructions to data analysis results<Instructions_results>`

Table of contents
=================
.. toctree::
   :maxdepth: 2

   Installation
   Manual
   API
   Appendix

------

**Getting help**

Any question on installation or runnning is appreciated, please mail
to 1410771@tongji.edu.cn.

**reference**
Manuscript is being submitted.

Others
======
* :ref:`genindex`
* :ref:`search`

.. * :ref:`modindex`

.. ChiLin documentation master file, created by
   sphinx-quickstart on Mon Jun 25 17:22:14 2012.
   You can adapt this file completely to your liking, but it should at least
   rep_control_output
   contain the root `toctree` directive.

The ChiLin --- ChIP-seq Data quality and Analysis Pipeline
===========================================================

Introduction
-----------------

* ChIP-seq experiment has been a mature and wide-spread technique for detecting the TF, histone modification and chromatin factors distribution on the genome scale.
* Along with the popularity of the technique and the increasingly huge number of highthroughput datasets, it may be confusing for biologists to get a quick and easy access to understand their biological meaning, and the same important thing is the unbiased judgement of the data quality. So, it's necessary for us to establish a universal and user-friendly ChIP-seq data analysis pipeline for biologists.
* For bioinformaticians, this is one of the most extensible and flexible ChIP-seq implemented with python so far. It support various genomic data format and integrate high-rated analysis toolbox.

**Features**

ChiLin_ provides a more flexible handle for understanding the ChIP-seq analysis workflow. ChiLin_ is part of the Cistrome_ project.
Fully tested on Ubuntu, CentOS, and Mac. and quality control nearly 9000 datasets on centos 6.0 with _slurm system.

Here is the PDF version of documentation, :download:`PDF version manual <_build/latex/ChiLin.pdf>`.

.. _ChiLin: https://github.com/cfce/chilin
.. _Cistrome: http://cistrome.org

.. _slurm: https://computing.llnl.gov/linux/slurm/

* :download:`QC report Summary Table instructions<../chilin2/modules/summary/instructions.pdf>`
* :ref:`Instructions to data analysis results<Instructions_results>`

Table of contents
-----------------------
.. toctree::
   :maxdepth: 2

   Installation
   Manual
   FAQ
   API
   Appendix


------

**Getting help**

Any question on installation or runnning is appreciated, please mail
to 1410771@tongji.edu.cn

**reference**
Manuscript are being submitted.

Others
-------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

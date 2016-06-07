===============
Installation
===============

The two clones are synchronized between https://github.com/cfce/chilin/. We have packaged all dependent software and species-specific data into **.tar.gz**.

Before installation, make sure that you have gcc, g++, make and java in place, we provide the installation for common system for these dependency.

.. _dependentsoft:

Dependent software list
=======================

.. note::

   python must be 2.7 version for macs2 support.

============================   ==============================      ======================================
Tool                           debian/centos/mac                     Usage for ChiLin
============================   ==============================      ======================================
python dev header              apt-get or yum or port              prerequisites
python setuptools              apt-get or yum or port              prerequisites 
python numpy package           apt-get or yum or port              prerequisites 
cython                         apt-get or yum or port              prerequisites 
`R`_                           apt-get or yum or manually          prerequisites 
java/gcc/g++                   apt-get/yum install/Xcode           prerequisites 
`ghostscript`_                 apt-get or yum or manually          prerequisites 
texlive-latex                  apt-get or yum or manually          prerequisites 
`ImageMagick`_                 apt-get or yum or manually          prerequisites 
`MACS2`_                       pypi                                peak calling
`seqtk`_                       built-in                            packaged into chilin
bx-python_                     built-in                            packaged into chilin
`FastQC`_                      built-in                            packaged into chilin
`BWA`_                         built-in                            packaged into chilin
`samtools`_                    built-in                            packaged into chilin
`bedtools`_                    built-in                            packaged into chilin
bedClip                        built-in                            `UCSC binary`_ (packaged into chilin)
bedGraphToBigWig               built-in                            `UCSC binary`_ (packaged into chilin)
wigCorrelate                   built-in                            `UCSC binary`_ (packaged into chilin)
wigToBigWig                    built-in                            `UCSC binary`_ (packaged into chilin)
mdseqpos                       built-in                             packaged into chilin
============================   ==============================      ======================================


.. _seqtk: https://github.com/lh3/seqtk
.. _mdseqpos: https://bitbucket.org/cistrome/cistrome-applications-harvard
.. _phantompeakqc: http://code.google.com/p/phantompeakqualtools/
.. _FastQC: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
.. _R: http://www.r-project.org/
.. _CEAS site: http://liulab.dfci.harvard.edu/CEAS/download.html
.. _pdflatex: http://www.tug.org/applications/pdftex/
.. _samtools: http://samtools.sourceforge.net/
.. _Bowtie: http://bowtie-bio.sourceforge.net/index.shtml
.. _BWA: https://github.com/lh3/bwa
.. _Bigwiggle: http://genome.ucsc.edu/goldenPath/help/bigWig.html
.. _SRA Toolkit: http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?cmd=show&f=software&m=software&s=software
.. _Cistrome: http://Cistrome.org
.. _MACS2: https://github.com/taoliu/MACS
.. _bedtools: http://code.google.com/p/bedtools/
.. _samtools: http://samtools.sourceforge.net/
.. _Bowtie: http://bowtie-bio.sourceforge.net/index.shtml
.. _R: http://www.r-project.org/
.. _FASTQ: http://en.wikipedia.org/wiki/FASTQ_format
.. _UCSC table browser: http://genome.ucsc.edu/cgi-bin/hgTables?command=start
.. _ENCODE track: https://sites.google.com/site/anshulkundaje/projects/blacklists
.. _UCSC binary: http://hgdownload.cse.ucsc.edu/admin/exe/
.. _bx-python: https://bitbucket.org/james_taylor/bx-python/wiki/Home
.. _ImageMagick: http://www.imagemagick.org/script/download.php
.. _ghostscript: http://www.ghostscript.com/download/gsdnld.html

--------

Ubuntu and debian installation
==============================
* If you are the administrator, use followings.

     .. code-block:: bash

                     sudo apt-get update
                     #python
                     sudo apt-get install python-dev python-numpy python-setuptools cython python-pip
                     #R
                     sudo apt-get install r-base
                     #java
                     sudo apt-get install default-jre
                     sudo apt-get install ghostscript
                     sudo apt-get install imagemagick --fix-missing
                     #Tex
                     sudo apt-get install texlive-latex-base

Then continue to install :ref:`chilin <install_chilin>`.

--------

Centos and Fedora installation
==============================
For centos, use:

     .. code-block:: bash

                     sudo yum install python-devel numpy python-setuptools python-pip 
                     rpm -Uvh http://mirror.chpc.utah.edu/pub/epel/6/x86_64/epel-release-6-8.noarch.rpm
                     sudo yum install tcl tcl-devel tk-devel
                     sudo yum install R
                     sudo yum install ImageMagick
                     sudo easy_install Cython
                     sudo yum install tetex

Then continue to install :ref:`chilin <install_chilin>`.

----------

Mac OS Installation
===================
- Install Xcode_ from MacStore.
- See the instruction of installing java_ for mac.
- Type `xcode-select --install` to install Command Line Tool.
- Install MacPorts_ to help install the dependent modules easier.

.. _Xcode: http://itunes.apple.com/us/app/xcode/id497799835?l=zh&mt=12

.. _java: http://www.java.com/en/download/apple.jsp

.. _MacPorts: http://www.macports.org/

For mac, we suggest using `macports`_, before install macport, user need to have Xcode and Java installed:

     .. code-block:: bash

                     ## download and install macport
                     # open https://distfiles.macports.org/MacPorts/ and download the right version
                     sudo port install py27-setuptools py27-pip py27-nose py27-cython py27-numpy @1.8.1  ## or use EPD to replace this
                     ## Install R manually from http://cran.cnr.berkeley.edu/bin/macosx/
                     ## For mac latex, install separately, download and click to install it
                     wget -c http://mirror.ctan.org/systems/mac/mactex/MacTeX.pkg

.. _macports: http://www.macports.org/

Install R_, MacTex, ImageMagick_ and ghostscript_ manually.

Then continue to install :ref:`chilin <install_chilin>`.

--------

.. _install_chilin:

Install ChiLin
==============

Test and install pipeline software
----------------------------------

- Type "which gcc g++ java make gs convert pdflatex R cython" to check the installation.

After solving the :ref:`dependent prerequisites<dependentsoft>`, install chilin as followings,

     .. code-block:: bash

                     git clone http://github.com/cfce/chilin/
                     cd chilin
                     python setup.py install -f


Then, check your installation::

                source chilin_env/bin/activate
                # check ChiLin dependent software and data
                python setup.py -l

If any software can not be installed, look into their official documentation. Most of time, see dependentsoft_ to check whether all prerequisites are installed or not, usually it's the problem of numpy, cython or gcc compiler problem, or R package `seqLogo` problem. Take a look at all software in the *software* directory to see what's going on, and try apt-get, yum, port and pypi to fix the issue.

Lastly, user may need to check the installation of mdseqpos dependency of R `seqLogo` package, open R console and install dependent R packages:

              .. code-block:: bash

                              R -e "source('http://bioconductor.org/biocLite.R');biocLite('seqLogo');library(seqLogo)"


Remember to source your `python virtual environment` "source ${ChiLin_ROOT}/chilin_env/bin/activate" everytime or put them into your ${HOME}/.bashrc or ${HOME}/.bash_profile.

.. note::

   After installation, the config file is auto-generated and set the species specific data directory default to `db` under the code root directory.


Download dependent data for hg38_, hg19_, mm9_, or mm10_
---------------------------------------------------------

under the ChiLin source code root directory,

  .. code-block:: bash

                  # download from our cistrome server
                  mkdir -p db

                  # change directory to db
                  cd db

                  # download the one you need, this would be over 10 GB, make sure your internet access is over 100k/s, or it's too slow..
                  # human
                  wget -c http://cistrome.org/chilin/_downloads/hg19.tgz
                  wget -c http://cistrome.org/chilin/_downloads/hg19.tgz.md5 ## check md5
                  #wget -c http://cistrome.org/chilin/_downloads/hg38.tgz
                  #wget -c http://cistrome.org/chilin/_downloads/hg38.tgz.md5

                  # mouse
                  #wget -c http://cistrome.org/chilin/_downloads/mm9.tgz
                  #wget -c http://cistrome.org/chilin/_downloads/mm9.tgz.md5
                  #wget -c http://cistrome.org/chilin/_downloads/mm10.tgz
                  #wget -c http://cistrome.org/chilin/_downloads/mm10.tgz.md5

                  # check the md5sum for completeness of hg19
                  md5sum -c hg19.tgz
                  tar xvfz hg19.tgz

                  # download mycoplasma that you are afraid of contaminating your samples
                  wget -c http://cistrome.org/chilin/_downloads/mycoplasma.tgz
                  wget -c http://cistrome.org/chilin/_downloads/mycoplasma.tgz.md5
                  md5sum -c mycoplasma.tgz.md5
                  tar xvfz mycoplasma.tgz

                  # change back
                  cd ..

                  # check your data and software installation, if download is ok
                  python setup.py -l


.. _hg38:  http://cistrome.org/chilin/_downloads/hg38.tgz
.. _mm10:  http://cistrome.org/chilin/_downloads/mm10.tgz
.. _hg19:  http://cistrome.org/chilin/_downloads/hg19.tgz
.. _mm9:   http://cistrome.org/chilin/_downloads/mm9.tgz

If you

* want to know more about dependent data
* want to prepare new version of reference data by yourself
* you have species assembly not in our list

see details about the :ref:`dependent data<dependentdata>`.

add species support in chilin
-----------------------------

After these preparation of software and reference data, if you are using our prepared hg38_, hg19_, mm10_, mm9_ dependent data, you can skip this part because setup.py already sets `chilin.conf.filled` for you.
If you have your own reference data, open your favorite text editor, appending section in `chilin.conf.filled` file add in species support like this.

fill in the section with your own data absolute path, then append filled following section

.. literalinclude:: hg19.conf
   :language: ini
   :linenos:

after conf, see details about this in :envvar:`[species] <[species]>`.

.. literalinclude:: ../chilin.conf


Test installation
-----------------
Test installation with demo data,

     .. code-block:: bash

                     # non cluster server
                     cd demo
                     bash foxa1

                     # if you are using slurm sytem
                     cd demo
                     # submit cluster script `foxa1`
                     sbatch foxa1

Check demo data results,

     .. code-block:: bash

                     du -h local/local.pdf ## quality report
                     # mac
                     open local/local.pdf
                     # linux
                     nautilus local/


For more options, see :ref:`Manual <Manual>`.

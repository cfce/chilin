===============
Installation
===============
:abbr:`ChiLin (ChIP-seq pipeline)`

ChiLin_ Theoretically, we support all species listed on UCSC, which have masked, non-masked genomics fasta, chromosome information file, and standard RefSeq files.
The species we have tested includes: hg19, hg38, mm9, mm10. We have successfully deployed on Ubuntu, CentOS, Mac, do not support for Windows.

.. _ChiLin: https://github.com/cfce/chilin

Dependent tools Prerequisites
==============================
Here is a list of prerequisites, we use python to develop chilin. All related python tool and packages should be installed under the same version of python (fully tested under 2.6/2.7), 
If you already have these software, skip to `virtualenv_installation`_.

============================   ==============================  ==================================
Tool Name                      Download source                 Usage
============================   ==============================  ==================================
git (optional)                 apt-get or yum install          for sync latest tool
mercurial (optional)           apt-get or yum install          for sync latest tool
python dev header              apt-get or yum install          for installing numpy
python setuptools              apt-get or yum install          python package distribution package
python numpy package           apt-get or yum install          for installing macs and bx-python
bx-python_                     pip                             parsing bigwig
cython                         apt-get or yum install          for installing macs and bx-python
R                              apt-get or yum install          for plotting
java/gcc/g++                   apt-get/yum install/Xcode       for FastQC/make/samtools etc.
ghostscript                    apt-get or yum install          for motif plotting
texlive-latex                  apt-get or yum install          for rendering latex template
ImageMagick                    apt-get or yum install          convert image
============================   ==============================  ==================================

* If you are the administrator, use followings.

For ubuntu or debian system, install as follows:

     .. code-block:: bash

		     sudo apt-get update
		     sudo apt-get install git
		     sudo apt-get install mercurial
		     sudo apt-get install python-dev python-numpy python-setuptools cython python-pip ## or install redhat EPD python
		     sudo apt-get install r-base
		     sudo apt-get install default-jre
		     sudo apt-get install ghostscript
                     sudo apt-get install imagemagick --fix-missing
		     sudo apt-get install texlive-latex-base
		     sudo pip install bx-python


For centos, use:

     .. code-block:: bash

		     sudo yum install git
		     sudo yum install mercurial
		     rpm -Uvh http://mirror.chpc.utah.edu/pub/epel/6/x86_64/epel-release-6-8.noarch.rpm
		     sudo yum install tcl tcl-devel tk-devel
		     sudo yum install R
		     sudo yum install python-devel numpy python-setuptools python-pip ## or use EPD package
                     sudo yum install ImageMagick
		     sudo easy_install Cython
		     sudo easy_install bx-python
		     sudo yum install tetex


For mac, we suggest using `macports`_, before install macport, user need to have Xcode and Java installed:

     .. code-block:: bash

		     sudo port install git
		     sudo port install mercurial
		     sudo port install py27-setuptools py27-pip py27-nose py27-cython py27-numpy @1.8.1  ## or use EPD to replace this
		     ## Install R manually from http://cran.cnr.berkeley.edu/bin/macosx/
		     ## For mac latex, install separately, download and click to install it
		     wget -c http://mirror.ctan.org/systems/mac/mactex/MacTeX.pkg


.. _macports: http://www.macports.org/

.. _virtualenv_installation:

* After solving the dependency need above, hopefully you can using the following command to install all dependent software in the virtualenv::

     git clone https://github.com/cfce/chilin
     cd chilin
     bash install_linux.sh or bash install_mac.sh

then see `mdseqpos`_ installation and download `dependentdata`_, currently mdseqpos cannot be installed on mac, this step can be :ref:`"skip" <skip>` with *--skip 12*


* Or follow the instruction for step-by-step installation.

Overall:
1. install `dependent software`_
2. download reference `dependentdata`_
3. fill in chilin.conf or chilin.conf.filled(this has high priority, do not use this two together), `fill_conf`_
4. `python setup.py install`

.. note::

   each time, users add a new species section or modify a species reference section in the chilin.conf, users need to reinstall chilin

.. _dependent software:

Dependent softwares 
========================
Following installation are tested under the shell under Ubuntu, CentOS, and Mac.

Here is the list of published software and `UCSC binary`_ we have fully tested with ChiLin:

============================      ====================================  ==============
Tool Name		          usage                                 version tested
============================      ====================================  ==============
`seqtk`_                            sub-sample FASTQ files              1.0
`FastQC`_                           sequence quality and gc contents    v0.10.1
`BWA`_                              mapping Fastq                       0.7.7
`samtools`_                         convert SAM to BAM                  0.1.19
`bedtools`_		            operate bed files                   v2.17.0
`MACS2`_                            peak calling                        2.0.10.2014XXXX
bedClip                             trim outlier of coordinates 
bedGraphToBigWig                    convert bedgraph to bigwig          v 4
wigCorrelate                        calculate reads count correlation
wigToBigWig                         convert wiggle to bigwig            v 4
`MDSeqPos`                          motif scan and SeqPos               v 2.01
============================      ====================================  ==============


.. _virtualenv: https://raw.githubusercontent.com/pypa/virtualenv/1.9.X/virtualenv.py

0 step is to setup `virtualenv`_, activate virtualenv each time before using `chilin`, use `macports` python2.7 as python virtualenv version:

     .. code-block:: bash

		     wget -c --no-check-certificate https://raw.githubusercontent.com/pypa/virtualenv/1.9.X/virtualenv.py
		     python virtualenv.py -p/opt/local/Library/Frameworks/Python.framework/Versions/2.7/Resources/Python.app/Contents/MacOS/Python --system-site-packages --distribute chilin_env
		     source chilin_env/bin/activate

First, seqtk for sampling fastq/fastq.gz, download from github:

     .. code-block:: bash
	
		     git clone https://github.com/lh3/seqtk
		     cd seqtk && make && chmod 755 seqtk && cp seqtk ../chilin_env/bin && cd ..

For UCSC utility, Linux `x86_64` users use the followings:

     .. code-block:: bash

		     wget -c http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig
		     wget -c http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bedClip
		     wget -c http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/wigCorrelate
		     wget -c http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/wigToBigWig
		     chmod 755 bedGraphToBigWig bedClip wigCorrelate wigToBigWig && cp bedGraphToBigWig bedClip wigCorrelate wigToBigWig chilin_env/bin


For mac `x86_64` users use the followings:

     .. code-block:: bash

		     wget -c http://hgdownload.cse.ucsc.edu/admin/exe/macOSX.x86_64/bedGraphToBigWig
		     wget -c http://hgdownload.cse.ucsc.edu/admin/exe/macOSX.x86_64/bedClip
		     wget -c http://hgdownload.cse.ucsc.edu/admin/exe/macOSX.x86_64/wigCorrelate
		     wget -c http://hgdownload.cse.ucsc.edu/admin/exe/macOSX.x86_64/wigToBigWig
		     cp bedGraphToBigWig bedClip wigCorrelate wigToBigWig chilin_env/bin
		     chmod 755 bedGraphToBigWig bedClip wigCorrelate wigToBigWig && cp bedGraphToBigWig bedClip wigCorrelate wigToBigWig chilin_env/bin

Alternatively, users can download and compile them yourselves from `UCSC utility`_:

.. _UCSC utility: http://genome-source.cse.ucsc.edu/gitweb/?p=kent.git;a=blob;f=src/userApps/README

Then for `FastQC`_, download from the link, and put the binary files into your virtualenv:

     .. code-block:: bash

		     wget -c http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.10.1.zip
		     unzip fastqc_v0.10.1.zip
		     chmod -R 755 FastQC
		     cp -r FastQC/* chilin_env/bin


Then for bwa, samtools and bedtools, download from the link and compile to put the binary into your virtualenv:

     .. code-block:: bash

		     wget -c --no-check-certificate https://github.com/lh3/bwa/archive/0.7.7.tar.gz
		     tar xvfz 0.7.7.tar.gz
		     cd bwa*
		     make && cp bwa ../chilin_env/bin && cd ..

		     wget -c --no-check-certificate http://sourceforge.net/projects/samtools/files/samtools/0.1.19/samtools-0.1.19.tar.bz2/download -O samtools.tar.bz2
		     tar xvfj samtools.tar.bz2
		     cd samtools* && make && cp samtools ../chilin_env/bin && cd ..

		     wget -c --no-check-certificate http://github.com/arq5x/bedtools/archive/v2.17.0.tar.gz
		     tar xvfz v2.17.0.tar.gz
		     cd bedtools* && make && cp bin/* ../chilin_env/bin && cd ..

Then peaks calling, we use high-rated macs2 for peaks calling, for version `2.0.10.2014XXXX`, python version needs to be 2.7 for installing macs2:

     .. code-block:: bash

		     wget -c --no-check-certificate https://github.com/taoliu/MACS/archive/de419af20b6b441f0c0b62a3c18a262228409c8a.zip -O macs_2.0.10.2014XXXX.zip 
		     unzip  macs_2.0.10.2014XXXX.zip 
		     cd MACS*
		     python setup_w_cython.py build_ext
		     python setup.py install
       

Then, Installation of `mdseqpos`_ needs R `seqLogo` package, open R console and install dependent R packages:

     .. code-block:: python

		     source("http://bioconductor.org/biocLite.R")
		     biocLite("seqLogo")
		     library(seqLogo)

Download `mdseqpos`_ dependent masked and non-masked sequence data. Refer to `genome`_ for filling settings.py.

at last, install `mdseqpos`_:

     .. code-block:: bash

		     wget -c --no-check-certificate https://bitbucket.org/cistrome/cistrome-applications-harvard/get/cfec0147b6f6.zip
		     unzip cfec0147b6f6.zip
		     cd cistrome*/mdseqpos
		     cp lib/settings.py.example lib/settings.py
		     python setup.py install
      
Take hg19 as example. Create /data/motif/assembly/hg19/raw and /data/motif/assembly/hg19/masked directory. Download raw and masked genome data. Modify mdseqpos/lib/settings.py.example in lines 55-64 according to where you put the sequence data. 

.. _mdseqpos document: https://bitbucket.org/cistrome/cistrome-applications-harvard/src/270e986a26a0ab3b143dbf9720d990f3654f677f/mdseqpos/?at=default

List of tested species sequence data, 

.. _genome:

============================   ==============================  ==================================
Genome version                  Raw genome sequence             Masked genome sequence
============================   ==============================  ==================================
hg19                           hg19_raw_		       hg19_mask_
hg38			       hg38_raw_		       hg38_mask_
mm9			       mm9_raw_  		       mm9_mask_
mm10			       mm10_raw_		       mm10_mask_
============================   ==============================  ==================================

uncompress the above *raw* and *masked* sequence into directory, rename *masked* *chrxxx.fa* to *chrxxx.fa.masked*, then fill in `settings.py` then install mdseqpos, See `mdseqpos document`_ if necessary.

.. _hg19_raw: http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/chromFa.tar.gz
.. _hg19_mask: http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/chromFaMasked.tar.gz
.. _hg38_raw: http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromFa.tar.gz
.. _hg38_mask: http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromFaMasked.tar.gz
.. _mm9_raw: http://hgdownload.cse.ucsc.edu/goldenpath/mm9/bigZips/chromFa.tar.gz
.. _mm9_mask: http://hgdownload.cse.ucsc.edu/goldenpath/mm9/bigZips/chromFaMasked.tar.gz
.. _mm10_raw: http://hgdownload.cse.ucsc.edu/goldenpath/mm10/bigZips/chromFa.tar.gz	       
.. _mm10_mask: http://hgdownload.cse.ucsc.edu/goldenpath/mm10/bigZips/chromFaMasked.tar.gz

	       
For other genomes, here is the major list we tested:
	
* Now fill in the chilin.conf in chilin :envvar:`[tool] section <[tool]>`::

    [tool]
    mdseqpos = path/MDSeqPos.py
    macs2 = path/macs2


.. _dependentdata:


Dependent data
==============

We have packaged reference data for hg19/hg38/mm9/mm10 on our server.

There are two categories of data.

First is large disk usage data:

============================   =====================  =========================
Data Name                       Used by                Data Source         
============================   =====================  =========================
:envvar:`genome_index`          bwa/bowtie/star        raw fasta indexed files
:envvar:`conservation`          conservation_plot.py   wiggle files
============================   =====================  =========================

Second is small pieces of reference files:

.. _reference_small:

.. code-block:: bash

		## hg19
		wget -c http://cistrome.org/~qqin/chilin/materials/hg19_chilin.tar.gz
		## hg38
		wget -c http://cistrome.org/~qqin/chilin/materials/hg38_chilin.tar.gz
		## mm9
		wget -c http://cistrome.org/~qqin/chilin/materials/mm9_chilin.tar.gz
		## mm10
		wget -c http://cistrome.org/~qqin/chilin/materials/mm10_chilin.tar.gz
		## contamination index
		wget -c http://cistrome.org/~qqin/chilin/materials/mycoplasma_bwa_index.tar.gz
                ## hg38 index
   		wget -c http://cistrome.org/~qqin/chilin/materials/hg38_bwa_index.tar.gz
                ## mm10 index
   		wget -c http://cistrome.org/~qqin/chilin/materials/mm10_bwa_index.tar.gz
                ## mm9 index
   		wget -c http://cistrome.org/~qqin/chilin/materials/mm9_bwa_index.tar.gz
                ## hg19 index
   		wget -c http://cistrome.org/~qqin/chilin/materials/hg19_bwa_index.tar.gz
                ## hg19 placetalmammals Phastcons conservation tracks
   		wget -c http://cistrome.org/~qqin/chilin/materials/hg19_placental_conservation.tar.gz
                ## mm9 placetalmammals Phastcons conservation tracks
                wget -c http://cistrome.org/~qqin/chilin/materials/mm9_placental_conservation.tar.gz

============================   =====================  ======================================
Data Name                       Used by                Data Source         
============================   =====================  ======================================
:envvar:`chrom_len`             samtools              `UCSC table browser`_  
:envvar:`chrom_bed`             bedtools              `UCSC table browser`_
:envvar:`regpotential`          RegPotential          `UCSC table browser`_
:envvar:`dhs`                   bedtools               Union DHS regions home processed
:envvar:`velcro`                bedtools	      `ENCODE track`_ blacklist regions
:envvar:`geneTable`             bedAnnotate           `UCSC table browser`_
:envvar:`ceas_exon`             bedtools	      `UCSC table browser`_ merged exons
:envvar:`ceas_intron`           bedtools	      `UCSC table browser`_ merged introns
:envvar:`ceas_intergenic`       bedtools	      `UCSC table browser`_ merged intergenic
:envvar:`ceas_promotor`         bedtools	      `UCSC table browser`_ merged promoter
:envvar:`[contamination]`	bwa                   `Mycoplasma`_
============================   =====================  ======================================


* Followings is how we generate these reference files, if you have any species other than hg19/hg38/mm9/mm10, you can find the reference files with the similar ways.

Mycoplasma genome
---------------------

It seems that Mycoplasma_ contamination would be a major source of contamination, so we recommended downloading the Mycoplasma fasta for indexing, data is in the link of the `mycoplasma genome`_.

.. _Mycoplasma: http://www.biodatamining.org/content/7/1/3/abstract?utm_campaign=22_05_14_BioDataMining_ArticleMailing_EBM_PA_REG_BMCUP&utm_content=8772920153&utm_medium=BMCemail&utm_source=Emailvision
.. _mycoplasma genome: http://mycoplasma.genome.uab.edu/genomes.asp

Then index with `bwa index -a is mycoplasma.fasta`.

BWA Index 
--------------------

download raw `genome`_ sequence data, and `tar xvfz` them and `cat *fa > genome.fa`. Use the following to index them:

.. code-block:: bash

		bwa index -a bwtsw genome.fasta

Then fill in species specific reference section `species_ref`_.


UCSC table browser
---------------------
* To get merged exons, introns, promoters and intergenic regions, open `UCSC table browser`_

Use Browser step by step
1. Go to the UCSC table browser.
2. Select desired species and assembly, such as hg19
3. Select group: Genes and Gene Prediction Tracks
4. Select track: Refseq
5. Select table: genes and prediction track
6. Select region: genome
7. Select output format: BED - browser extensible data
8. Enter output file: UCSC_exon/gene.bed
9. Hit the 'get output' button
10. A second page of options relating to the BED file will appear.
11. Under 'create one BED record per:'. Select 'Whole Genome/Exons Plus'
12. Hit the 'get BED' option.

If you have mysql configured properly, auxiliary scripts modified from biostar may be useful, for debian and for hg19, or use `reference_small`_ data above::

  sudo apt-get install mysql-client-core-5.5  

  1.bash chilin2/modules/ceas/chromInfo.sh hg19  ## get chrom_len and chrom_bed
  2.bash chilin2/modules/ceas/meta_info.sh gene.bed exon.bed 2000 chromInfo_file(from 1) ## get ceas_exon, ceas_intron, ceas_intergenic, ceas_promotor
  3.bash chilin2/modules/ceas/ceas_reg.sh hg19  ## get regpotential reference file


Conservation score
----------------------

* (Optional) get Phaston conservation, for most common species version, such as hg19, get from `hg19_conserv`_ or `mm9_conserv`_, and use wigToBigWig to convert them into bigwig, we provide hg19/mm9 conservation score on our server, for other species, just left the chilin.conf conservation section blank. Take hg19 as an example::

                  wget -r -np -nd --accept=gz http://hgdownload-test.cse.ucsc.edu/goldenPath/hg19/phastCons46way/placentalMammals/
		  for c in chr*wig*gz
		  do
          	  bw=${c%phastCons46way.placental.wigFix.gz}bw
		  echo $bw
		  gunzip -c $c | wigToBigWig stdin chrom_len $bw  ## chrom_len is where you put your reference chromosome information file
		  done


.. _hg19_conserv: http://hgdownload-test.cse.ucsc.edu/goldenPath/hg19/phastCons46way/placentalMammals/
.. _mm9_conserv: http://hgdownload-test.cse.ucsc.edu/goldenPath/mm9/phastCons30way/placental/

.. _UCSC table browser: http://genome.ucsc.edu/cgi-bin/hgTables?command=start

Test
------------
After downloading, users can use our pre-prepared small size datasets for testing:

  .. code-block:: bash

		  wget -c cistrome.org/~qqin/test_data/foxa1_t1.fastq.gz
		  wget -c cistrome.org/~qqin/test_data/foxa1_c1.fastq.gz
                  wget -c http://cistrome.org/~qqin/test_data/chilin_test.tar.gz ## reference test data


.. _species_ref:

Generating Species References
-------------------------------
After these preparation of software and reference data, you need to fill in chilin.conf,

.. _fill_conf:

.. literalinclude:: ../chilin.conf

Take hg19 as an example,

.. literalinclude:: hg19.conf

:language: ini
:linenos:

You can append as many species as you want.

The last step of installation is running `python setup.py install`.

Now try :ref:`"quick start" <quick-start>` to run ChiLin.


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

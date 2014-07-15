===============
Installation
===============
:abbr:`ChiLin (ChIP-seq pipeline)`
   
ChiLin_ Theoretically, we support all species listed on UCSC, which have masked, non-masked genomics fasta file, and standard RefSeq files.
The species we have tested includes: hg19, hg38, mm9, mm10. We have successfully deployed on Ubuntu, CentOS, Mac.

.. _ChiLin: https://bitbucket.org/Alvin_Qin/chilin/overview

Here is the instruction for step-by-step installation.

Prerequisites
-----------------

Here is a list of prerequisites:

============================   =======================  ==================================
Tool Name                      Download source           Usage
============================   =======================  ==================================
git (optional)                 apt-get or yum install    for sync latest tool
mercurial (optional)           apt-get or yum install    for sync latest tool
python dev header              apt-get or yum install    for installing numpy
python numpy package           apt-get or yum install    for installing macs and bx-python
cython                         apt-get or yum install    for installing macs and bx-python
R                              apt-get or yum install    for plotting
java                           apt-get or yum install    for FastQC
ghostscript                    apt-get or yum install    for motif plotting
texlive-latex                  apt-get or yum install    for rendering latex template
============================   =======================  ==================================

Most of the softwares are under source control system, you'd better install, on ubuntu or debian system, install as follows::

        sudo apt-get update
	sudo apt-get install git
	sudo apt-get install mercurial

Alternatively, you can download this source control software by hand.  

Install python develop header files and numpy::
       
        sudo apt-get install python-dev python-numpy cython

Install R for conservation plot and Motif analysis::

	sudo apt-get install r-base

Install java for FastQC::

        sudo apt-get install default-jre


Install python package management system files and numpy::


	sudo apt-get install python-pip


..
   Install dependency C++ package boost for NSC/RSC/Qtag calculation::

	   sudo apt-get install libboost-all-dev

.. _boost:	


Install prerequisites for motif::

        sudo apt-get install ghostscript

Install *pdflatex* for rending pdf::

        sudo apt-get install texlive-latex-base


Dependent softwares
--------------------
Following installation are tested under the bash shell and zshell under Ubuntu, CentOS, and Mac.

Here is the list of the software we integrate in ChiLin:

============================      ====================================
Tool Name		          usage
============================      ====================================
`seqtk`_                            sub-sample FASTQ files
`FastQC`_                           sequence quality and gc contents
`MACS2`_                            peak calling
`bedtools`_		            operate bed files
`samtools`_                         convert SAM to BAM
`BWA`_                              mapping Fastq
`Bowtie`_	                    library contamination
============================      ====================================

It seems that Mycoplasma_ contamination would be a major source of contamination, so we recommended downloading the Mycoplasma fasta for indexing, data is in the link of [[http://mycoplasma.genome.uab.edu/genomes.asp][mycoplasma]] ::

.. _Mycoplasma: http://www.biodatamining.org/content/7/1/3/abstract?utm_campaign=22_05_14_BioDataMining_ArticleMailing_EBM_PA_REG_BMCUP&utm_content=8772920153&utm_medium=BMCemail&utm_source=Emailvision

first, we use bwa for mapping and samtools for post-mapping filtering::

	sudo apt-get install bwa # or 	sudo apt-get install bowtie
	sudo apt-get install samtools
  

Second, install *bedtools* for bed and bam comparison::

        sudo apt-get install bedtools

Then, we can install bx-python_  for plotting conservation distribution::

        sudo pip install bx-python


0 step is to setup a PATH bin for your local installation::

                    mkdir ~/bin
                    echo -e "export PATH=${PATH}:~/bin" >> ~/.bashrc

First, seqtk for sampling fastq/fastq.gz, download from github:

     .. code-block:: bash
	
		     git clone https://github.com/lh3/seqtk
		     cd seqtk
		     make
	             cp seqtk ~/bin	

Then for `FastQC`_, download from the link, and put the binary files into your PATH::

  
     wget -c http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.10.1.zip
     unzip fastqc_v0.10.1.zip
     chmod -R 755 FastQC
     cp -r FastQC/* ~/bin/
  

..
   Then for NSC/RSC/Qtag, we use phantompeakqc packages derived from spp peaks caller, `phantompeakqc`_ , this needs `boost`_ for compiling:

   - first, spp dependency, enter R console::

	install.packages("caTools")


   - second, install spp.

       .. code-block:: bash

			wget -c http://phantompeakqualtools.googlecode.com/files/ccQualityControl.v.1.1.tar.gz
			tar xvfz ccQualityControl.v.1.1.tar.gz
			cd phantompeakqualtools/
			R CMD INSTALL spp_1.10.1.tar.gz

Then peaks calling, we use high-rated macs2 for peaks calling::

    sudo pip install MACS2


Then, install `mdseqpos`_ and required dependency, refer to `genome`_ for filling settings.py::

    hg clone https://bitbucket.org/cistrome/cistrome-applications-harvard
    cd cistrome-applications-harvard/mdseqpos
    ## refer the following for other steps

	
- fill in the :envvar:`[tool] section <[tool]>`::
  
    [tool]
    mdseqpos = path/MDSeqPos.py
    macs2 = path/macs2
    bedannotate = path/bedAnnotate.py
  
  
Open R console and install dependent R packages::

   source("http://bioconductor.org/biocLite.R")
   biocLite("seqLogo")
   library(seqLogo)


Installation of MDSeqpos needs dependent masked and non-masked sequence data.
  
ChiLin involve `UCSC binary`_ for converting files, download them and put them into PATH:

============================   ==================================
Tool Name                      from        
============================   ==================================
bedClip                        trim outlier of coordinates
bedGraphToBigWig	       convert bedgraph to bigwig
wigCorrelate                   calculate reads count correlation
wigToBigWiggle                 convert wiggle to bigwig
============================   ==================================

For linux x86_64, use the followings::

   wget -c http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig
   wget -c http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bedClip
   wget -c http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/wigCorrelate
   wget -c http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/wigToBigWiggle

  
.. _dependentdata:

Dependent data
==============


Take hg19 as example, we support all species with following annotations.

* Before preparing dependent data, use a uniform data::

    mkdir  ~/data


.. _genome:

* First, preparing raw genome fasta for getting genome index for bwa and bowtie, and masked and raw genome sequence for MDSeqPos sequence database::

    use bowtie source codes contained e_coli for testing your library purity

    > wget -c ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz
    > wget -c ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFaMasked.tar.gz
    > tar xvfz chromFaMasked.tar.gz
    > tar xvfz chromFa.tar.gz
    # cat all base chromosome for indexing genome into hg19.fasta
    > bwa index -p hg19 -a bwtsw hg19.fasta

..
   * (Optional), Motif analysis genome masked and raw genome, such as hg19,

       .. code-block:: bash


	   wget -c http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz

       .. code-block:: bash


	   wget -c http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFaMasked.tar.gz

.. note::

   Our pipeline starts from fastq or fastq.gz file, currently, we do not support bam, bed file, and ab solid fastq file.
  
..
   * if platform is ab solid colorspace sequence, latest bwa may not support it, used version before 0.6.0 instead::

       > bwa index -c -p hg19_color -a bwtsw hg19.fasta

..
   * if platform is ab solid colorspace sequence, latest bwa may not support it, used bowtie instead::
       
       > bwa index -c -p hg19_color -a bwtsw hg19.fasta



* Second, get chromosome information and chromosome bed files, we have a auxiliary script, user needs to install mysql ::

    sudo apt-get install mysql-client-core-5.5
    bash chilin2/modules/ceas/chromInfo.sh hg19

* If you do not have mysql you can download all chromosome information files from the link::
  
   > wget -c from daisy

* To get merged exons, introns, promoters and intergenic regions, open `UCSC table browser`_::

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

    after that, run:
    >bash chilin2/modules/ceas/meta_info.sh gene.bed exon.bed 2000 chromInfo


* get bed files for computing regulatory potentials::

    > bash chilin2/modules/ceas/mysql_reg.sh hg19

* (Optional) get Phaston conservation, for most common species version, such as hg19::

   > wget -c http://hgdownload-test.cse.ucsc.edu/goldenPath/hg19/phastCons46way/placentalMammals/

   > wigToBigWiggle to convert them into bigwiggle, download only the normal chromosome wig.gz files and convert them.

The ChiLin package includes all the build-in data for hg19 and mm9. For other species, you may need to download these data from data source or custom it yourself.


============================   ============  =====================  =========  
Data Name                       Used by       Data Source           Format     
============================   ============  =====================  =========  
Chromesome length              samtools      `UCSC table browser`_  2-column   
DHS region                     bedtools      Custom                 BED
Velcro region                  bedtools	     `ENCODE track`_        BED
Motif database                 MDSeqPos      `mdseqpos`_            xml
FastQC result database         QCreport      Custom                 bed
Data summary database          QCreport      Custom                 bed
============================   ============  =====================  =========

Test
------------
After downloading, users can use our pre-prepared small size datasets for testing::

   > wget -c cistrome.org/~qqin/test_data/foxa1_t1.fastq
   > wget -c cistrome.org/~qqin/test_data/foxa1_c1.fastq

ChiLin Installation
---------------------
To install the tool, you will first have to cp chilin.conf to chilin.conf.filled
and modify chilin.conf to define site-wide defaults for your system.

In chilin.conf, there are three main sections-- tools, references, and params.

In the tools section, you will have to define the absolute paths to the various
sub-tools or modules in the analysis pipeline.

In the references section, you will have to define the absolute paths to the
reference files (e.g. assembly, DHS file, gene regions, etc.) for the species
you wish to support.  Please read 'Generating Species References' below for
information on how to generate files for your species of interest.

In the params section, you will have the ability to fine-tune and define
site-wide defaults for some of the sub-tools/modules.  NOTE: the defaults
should work for almost all installations--only make modifications only if you
are sure!

After tailoring chilin.conf to suit your system's needs, simple type:
python setup.py install

NOTE: if you are installing system-wide, you may want to add 'sudo' in front
of that command

Generating Species References
-------------------------------
Take hg19 as an example,

.. literalinclude:: hg19.conf
   :language: ini
   :linenos:


use the standard python package installation method::

  hg clone https://bitbucket.org/Alvin_Qin/chilin
  cd chilin
  ## support python >= 2.6
  sudo python setup.py install


After these installation, you could try :ref:`"quick start" <quick-start>` to run ChiLin.


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

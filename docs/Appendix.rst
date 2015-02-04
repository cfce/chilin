========================
Appendix: Dependent data
========================

.. _dependentdata:

Get dependent data
------------------

ChiLin support all species listed on UCSC website, which includes dependent data as we list in :envvar:`[species] <[species]>`:

* (Must) genome index for your species, we recommended bwa index.

* (Must) chromosome length information text file.

* (Must) standard *RefSeq* files.

* (Optionally)PhastCons conservation bigwiggle files.

* (Optionally) genome directory containing chromosome separated sequence fasta files

* (Optionally) Union DHS and blacklist regions

We have packaged all dependent data for hg19_, hg38_, mm9_, mm10_.

.. _hg38:  http://cistrome.org/chilin/_downloads/hg38.tgz
.. _mm10:  http://cistrome.org/chilin/_downloads/mm10.tgz
.. _hg19:  http://cistrome.org/chilin/_downloads/hg19.tgz
.. _mm9:  http://cistrome.org/chilin/_downloads/mm9.tgz

.. _ChiLin: https://cistrome.org/chilin

Data details
---------------------------------

* First is large disk usage data:

============================   =====================  =========================
Data Name                       Used by                Data Source
============================   =====================  =========================
:envvar:`genome_index`          bwa/bowtie/star        raw fasta indexed files
:envvar:`genome_dir`            bwa/bowtie/star        genome fasta files
:envvar:`conservation`          conservation_plot.py   wiggle files
============================   =====================  =========================

.. _genome:

============================   ==============================  ==================================
Genome version                  Raw genome sequence             Masked genome sequence
============================   ==============================  ==================================
hg19                           hg19_raw_                       hg19_mask_
hg38                           hg38_raw_                       hg38_mask_
mm9                            mm9_raw_                        mm9_mask_
mm10                           mm10_raw_                       mm10_mask_
============================   ==============================  ==================================

.. _hg19_raw: http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/chromFa.tar.gz
.. _hg19_mask: http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/chromFaMasked.tar.gz
.. _hg38_raw: http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromFa.tar.gz
.. _hg38_mask: http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromFaMasked.tar.gz
.. _mm9_raw: http://hgdownload.cse.ucsc.edu/goldenpath/mm9/bigZips/chromFa.tar.gz
.. _mm9_mask: http://hgdownload.cse.ucsc.edu/goldenpath/mm9/bigZips/chromFaMasked.tar.gz
.. _mm10_raw: http://hgdownload.cse.ucsc.edu/goldenpath/mm10/bigZips/chromFa.tar.gz
.. _mm10_mask: http://hgdownload.cse.ucsc.edu/goldenpath/mm10/bigZips/chromFaMasked.tar.gz

* Second is small pieces of reference files:

============================   =====================  ==========================================
Data Name                       Used by                Data Source
============================   =====================  ==========================================
:envvar:`chrom_len`             samtools              `UCSC table browser`_
:envvar:`dhs`                   bedtools               Union DHS regions from Cistrome DB
:envvar:`velcro`                bedtools               blacklist regions
:envvar:`geneTable`             bedAnnotate           `UCSC table browser`_
:envvar:`[contamination]`	bwa                   `Mycoplasma`_ genome index(set by --mapper)
============================   =====================  ==========================================

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

UCSC table browser
---------------------
Use Browser step by step

* To get refseq files, open `UCSC table browser`_
* Go to the UCSC table browser.
* Select desired species and assembly, such as hg19
* Select group: Genes and Gene Prediction Tracks
* Select track: RefSeq Genes
* Select table: refGene
* Select region: genome
* Select output format: all fields from selected table
* Enter output file: *species*.refgene
* Hit the 'get output' button
* d*ownload and remove the header line with command,

.. code-block:: bash

                sed 1d species.refgene > sp.refgene

Conservation score
----------------------

* (Optional) get Phaston conservation, for most common species version, hg19_conserv_, hg38_conserv_, mm10_conserv_, mm9_conserv_ and use wigToBigWig to convert them into bigwig, we provide hg19/mm9 conservation score on our server, for other species, just left the chilin.conf conservation section blank. Take hg19 as an example:

.. code-block:: bash

                  wget -r -np -nd --accept=gz http://hgdownload-test.cse.ucsc.edu/goldenPath/hg19/phastCons46way/placentalMammals/
                  for c in chr*wig*gz
                  do
                  bw=${c%phastCons46way.placental.wigFix.gz}bw
                  echo $bw
                  gunzip -c $c | wigToBigWig stdin chrom_len $bw  ## chrom_len is where you put your reference chromosome information file
                  done

.. _hg19_conserv: http://hgdownload-test.cse.ucsc.edu/goldenPath/hg19/phastCons46way/placentalMammals/
.. _mm9_conserv: http://hgdownload-test.cse.ucsc.edu/goldenPath/mm9/phastCons30way/placental/
.. _hg38_conserv: http://hgdownload-test.cse.ucsc.edu/goldenPath/hg38/phastCons7way/
.. _mm10_conserv: http://hgdownload-test.cse.ucsc.edu/goldenPath/mm10/phastCons60way/
.. _UCSC table browser: http://genome.ucsc.edu/cgi-bin/hgTables?command=start

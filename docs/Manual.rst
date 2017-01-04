=================
Instructions
=================


Through the pipeline, several temporary files will be generated, some of them are only used for settings
and transitions, others for continuing the next step, the rest for publishing and interpreting a biological
story.

.. note::

   a better way to organize your ChIP-seq project

Raw Data
========

.. _Raw Data:

========  =======  ======================================================
Format    Type     instruction
========  =======  ======================================================
FASTQ     text     single end
FASTQ.gz  gz       single end
FASTQ     text     pair end
FASTQ.gz  gz       pair end
========  =======  ======================================================

-----

.. _Manual:

Instructions to usage
=========================

Simple mode (The major mode)
---------------------------------------------

Demo data command is as follows::

   chilin  simple -p narrow -t foxa1_t1.fastq  -c foxa1_c1.fastq -i local -o local -s hg19  --skip 10,12 --dont_remove

See skip_ option for details.

This is major and the easiest mode to run ChiLin for single end data with default bwa mapper, for single end data using comma to separate sample replicates for IP and input ChIP-seq sample::

  chilin  simple -u your_name -s your_species --threads 8 -i id -o output -t treat1.fastq,treat2.fastq -c control1.fastq,control2.fastq  -p narrow -r tf

For pair end data, use semicolon to separate sample replicates, use comma to separate pairs, do not forget to add `quotes(")` of your sample file path::

  chilin simple --threads 8 -i H3K27me3_PairEnd -o H3K27me3_PairEnd -u you -s mm9 -t "GSM905438.fastq_R1.gz,GSM905438.fastq_R2.gz" -c "GSM905434.fastq_R1.gz,GSM905434.fastq_R2.gz;GSM905436.fastq_R1.gz,GSM905436.fastq_R2.gz" -p both --pe

See more options about *simple* by::

  chilin simple -h

.. _simple-mode:

simple mode useful option
-------------------------

* -t In simple mode, this is the options for specifying path to treatment.
* -c In simple mode, this is the options for specifying path to control.
* -p peaks calling type, narrow or broad or both, #e.g., If your factor is transcription factor, and you want to call narrow peaks only, or broad histone mark with broad peak calling
* -i prefix of the output name
* -o output results directory
* -s species, must be filled, the -s option is corresponding to your filled python config file [section], see section for details, :envvar:`[species] <[species]>`.
* -u user
* --pe pair end mode
* --pe pair end mode
* --maper to choose your mapping tool
* --threads threads number for mapping tool
* -r factor type, some default settings for three factor types, tf,histone,dnase

(optional) gen mode
------------------------------

This mode is to generate config file for `run-mode`_. A config file is look like this,

- The major section user needs to fill is the :envvar:`[basics] <[basics]>` section.

.. code-block:: bash

                chilin gen -o test.conf

.. code-block:: ini

		[basics]
		user = anonymous
		id = local
		time = 2014-05-09
		species = hg19
		factor = tf
		treat = foxa1_t1.fastq,foxa1_t1.fastq
		cont = foxa1_c1.fastq,foxa1_c1.fastq
		output = output_directory
		version = 2.0.0

run mode usage
------------------------------

.. _run-mode:

After configurating the config files above, you could use run mode with a single command::

  chilin  gen -o my_config
  ## modify tool parameters and run
  chilin  run -c my_config

batch mode usage
------------------------------

This mode help user run dataset one by one with one process.

After configurating a batch of the config files above, such as e.g. 1.conf, 2.conf, 3.conf, then you fill in a file called *batch.conf*::

  1.conf
  2.conf
  3.conf

you could use batch mode with a single command::

  chilin  batch -b batch.conf

.. _skip:

Step control, mapper choice and threads control, mimic run
----------------------------------------------------------

Common options can be used for simple mode, run and batch modes.
Each step control is tolerant, continue running even tool failed processing.

- --skip, step control, e.g::

    chilin  simple -s hg19 -i id -p narrow -o output -u user --skip 1,3,5,9,10,11 -t treat1.fastq,treat2.fastq`

  - step 1(*can skip*): FastQC sequence quality evaluation, reads GC contents evaluation and library contamination, this step can be skipped.
  - step 2: bwa(default), bowtie or star mapping, this step cannot be skipped, because this provided necessary BAM files.
  - step 3: sub-sample bam files and do macs2 fragment size estimation.
  - step 4: sub-sample bam files(if step 3 is run, skip) and do PBC evaluation
  - step 5: call peak for replicates samples, and do replicates peak overlap/correlation analysis
  - step 6: call peak for merged bam file, this step cannot be skipped, because this provided peak for annotation step
  - step 7: calculate FRiP scores for each sample.
  - step 8: use bedAnnotate.py script to evaluate merged peak calling meta regions distribution(promoters, exons, introns, intergenic, dhs, black list regions)
  - step 9: sub-sample bam files(if step 3 is run, skip) and calculate reads ratio in meta regions
  - step 10(*can skip*): draw Phastcon scores distribution around peak call summits, if you do not have Phastcon score bigwig files, use --skip 10 or leave chilin.conf blank for that reference
  - step 11(*can skip*): Regulatory potential score calculation on top 10k peaks
  - step 12(*can skip*): use MDSeqPos to perform motif analysis, for peaks number less than 200, ChiLin skip this step automatically
  - step 13(*can skip*): generate report, use explicitly skip to skip 10-12 if necessary.

- --dont_resume, by default, each re-run would use previous temporary files to resume from the step it crashed. When dont_resume is on, ChiLin would start from first step, so user do not to clean up the work directory.
- --dont_remove, keep temporary files
- --dry-run, mimic run chilin command
- --threads, BWA, Bowtie and FastQC multithreads options.
- --mapper, to choose mapping tools, should match your genome index in :envvar:`genome index<genome_index>`

Instructions to config file
==============================

basics
----------------------------

.. envvar:: [basics]

    Lists all the meta-data of current workflow.
    Consist of the following options:

    .. envvar:: user

	user name

    .. envvar:: time

	time you start to run


    .. envvar:: species

        The name of species, written to the QCreport and log
        Limit: a string (1) consist of ``numbers``, ``alphabets`` or ``'_'`` (2) shorter than 20 characters

    .. envvar:: id

        This is used as output prefix, such as input id: test, output file would be: test_treat.bam

    .. envvar:: factor

        The name of species, writen to DC summary and QCreport, log
        Limit: a string (1) come from GO standard term


    .. envvar:: treat

       The paths of treatment files
       Limit: absolute ``path`` of files in :ref:`supported formats<Raw Data>`
		

    .. envvar:: cont

       The paths of control files
       Limit: absolute or relative ``path`` of files in :ref:`supported formats<Raw Data>`

    .. envvar:: output

       The paths of output directory.
		

tool
----------------------------
The tool section is like this:

.. literalinclude:: ../chilin.conf
   :language: ini
   :lines: 8-16
   :emphasize-lines: 13-15
   :linenos:

.. envvar:: [tool]

    Lists all the meta-data of current workflow.
    Consist of the following options:

    .. envvar:: mdseqpos

       absolute path to ``MDSeqPos.py``

    .. envvar:: macs2

       absolute path to ``macs2``


species
----------------------------

You can add as many species as possible. To add species, first you need to read :ref:`dependent data <dependentdata>` section to fill the following.
Then, you should fill the config files species section, the rule is like follows, e.g. hg19 assembly.

.. literalinclude:: hg19.conf
   :language: ini
   :lines: 1-10
   :emphasize-lines: 1-10
   :linenos:

And mm9 assembly,

.. literalinclude:: hg19.conf
   :language: ini
   :lines: 12-21
   :emphasize-lines: 12-21
   :linenos:

And hg38 assembly,

.. literalinclude:: hg19.conf
   :language: ini
   :lines: 23-32
   :emphasize-lines: 23-32
   :linenos:

.. envvar:: [species]

    specific species assembly version you want to analyze
    Consist of the following options:

    .. envvar:: genome_index

       absolute path to corresponding mappers genome index, if you use default bwa, this should be
       bwa index.

    .. envvar:: genome_dir

       absolute path to genome fasta files, separated by chromosome, like `chr1.fa, chr2.fa, chr3.fa ...`

    .. envvar:: chrom_len

       absolute path to chromosome length text file

    .. envvar:: dhs

       absolute path union DHS regions

    .. envvar:: velcro

       absolute path black list regions

    .. envvar:: conservation

       absolute path to the directory containing UCSC Phastcon score bigwig files

    .. envvar:: geneTable

       standard refSeq annotation table from UCSC table browser


contamination
----------------------------

you can add all species you are suspicious of sampling swap or library contamination.

.. literalinclude:: ../chilin.conf
   :language: ini
   :lines: 52-61
   :emphasize-lines: 52-61
   :linenos:

.. envvar:: [contamination]

   specific species assembly path that you want to screen.


software options
----------------------------
ChiLin has some user-defined parameters for macs2, regulatory potential, conservation score and motif analysis.

.. literalinclude:: ../chilin.conf
   :language: ini
   :lines: 22-51
   :emphasize-lines: 23-50
   :linenos:
   
.. envvar:: [macs2]

   macs2 parameters
	    
    .. envvar:: extsize

		fixed extension size for macs2 peak calling
       
    .. envvar:: type [for macs2]

		peak calling types, user can choose narrow, broad and both, we suggest user use narrow for TF and active histone marks, use broad for broad histone marks, use both for chromatin regulators.
    

    .. envvar:: fdr

		FDR cutoff for macs2 peak calling

    .. envvar:: keep_dup

		duplicates level, suggest 1 for removing redundancy, or all for preserving all redudancy for DNase-seq
       

.. envvar:: [reg]

   specific species assembly version you want to analyze
	    
    .. envvar:: peaks [for reg]

                top peaks number(sorted by macs2 score) for estimating regulatory score.
       
    .. envvar:: dist

                distance to TSS cutoff when calculating regulatory potential score


.. envvar:: [conservation]

   specific species assembly version you want to analyze
	    
    .. envvar:: type [for conservation]

                transcription factor or histone mark
       
    .. envvar:: peaks

		top peaks (sorted by macs2 score) for plotting conservation distribution
       
    .. envvar:: width

		window width around peaks summits for plotting conservation

.. envvar:: [seqpos]

   specific species assembly version you want to analyze
	    
    .. envvar:: peaks [for conservation]

                top peaks number for seqpos (search in the motif database)
       
    .. envvar:: mdscan_width

                motif scan window width around peak summit
       
    .. envvar:: mdscan_top_peaks

		top peaks for denovo motif scan

    .. envvar:: width

		seqpos width

    .. envvar:: seqpos_mdscan_top_peaks_refine

                seqpos and mdscan top peaks refine, see mdseqpos

    .. envvar:: db

                choose mdseqpos motif database, default cistrome.xml

    .. envvar:: pvalue_cutoff

                cutoff for motif analysis

.. _Instructions_table:
.. _Instructions_results:


Instructions to results
=========================

The output prefix is from:

* `simple-mode`_ -i id specified, or `run-mode`_ filled in :envvar:`[basics] <[basics]>` section :envvar:`id <id>` part.
* The output directory is `simple mode` -o output specified or 2. `run-mode`_  filled in :envvar:`[basics] <[basics]>` section :envvar:`output <output>` part.
* For a fully test dataset with replicates of treatments and replicates of controls, the results folder are like following, which are generated with *-dont_remove* option, to preserve all temporary files, use *-dont_remove* option.

.. _cleaner:

Without `--dont_remove` option, the work directory would be cleaned up::

  id
  |-- attic
  |   |-- json
  |   |   |-- id_conserv.json
  |   |   |-- id_contam.json
  |   |   |-- id_dhs.json
  |   |   |-- id_enrich_meta.json
  |   |   |-- id_fastqc.json
  |   |   |-- id_frag.json
  |   |   |-- id_frip.json
  |   |   |-- id_macs2.json
  |   |   |-- id_macs2_rep.json
  |   |   |-- id_map.json
  |   |   |-- id_meta.json
  |   |   |-- id_pbc.json
  |   |   `-- id_rep.json
  |   |-- id_conserv.pdf
  |   |-- id_control.bam
  |   |-- id_control_rep1.bam
  |   |-- id_control_rep2.bam
  |   |-- id_gene_score.txt
  |   |-- id_treat_rep1.bam
  |   |-- id_treat_rep2.bam
  |   `-- id_treatment.bam
  |-- id_report.pdf
  |-- id_control.bw
  |-- id_peaks.xls
  |-- id_sort_peaks.narrowPeak
  |-- id_sort_summits.bed
  |-- id_treat.bw
  |-- id_treat_rep1_control.bw
  |-- id_treat_rep1_peaks.xls
  |-- id_treat_rep1_sort_peaks.narrowPeak
  |-- id_treat_rep1_treat.bw
  |-- id_treat_rep2_control.bw
  |-- id_treat_rep2_peaks.xls
  |-- id_treat_rep2_sort_peaks.narrowPeak
  `-- id_treat_rep2_treat.bw


With `--dont_remove` option,

   .. code-block:: bash

                   output
                   |-- json  ## qc statistics
                   |   |-- id_conserv.json  ## conservation scores
                   |   |-- id_contam.json   ## library contamination evaluation
                   |   |-- id_dhs.json      ## union dhs overlap
                   |   |-- id_enrich_meta.json ## meta regions reads ratio
                   |   |-- id_fastqc.json      ## fastqc evaluation
                   |   |-- id_frag.json        ## fragment size evaluation
                   |   |-- id_frip.json        ## FRiP scores
                   |   |-- id_macs2.json       ## merged macs2 peak calling number
                   |   |-- id_macs2_rep.json   ## macs2 replicates peaks number
                   |   |-- id_map.json         ## mapping ratio statistics
                   |   |-- id_meta.json        ## peak meta regions distribution
                   |   |-- id_pbc.json         ## PBC score
                   |   `-- id_rep.json         ## replicates consistency
                   |-- latex  ## rendered latex document
                   |   |-- id_begin.tex       
                   |   |-- id_conserv.tex
                   |   |-- id_contam.tex
                   |   |-- id_end.tex
                   |   |-- id_fastqc.tex
                   |   |-- id_fastqc_gc.tex
                   |   |-- id_frip.tex
                   |   |-- id_map.tex
                   |   `-- id_summary_table.tex
                   |-- id.aux ## latex log file
                   |-- id.cor ## correlation analysis temporary file
                   |-- id.dhs ## dhs overlap analysis temporary file
                   |-- id.log ## latex log file
                   |-- id.meta ## meta regions peak distribution temporary file
                   |-- id.out  ## latex log file
                   |-- id_report.pdf  ## pdf document generated
                   |-- id.tex  ## file latex file
                   |-- id_0_1.overlap  ## replicates peak overlap
                   |-- id_bwa_compare.R ## R script for comparing new data to historic data
                   |-- id_bwa_compare.pdf ## pdf generated by R script above
                   |-- id_conserv.R  ## conservation plot R code
                   |-- id_conserv.pdf ## pdf generated by R script above
                   |-- id_conserv.txt ## 7 or 5 point conservation scores around summits
                   |-- id_conserv_cluster.R ## conservation scores clustering plot
                   |-- id_conserv_compare.pdf  ## conservation pdf generated by R script above
                   |-- id_conserv_img.pdf  ## low resolution image of conservation plot
                   |-- id_control.bam      ## merged control bam files
                   |-- id_control.bw       ## control bigwiggle file
                   |-- id_control_lambda.bdg  ## control bedgraph file
                   |-- id_control_lambda.bdg.tmp  ## bedClip filtered bedgraph file
                   |-- id_control_rep1.bam   ## sorted, mapping quality 1 filtered replicate 1st bam file 
                   |-- id_control_rep1.enrich.dhs   ## reads ratio in DHS regions
                   |-- id_control_rep1.enrich.exon  ## reads ratio in exon regions
                   |-- id_control_rep1.enrich.promoter  ## reads ratio in promoter regions
                   |-- id_control_rep1.fastq  ## copied fastq file
                   |-- id_control_rep1.frip   ## FRiP score from replicate control 1st
                   |-- id_control_rep1.hist   ## read locations histogram of replicate control 1st 
                   |-- id_control_rep1.nochrM  ## chromosome information without chrM
                   |-- id_control_rep1.pbc  ## bwa PBC score
                   |-- id_control_rep1.sai  ## bwa sai file
                   |-- id_control_rep1.sam  ## bwa sam file
                   |-- id_control_rep1.tmp.bam     ## mapping quality filtered bam files, without sorting
                   |-- id_control_rep1_100k.fastq  ## subsampled fastq reads
                   |-- id_control_rep1_100k_fastqc ## fastqc temporary results
                   |   |-- Icons
                   |   |   |-- error.png
                   |   |   |-- fastqc_icon.png
                   |   |   |-- tick.png
                   |   |   `-- warning.png
                   |   |-- Images
                   |   |   |-- duplication_levels.png
                   |   |   |-- kmer_profiles.png
                   |   |   |-- per_base_gc_content.png
                   |   |   |-- per_base_n_content.png
                   |   |   |-- per_base_quality.png
                   |   |   |-- per_base_sequence_content.png
                   |   |   |-- per_sequence_gc_content.png
                   |   |   |-- per_sequence_quality.png
                   |   |   `-- sequence_length_distribution.png
                   |   |-- fastqc_data.txt
                   |   |-- fastqc_report.html
                   |   `-- summary.txt
                   |-- id_control_rep1_100k_fastqc.zip 
                   |-- id_control_rep1_4000000.bam   ## subsampled 4M reads bam file
                   |-- id_control_rep1_4000000_nochrM.bam  ## subsampled non-chrM 4M reads bam file
                   |-- id_control_rep1_mapped.bwa  ## replicate control 1st mapped reads statistics
                   |-- id_control_rep1_nochrM.bam  ## sorted, mapping quality filtered bam file
                   |-- id_control_rep1_nochrM.sam  ## mapped sam files without chrM
                   |-- id_control_rep1_nochrM.sam.4000000  ## subsampled 4M reads without chrM
                   |-- id_control_rep1_total.bwa  ## total reads statistics from bwa
                   |-- id_control_rep1_u.sam  ## unique reads SAM file
                   |-- id_control_rep1_u.sam.4000000  ## subsampled unique reads SAM file
                   |-- id_control_rep1mbr.bam  ## cross species mapping to mbr, or species you specified
                   |-- id_control_rep1mbr.sai  
                   |-- id_control_rep1mbr.sam
                   |-- id_control_rep1mbr.tmp.bam
                   |-- id_control_rep1mbr_mapped.bwa
                   |-- id_control_rep1mbr_total.bwa
                   |-- id_control_rep2.bam   ## control replicates 2nd bam file
                   |-- id_control_rep2.enrich.dhs 
                   |-- id_control_rep2.enrich.exon
                   |-- id_control_rep2.enrich.promoter
                   |-- id_control_rep2.fastq
                   |-- id_control_rep2.frip
                   |-- id_control_rep2.hist
                   |-- id_control_rep2.nochrM
                   |-- id_control_rep2.pbc
                   |-- id_control_rep2.sai
                   |-- id_control_rep2.sam
                   |-- id_control_rep2.tmp.bam
                   |-- id_control_rep2_100k.fastq
                   |-- id_control_rep2_100k_fastqc
                   |   |-- Icons
                   |   |   |-- error.png
                   |   |   |-- fastqc_icon.png
                   |   |   |-- tick.png
                   |   |   `-- warning.png
                   |   |-- Images
                   |   |   |-- duplication_levels.png
                   |   |   |-- kmer_profiles.png
                   |   |   |-- per_base_gc_content.png
                   |   |   |-- per_base_n_content.png
                   |   |   |-- per_base_quality.png
                   |   |   |-- per_base_sequence_content.png
                   |   |   |-- per_sequence_gc_content.png
                   |   |   |-- per_sequence_quality.png
                   |   |   `-- sequence_length_distribution.png
                   |   |-- fastqc_data.txt
                   |   |-- fastqc_report.html
                   |   `-- summary.txt
                   |-- id_control_rep2_100k_fastqc.zip
                   |-- id_control_rep2_4000000.bam
                   |-- id_control_rep2_4000000_nochrM.bam
                   |-- id_control_rep2_mapped.bwa
                   |-- id_control_rep2_nochrM.bam
                   |-- id_control_rep2_nochrM.sam
                   |-- id_control_rep2_nochrM.sam.4000000
                   |-- id_control_rep2_total.bwa
                   |-- id_control_rep2_u.sam
                   |-- id_control_rep2_u.sam.4000000
                   |-- id_control_rep2mbr.bam
                   |-- id_control_rep2mbr.sai
                   |-- id_control_rep2mbr.sam
                   |-- id_control_rep2mbr.tmp.bam
                   |-- id_control_rep2mbr_mapped.bwa
                   |-- id_control_rep2mbr_total.bwa
                   |-- id_gene_score.txt  ## regulatory potential for top 10000 peaks
                   |-- id_peaks.narrowPeak  ## merged peak call for narrowPeak or broadPeak
                   |-- id_peaks.xls ## macs2 excel file 
                   |-- id_peaks_top_conserv.bed  ## top peaks for conservation plot
                   |-- id_peaks_top_reg.bed  ## top peaks for regulatory potential score calculation
                   |-- id_raw_sequence_qc.R  ## median raw sequence quality plot
                   |-- id_raw_sequence_qc.pdf
                   |-- id_sort_peaks.narrowPeak ## sorted merged peak calling
                   |-- id_sort_summits.bed  ## sorted summits of peaks
                   |-- id_summary.txt  ## plain text for qc summary
                   |-- id_summits.bed  ## merged peak calling summits file
                   |-- id_treat.bw  ## merged pileup treatment bigwiggle file
                   |-- id_treat_pileup.bdg  ## merged pileup treatment bedgraph file
                   |-- id_treat_pileup.bdg.tmp  ## merged pileup treatment bedgraph temporary file
                   |-- id_treat_rep1   ## MACS2 predictd R script
                   |-- id_treat_rep1.bam  ## bam file generated by bwa and samtools
                   |-- id_treat_rep1.enrich.dhs
                   |-- id_treat_rep1.enrich.exon
                   |-- id_treat_rep1.enrich.promoter
                   |-- id_treat_rep1.fastq
                   |-- id_treat_rep1.frip
                   |-- id_treat_rep1.hist
                   |-- id_treat_rep1.nochrM
                   |-- id_treat_rep1.pbc
                   |-- id_treat_rep1.sai
                   |-- id_treat_rep1.sam
                   |-- id_treat_rep1.tmp.bam
                   |-- id_treat_rep1_100k.fastq
                   |-- id_treat_rep1_100k_fastqc
                   |   |-- Icons
                   |   |   |-- error.png
                   |   |   |-- fastqc_icon.png
                   |   |   |-- tick.png
                   |   |   `-- warning.png
                   |   |-- Images
                   |   |   |-- duplication_levels.png
                   |   |   |-- kmer_profiles.png
                   |   |   |-- per_base_gc_content.png
                   |   |   |-- per_base_n_content.png
                   |   |   |-- per_base_quality.png
                   |   |   |-- per_base_sequence_content.png
                   |   |   |-- per_sequence_gc_content.png
                   |   |   |-- per_sequence_quality.png
                   |   |   `-- sequence_length_distribution.png
                   |   |-- fastqc_data.txt
                   |   |-- fastqc_report.html
                   |   `-- summary.txt
                   |-- id_treat_rep1_100k_fastqc.zip
                   |-- id_treat_rep1_4000000.bam
                   |-- id_treat_rep1_4000000_nochrM.bam
                   |-- id_treat_rep1_control.bw
                   |-- id_treat_rep1_control_lambda.bdg
                   |-- id_treat_rep1_control_lambda.bdg.tmp
                   |-- id_treat_rep1_frag_sd.R  ## fragment analysis script for parsing macs2 R script
                   |-- id_treat_rep1_mapped.bwa
                   |-- id_treat_rep1_model.R  ## MACS2 R script for analyzing fragment size
                   |-- id_treat_rep1_nochrM.bam
                   |-- id_treat_rep1_nochrM.sam
                   |-- id_treat_rep1_nochrM.sam.4000000
                   |-- id_treat_rep1_peaks.narrowPeak  ## replicate 1 peak calling
                   |-- id_treat_rep1_peaks.xls
                   |-- id_treat_rep1_sort_peaks.narrowPeak
                   |-- id_treat_rep1_summits.bed
                   |-- id_treat_rep1_total.bwa
                   |-- id_treat_rep1_treat.bw
                   |-- id_treat_rep1_treat_pileup.bdg
                   |-- id_treat_rep1_treat_pileup.bdg.tmp
                   |-- id_treat_rep1_u.sam  ## uniquely mapping sam file, defined by mapping quality above 1
                   |-- id_treat_rep1_u.sam.4000000
                   |-- id_treat_rep1mbr.bam
                   |-- id_treat_rep1mbr.sai
                   |-- id_treat_rep1mbr.sam
                   |-- id_treat_rep1mbr.tmp.bam
                   |-- id_treat_rep1mbr_mapped.bwa
                   |-- id_treat_rep1mbr_total.bwa
                   |-- id_treat_rep2
                   |-- id_treat_rep2.bam
                   |-- id_treat_rep2.enrich.dhs
                   |-- id_treat_rep2.enrich.exon
                   |-- id_treat_rep2.enrich.promoter
                   |-- id_treat_rep2.fastq
                   |-- id_treat_rep2.frip
                   |-- id_treat_rep2.hist
                   |-- id_treat_rep2.nochrM
                   |-- id_treat_rep2.pbc
                   |-- id_treat_rep2.sai
                   |-- id_treat_rep2.sam
                   |-- id_treat_rep2.tmp.bam
                   |-- id_treat_rep2_100k.fastq
                   |-- id_treat_rep2_100k_fastqc
                   |   |-- Icons
                   |   |   |-- error.png
                   |   |   |-- fastqc_icon.png
                   |   |   |-- tick.png
                   |   |   `-- warning.png
                   |   |-- Images
                   |   |   |-- duplication_levels.png
                   |   |   |-- kmer_profiles.png
                   |   |   |-- per_base_gc_content.png
                   |   |   |-- per_base_n_content.png
                   |   |   |-- per_base_quality.png
                   |   |   |-- per_base_sequence_content.png
                   |   |   |-- per_sequence_gc_content.png
                   |   |   |-- per_sequence_quality.png
                   |   |   `-- sequence_length_distribution.png
                   |   |-- fastqc_data.txt
                   |   |-- fastqc_report.html
                   |   `-- summary.txt
                   |-- id_treat_rep2_100k_fastqc.zip
                   |-- id_treat_rep2_4000000.bam
                   |-- id_treat_rep2_4000000_nochrM.bam
                   |-- id_treat_rep2_control.bw
                   |-- id_treat_rep2_control_lambda.bdg
                   |-- id_treat_rep2_control_lambda.bdg.tmp
                   |-- id_treat_rep2_frag_sd.R
                   |-- id_treat_rep2_mapped.bwa
                   |-- id_treat_rep2_model.R
                   |-- id_treat_rep2_nochrM.bam
                   |-- id_treat_rep2_nochrM.sam
                   |-- id_treat_rep2_nochrM.sam.4000000
                   |-- id_treat_rep2_peaks.narrowPeak
                   |-- id_treat_rep2_peaks.xls
                   |-- id_treat_rep2_sort_peaks.narrowPeak
                   |-- id_treat_rep2_summits.bed
                   |-- id_treat_rep2_total.bwa
                   |-- id_treat_rep2_treat.bw
                   |-- id_treat_rep2_treat_pileup.bdg
                   |-- id_treat_rep2_treat_pileup.bdg.tmp
                   |-- id_treat_rep2_u.sam   ## uniquely mapping sam file, defined by mapping quality above 1
                   |-- id_treat_rep2_u.sam.4000000
                   |-- id_treat_rep2mbr.bam
                   |-- id_treat_rep2mbr.sai
                   |-- id_treat_rep2mbr.sam
                   |-- id_treat_rep2mbr.tmp.bam
                   |-- id_treat_rep2mbr_mapped.bwa
                   |-- id_treat_rep2mbr_total.bwa
                   `-- id_treatment.bam  ## samtools merged filtered bam files



Built-in tools
=================
* conservation_plot.py for generating conservation profiles
* bedAnnotate.py is used to calculate meta gene distribuiton.

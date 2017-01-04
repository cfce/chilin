Cistrome ChiLin
================
It is a python package for one-in-all solution of processing ChIP-seq and DNase-seq data.

Quick Start
===============

See if you have gcc, g++, java, R, python-dev installed (http://cistrome.org/chilin/Installation.html#dependent-software-list).

First, clone:

``` sh
git clone https://github.com/cfce/chilin && cd chilin
```

then install through:

``` sh
python setup.py clean && python setup.py install -f
```

source virtual environment and use:

``` sh
source chilin_env/bin/activate
chilin -h
```

fetch `hg19` reference data, and test on `demo` data:

under the root of the chilin code.

``` sh
# change to default directory
mkdir -p db

cd db

# all hg19 reference data
wget -c http://cistrome.org/chilin/_downloads/hg19.tgz
wget -c http://cistrome.org/chilin/_downloads/hg19.tgz.md5 ## check md5
md5sum -c hg19.tgz
tar xvfz hg19.tgz
# download mycoplasma for judgement of contamination in your samples
wget -c http://cistrome.org/chilin/_downloads/mycoplasma.tgz
wget -c http://cistrome.org/chilin/_downloads/mycoplasma.tgz.md5
md5sum -c mycoplasma.tgz.md5
tar xvfz mycoplasma.tgz

# check all installation
cd .. && python setup.py -l
cd demo && bash foxa1
```

Usage
==============================

Demo data command is as follows:

``` sh
   chilin  simple -p narrow -t foxa1_t1.fastq  -c foxa1_c1.fastq -i local -o local -s hg19  --skip 10,12 --dont_remove
```

See skip_ option for details.

This is major and the easiest mode to run ChiLin for single end data with default bwa mapper, for single end data using comma to separate sample replicates for IP and input ChIP-seq sample:

``` sh

  chilin  simple -u your_name -s your_species --threads 8 -i id -o output -t treat1.fastq,treat2.fastq -c control1.fastq,control2.fastq  -p narrow -r tf
```

For pair end data, use semicolon to separate sample replicates, use comma to separate pairs, do not forget to add `quotes(")` of your sample file path:

``` sh
  chilin simple --threads 8 -i H3K27me3_PairEnd -o H3K27me3_PairEnd -u you -s mm9 -t "GSM905438.fastq_R1.gz,GSM905438.fastq_R2.gz" -c "GSM905434.fastq_R1.gz,GSM905434.fastq_R2.gz;GSM905436.fastq_R1.gz,GSM905436.fastq_R2.gz" -p both --pe
```

Currently, only bwa support pair end processing. bwa supports both fastq.gz and fastq file, bowtie only support fastq file, the pipeline should use the corresponding aligner's genome index configured in the [configuration files](http://cistrome.org/chilin/Manual.html#species).

Update the configuration
==============================
If you modify the code or update any part of the configuration file *chilin.conf.filled*, such as different aligner's genome index, union DHS BED file, reinstall the package itself only.

``` sh
source chilin_env/bin/activate && python setup.py install 
```

Uninstall
===============

``` sh
python setup.py clean
deactivate
```

Troubleshooting
==================

If any error of the dependent software occur, try to upgrade the corresponding software. 
Those *warnings* generated in *pdflatex* step is ok.
There is one known issue of mm9 chrom-info in CentOS. ChiLin is suggested to be used under Ubuntu.
If *sys_platform* error occurs, uninstall the system setuptools and install the latest setuptools manually.

Documentation
================
full documentation: http://cistrome.org/chilin

github wiki: https://github.com/cfce/chilin/wiki


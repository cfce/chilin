Cistrome ChiLin
================
It is a python package for one-in-all solution of processing ChIP-seq and DNase-seq data.

Quick Start
===============

See if you have gcc, g++, java, R, python-dev installed (http://cistrome.org/chilin/Installation.html#dependent-software-list).

then install through:

``` sh
git clone https://github.com/cfce/chilin && cd chilin
python setup.py clean && python setup.py install -f
```

source virtual environment and use:

``` sh
source chilin_env/bin/activate
chilin -h
```

fetch `hg19` reference data, and test on `demo` data:

``` sh
# change to default directory
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

For details about the [dependency data](https://github.com/cfce/chilin/wiki/Appendix#data-details) and [software](http://cistrome.org/chilin/Installation.html#dependent-software-list)

Documentation
================
full documentation: http://cistrome.org/chilin
github wiki: https://github.com/cfce/chilin/wiki


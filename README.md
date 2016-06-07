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

Uninstall
===============

``` sh
python setup.py clean
deactivate
```

Documentation
================
full documentation: http://cistrome.org/chilin

github wiki: https://github.com/cfce/chilin/wiki


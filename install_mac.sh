#!/bin/bash

wget -c --no-check-certificate https://raw.githubusercontent.com/pypa/virtualenv/1.9.X/virtualenv.py
python virtualenv.py -p/opt/local/Library/Frameworks/Python.framework/Versions/2.7/Resources/Python.app/Contents/MacOS/Python --system-site-packages --distribute chilin_env
source chilin_env/bin/activate

git clone https://github.com/lh3/seqtk seqtk &>/dev/null
cd seqtk && make && chmod 755 seqtk && cp seqtk ../chilin_env/bin && cd ..

which bedGraphToBigWig || wget -c http://hgdownload.cse.ucsc.edu/admin/exe/macOSX.x86_64/bedGraphToBigWig &>/dev/null
wget -c http://hgdownload.cse.ucsc.edu/admin/exe/macOSX.x86_64/bedClip
wget -c http://hgdownload.cse.ucsc.edu/admin/exe/macOSX.x86_64/wigCorrelate
wget -c http://hgdownload.cse.ucsc.edu/admin/exe/macOSX.x86_64/wigToBigWig
chmod 755 bedGraphToBigWig bedClip wigCorrelate wigToBigWig && cp bedGraphToBigWig bedClip wigCorrelate wigToBigWig chilin_env/bin

wget -c http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.10.1.zip
unzip fastqc_v0.10.1.zip
chmod -R 755 FastQC
cp -r FastQC/* chilin_env/bin

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

wget -c --no-check-certificate https://github.com/taoliu/MACS/archive/de419af20b6b441f0c0b62a3c18a262228409c8a.zip -O macs_2.0.10.2014XXXX.zip 
unzip  macs_2.0.10.2014XXXX.zip 
cd MACS*
python setup_w_cython.py build_ext
python setup.py install
R -e "source(\"http://bioconductor.org/biocLite.R\"); biocLite(\"seqLogo\"); library(seqLogo) "

wget -c --no-check-certificate https://bitbucket.org/cistrome/cistrome-applications-harvard/get/cfec0147b6f6.zip
unzip cfec0147b6f6.zip
cd cistrome*/mdseqpos
cp lib/settings.py.example lib/settings.py

echo "please install mdseqpos manually, https://bitbucket.org/cistrome/cistrome-applications-harvard/src/270e986a26a0ab3b143dbf9720d990f3654f677f/mdseqpos/?at=default"



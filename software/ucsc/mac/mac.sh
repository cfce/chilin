#!/bin/bash

which bedGraphToBigWig || wget -c http://hgdownload.cse.ucsc.edu/admin/exe/macOSX.x86_64/bedGraphToBigWig &>/dev/null
wget -c http://hgdownload.cse.ucsc.edu/admin/exe/macOSX.x86_64/bedClip
wget -c http://hgdownload.cse.ucsc.edu/admin/exe/macOSX.x86_64/wigCorrelate
wget -c http://hgdownload.cse.ucsc.edu/admin/exe/macOSX.x86_64/wigToBigWig



#!/bin/bash -ex
#wget -c https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
__conda_setup="$('/home/qinqian/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/home/qinqian/miniconda3/etc/profile.d/conda.sh" ]; then
        . "/home/qinqian/miniconda3/etc/profile.d/conda.sh"
    else
        export PATH="/home/qinqian/miniconda3/bin:$PATH"
    fi
fi
unset __conda_setup

conda create -n chilin_env python=2
conda activate chilin_env

conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
conda install bwa=0.7.13 samtools=0.1.19 bedtools=2.17.0 seqtk ucsc-bedclip ucsc-bedgraphtobigwig ucsc-wigcorrelate ucsc-wigtobigwig fastqc numpy macs2=2.1.0 bioconductor-seqlogo 

easy_install pip

sudo apt install cm-super
updmap

### UCSC hg38
#wget -c http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
#wget -c http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.masked.gz
#
### NCBI GRCH38
### http://asia.ensembl.org/info/data/ftp/index.html
wget -c ftp://ftp.ensembl.org/pub/release-97/fasta/homo_sapiens/dna//Homo_sapiens.GRCh38.dna.toplevel.fa.gz

## bwakit similar to UCSC version..
#wget -c ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_analysis_set.fna.gz


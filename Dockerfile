FROM continuumio/miniconda3

RUN apt-get update -qq

RUN apt-get install -y texlive-latex-base
RUN apt-get install -y texlive-latex-extra # --no-install-recommends
RUN apt-get install -y texlive-fonts-recommended # --no-install-recommends
RUN apt-get install -y texlive-plain-generic # --no-install-recommends
RUN apt-get clean -y && apt-get autoremove -y

RUN conda config --add channels conda-forge
RUN conda config --add channels bioconda
RUN conda install mamba -c conda-forge/label/mamba-alpha -c conda-forge
RUN mamba create -y -q -n chilin python=2.7.15 
SHELL ["conda", "run", "-n", "chilin", "/bin/bash", "-c"]
RUN mamba install -y -q -c anaconda cython
SHELL ["conda", "run", "-n", "chilin", "/bin/bash", "-c"]
RUN mamba install -y -q -c bioconda bx-python bwa=0.7.13 samtools=0.1.19 bedtools=2.17.0 seqtk ucsc-bedclip ucsc-bedgraphtobigwig ucsc-wigcorrelate ucsc-wigtobigwig fastqc numpy macs2=2.1.0 bioconductor-seqlogo

SHELL ["conda", "run", "-n", "chilin", "/bin/bash", "-c"]
RUN git clone https://github.com/cfce/chilin && cd chilin && python setup.py install && cp chilin.conf.filled /opt/conda/envs/chilin/lib/python2.7/site-packages/chilin2-2.0.0-py2.7.egg/chilin2/modules/config/chilin.conf && cp chilin.conf.filled /opt/conda/envs/chilin/lib/python2.7/site-packages/chilin2-2.0.0-py2.7.egg/chilin2/modules/config/chilin.conf.filled


RUN apt-get install -y ghostscript
RUN apt-get install -y imagemagick --fix-missing
RUN apt install -y ncbi-entrez-direct

RUN esearch -db nucleotide -query "NC_000908.2" | efetch -format fasta > mycoplasma.fasta && bwa index -a bwtsw mycoplasma.fasta
RUN esearch -db nucleotide -query "NC_000913.3" | efetch -format fasta > ecoli.fasta && bwa index -a bwtsw ecoli.fasta

WORKDIR /

ENTRYPOINT [ "conda", "run", "-n", "chilin", "chilin" ]

#!/bin/bash -ex

chilin=`pwd`/
#wget -c http://cistrome.org/chilin/_downloads/hg38.tgz
#tar xvfz hg38.tgz

cd hg38
wget -c https://www.ncbi.nlm.nih.gov/nuccore/108885074?report=fasta&log$=seqview&format=text

cd demo
docker run --memory="8000m" --cpus="1" -v ${chilin}:/chilin/db/ chilin:1.2 simple -u qinq -s hg38 --threads 1 -i local -o /chilin/ -t /chilin/db/demo/foxa1_t1.fastq.gz,/chilin/db/demo/foxa1_t1.fastq.gz -c /chilin/db/demo/foxa1_c1.fastq.gz,/chilin/db/demo/foxa1_c1.fastq.gz -p narrow -r tf --dont_remove

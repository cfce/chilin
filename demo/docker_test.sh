#!/bin/bash -ex

wget -c http://cistrome.org/chilin/_downloads/hg38.tgz
tar xvfz hg38.tgz

chilin=`pwd`/
cd demo

docker run -it --memory="10g" --cpus="2" -v ${chilin}:/chilin/db/ chilin:1.0 simple -u qinq -s hg38 --threads 2 -i local -o /chilin/ -t /chilin/db/demo/foxa1_t1.fastq.gz,/chilin/db/demo/foxa1_t1.fastq.gz -c /chilin/db/demo/foxa1_c1.fastq.gz,/chilin/db/demo/foxa1_c1.fastq.gz -p narrow -r tf --dont_remove

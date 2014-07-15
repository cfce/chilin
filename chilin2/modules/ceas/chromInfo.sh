#!/bin/bash
# this script prepare data for samtools and bedtools usage
ORG=$1
mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e \
        "select chrom, size from ${ORG}.chromInfo" > ${ORG}.genome

mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e \
        "select chrom, size from ${ORG}.chromInfo" | awk '{OFS="\t"; print $1,1,$2}' | sort -V -k1,1 -k2,2n > ${ORG}.bed

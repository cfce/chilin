#!/bin/env sh

## this script prepare meta information for regulatory potential score and ceas 

ORG=$1

### suppress header, for regulatory potential script
mysql -N -B --user=genome --host=genome-mysql.cse.ucsc.edu -A -D $ORG -e "SELECT name, chrom, strand, txStart, txEnd, name2 from refGene;" > ${ORG}.reg

## for ceas bedannotate script
mysql -N -B --user=genome --host=genome-mysql.cse.ucsc.edu -A -D $ORG -e "SELECT * from refGene;"  > ${ORG}.refGene

#!/bin/bash -ex
## exon.bed and gene.bed is parsed from standard UCSC RefSeq file
## this script prepare annotation for reads ratio in meta regions evaluation

# Start with the exons and genes BED file from UCSC(http://www.biostars.org/p/17985/)
# 1. Use BEDTools to merge overlapping exons in this file. The resulting merged exons for each gene make a sort of 'fake' transcript. Some people refer to these as the 'exon content' of a gene, the 'squashed' transcriptome, etc. Think about how you want to deal with transcripts that overlap on opposite strands at this point.
# exons and introns bed could be obtained by(https://www.biostars.org/p/13290/):
# Go to the UCSC table browser.
# Select desired species and assembly
# Select group: Genes and Gene Prediction Tracks
# Select track: RefSeq Genes(UCSC Genes, Ensembl, etc.)
# Select table: refGene
# Select region: genome (or you can test on a single chromosome or smaller region)
# Select output format: BED - browser extensible data
# Enter output file: hg19_Exons.bed(hg19_Genes.bed)
# Select file type returned: gzip compressed
# Hit the 'get output' button
# A second page of options relating to the BED file will appear.
# Under 'create one BED record per:'. Select 'Introns plus'(or Whole Gene)
# leave as 0 to get just the introns(or Whole Gene)
# Hit the 'get BED' option

# 2. Now you can extract intron coordinates for each gene using the merged/squashed exons. You should also be able to do this with BEDTools where you have one BED file containing your squashed exons and another containing the outer boundaries of each gene.
# 3. for chip-seq merging exons, we do not consider strand
# 4. get chromInfo.txt from chromInfo.sh
# 5. run this script

if [ $# -lt 4 ];
then
    echo "$0 gene.bed exon.bed size <limit>"
    echo "example: $0 gene.bed exon.bed 2000 hg19_chromInfo.txt"
    exit 1
fi

genes=$1
exon=$2
size=$3
limit=$4

## get promoter +- size kb, consider strand because of transcript strandness
awk -v s=$size 'BEGIN{OFS="\t"} {if ($6 == "+"){
n+=1
if ($2-s>=0)  {
    start=$2-s
}
else{
    start=0
}
print $1,start,$2+s,n,0,$6
}
else {
if ($3-s>=0)  {
    start=$3-s
}
else {
    start=0
}
print $1,start,$3+s,n,0,$6
}
}' $genes | cut -f 1,2,3 | sort -k1,1 -k2,2n | uniq | bedtools merge -i - > ${genes}.tmp
bedClip ${genes}.tmp $limit ${genes}.filter
mv ${genes}.filter ${genes}_promoter
rm ${genes}.tmp 

## exon and intron remove promoter regions
cut -f 1,2,3 $exon | sort -k 1,1 -k 2,2n - | uniq | bedtools merge -i - > ${exon}.merged

cut -f 1,2,3 $genes | sort -k 1,1 -k 2,2n - | uniq | bedtools merge -i - > ${genes}.merged

## remove exons from genes to get artificial introns 
bedtools subtract -a ${genes}.merged -b ${exon}.merged  > ${genes}_intron.tmp

## when merging meta introns, remove promoter 
bedtools subtract -a ${genes}_intron.tmp -b ${genes}_promoter | sort -k1,1 -k2,2n | uniq | bedtools merge -i - > ${genes}_intron

## remove promoters from exons
bedtools subtract -a ${exon}.merged -b ${genes}_promoter > ${genes}_exon
rm ${exon}.merged ${genes}_intron.tmp

## get intergenic regions, do not consider strand
bedtools complement -i ${genes}.merged -g $limit | bedtools merge -i - > ${genes}_intergenic.tmp

bedtools subtract -a ${genes}_intergenic.tmp -b ${genes}_promoter > ${genes}_intergenic

rm ${genes}_intergenic.tmp

#!/bin/bash 

#SBATCH -n 12  #Number of cores
#SBATCH -N 1   #node number

#SBATCH -t 3000  #Runtime in minutes 

#SBATCH -p general  #Partition to submit to 

#SBATCH --mem=10600  #Memory per node in MB (see also --mem-per-cpu)

#SBATCH -o S4-%A_%a.out
#SBATCH -J S4-%A_%a.err

source ../chilin_env/bin/activate
chilin simple -u qinq -s mm10 --threads 8 -i H3K4me3_S4 -o S4 -t /mnt/Storage/home/qinq/01_Projects/Programming/chilin/demo/GSM487548/SRR036173.fastq,/mnt/Storage/home/qinq/01_Projects/Programming/chilin/demo/GSM487550/SRR036175.fastq -p narrow -r histone --dont_remove

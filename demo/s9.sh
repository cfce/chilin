#!/bin/bash 

#SBATCH -n 12  #Number of cores
#SBATCH -N 1   #node number

#SBATCH -t 3000  #Runtime in minutes 

#SBATCH -p general  #Partition to submit to 

#SBATCH --mem=10600  #Memory per node in MB (see also --mem-per-cpu)

#SBATCH -o S9-%A_%a.out
#SBATCH -J S9-%A_%a.err

source ../chilin_env/bin/activate
chilin simple -u qinq -s hg38 --threads 8 -i AR_S9 -o S9 -t GSM980657/SRR531808.fastq,GSM980658/SRR531809.fastq -c GSM980659/SRR531810.fastq,GSM980661/SRR531812.fastq -p narrow -r tf --dont_remove

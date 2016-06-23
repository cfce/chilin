#!/bin/bash 

#SBATCH -n 12  #Number of cores
#SBATCH -N 1   #node number

#SBATCH -t 3000  #Runtime in minutes 

#SBATCH -p general  #Partition to submit to 

#SBATCH --mem=10600  #Memory per node in MB (see also --mem-per-cpu)

#SBATCH -o S7-%A_%a.out
#SBATCH -J S7-%A_%a.err

source ../chilin_env/bin/activate
chilin simple -u qinq -s hg38 --threads 8 -i FOXA1_S7 -o S7 -t GSM759662/SRR309194.fastq,GSM759661/SRR309193.fastq -c GSM759671/SRR309203.fastq -p narrow -r tf --dont_remove

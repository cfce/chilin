#!/bin/bash 

#SBATCH -n 12  #Number of cores
#SBATCH -N 1   #node number

#SBATCH -t 3000  #Runtime in minutes 

#SBATCH -p general  #Partition to submit to 

#SBATCH --mem=10600  #Memory per node in MB (see also --mem-per-cpu)

#SBATCH -o S10-%A_%a.out
#SBATCH -J S10-%A_%a.err

source ../chilin_env/bin/activate
chilin simple -u qinq -s mm10 --threads 8 -i RAD21_S10 -o S10 -t GSM672403/SRR099837_1.fastq,GSM672403/SRR099837_2.fastq  -p narrow -r tf --dont_remove --pe

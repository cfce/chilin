#!/bin/bash 

#SBATCH -n 12  #Number of cores
#SBATCH -N 1   #node number

#SBATCH -t 3000  #Runtime in minutes 

#SBATCH -p general  #Partition to submit to 

#SBATCH --mem=10600  #Memory per node in MB (see also --mem-per-cpu)

#SBATCH -o S11-%A_%a.out
#SBATCH -J S11-%A_%a.err

source ../chilin_env/bin/activate
chilin simple -u qinq -s mm10 --threads 8 -i RAG2_S11 -o S11 -t GSM530318/SRR039351.fastq -p narrow -r tf --dont_remove

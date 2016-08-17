# Download fastq-dump here http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software
# run the examples in the paper

if ! which fastq-dump.2.4.5 &>/dev/null
then
   echo 'need fastq-dump.2.4.5 in PATH'
   exit 1
fi
echo "Downloading..."
# File S1
fastq-dump.2.4.5 SRR038978 -O GSM392049
fastq-dump.2.4.5 SRR038979 -O GSM392050

# File S2
fastq-dump.2.4.5 SRR014982 -O GSM325898
fastq-dump.2.4.5 SRR014983 -O GSM325898
fastq-dump.2.4.5 SRR014984 -O GSM325898
fastq-dump.2.4.5 SRR014985 -O GSM325898

# File S3
wget -c http://hgdownload.cse.ucsc.edu/goldenPath/mm9/encodeDCC/wgEncodeLicrHistone/wgEncodeLicrHistoneBmarrowH3k4me3MAdult8wksC57bl6StdRawDataRep1.fastq.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenPath/mm9/encodeDCC/wgEncodeLicrHistone/wgEncodeLicrHistoneBmarrowH3k4me3MAdult8wksC57bl6StdRawDataRep2.fastq.gz

# File S4
fastq-dump.2.4.5 SRR036173 -O GSM487548	
fastq-dump.2.4.5 SRR036175 -O GSM487550

# File S5
wget -c http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHistone/wgEncodeBroadHistoneGm12878H2azStdRawDataRep1.fastq.gz
wget -c http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHistone/wgEncodeBroadHistoneGm12878H2azStdRawDataRep2.fastq.gz

 File S6
fastq-dump.2.4.5 SRR060525 -O GSM566158
fastq-dump.2.4.5 SRR060531 -O GSM566164
fastq-dump.2.4.5 SRR060535 -O GSM566168

# File S7
fastq-dump.2.4.5 SRR309194 -O GSM759662
fastq-dump.2.4.5 SRR309193 -O GSM759661
fastq-dump.2.4.5 SRR309203 -O GSM759671

# File S8
fastq-dump.2.4.5 SRR039652 -O GSM445803

# File S9
fastq-dump.2.4.5 SRR531808 -O GSM980657
fastq-dump.2.4.5 SRR531809 -O GSM980658
fastq-dump.2.4.5 SRR531810 -O GSM980659
fastq-dump.2.4.5 SRR531812 -O GSM980661

# File S10
fastq-dump.2.4.5 --split-3 SRR099837 -O GSM672403

# File S11
fastq-dump.2.4.5 SRR039351 -O GSM530318

# File S12
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR023/ERR023718/ERR023718.fastq.gz -O ERX009561.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/ERA295/ERA295486/fastq/s_4.fastq.gz -O ERX009560.fastq.gz

echo "running..."
for i in *bwa
do
  echo $i
  bash $i
done

# for bowtie aligner, should uncompress the .fastq.gz to .fastq
for i in *bowtie
do
  echo $i
  bash $i
done

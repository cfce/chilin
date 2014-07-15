"""
a sampling by reads ratio module, reads number (e.g. 4M) / total reads
"""



def sampling(orig, dest, rand, format, conf): # call fastq_sampling
    """
    prepare sampling fastq files for library contamination and fastqc
    rand: the number of random selected fastq reads
    use lh3's https://github.com/lh3/seqtk/ to sample fastq and fastq.gz
    """
    ## samtools sampling
    ## add judge condition
    return ShellCommand("""
                        count=$({tool} view -Sc {input[sam]})
                        ## judge mapped reads number less than sampling number
                        if [ $count -le {param[random_number]} ]
                        then
                            ln -f {input[sam]} {input[sam]}.{param[random_number]}
                            {tool} view -bS {input[sam]}.{param[random_number]} > {output[samp]}
                        else
                            sampling_pe_sam.py {input[sam]} {param[random_number]}
                            {tool} view -bS {input[sam]}.{param[random_number]} > {output[samp]}
                        fi
                        """,
                        tool = "samtools",
                        input={"sam": orig},
                        output={"samp": dest},
                        param={"random_number": rand},
                        name = "sampling bam")


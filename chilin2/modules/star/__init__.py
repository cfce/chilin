from samflow.command import ShellCommand
from samflow.workflow import attach_back

from chilin2.modules.bwa.qc import stat_bwa
from chilin2.modules.bwa.tex import tex_bwa

def star(workflow, conf):   # Mapping
    """
    Use star to map reads to genome,
    call __bwa_sam2bam to convert sam to bam
    :param workflow: samflow defined class
    :param conf: parsed config files
    :return: void
    """
    for target in conf.sample_targets:
        star = attach_back(workflow,
                          ShellCommand(
                              "{tool} --genomeDir {param[index]} --runThreadN {param[NUM_THREADS]} --readFilesIn {input[fastq]} --outFileNamePrefix {param[prefix]}",
                              tool = "STAR",
                              input = {"fastq": target + ".fastq"},
                              output = {"sam": target + "Aligned.out.sam"},
                              param = {"NUM_THREADS": conf.threads,
                                       "prefix": target,
                                       ## judge chosen species from basics section
                                       "index": conf.get_path(conf.get("basics", "species"), "genome_index")},
                              name = "star aln"))
        star.update(param = conf.items("bowtie"))

    _star_sam2bam(workflow, conf)

    ## QC part--NOTE keeping the bwa legacy code!
    stat_bwa(workflow, conf)
    if conf.long:
        tex_bwa(workflow, conf)


def _star_sam2bam(workflow, conf):  # SAM -> BAM
    """
    convert SAM to BAM and use mapping quality as cutoff
    :param workflow: samflow defined class
    :param conf: parsed config file
    :return: void
    """
    for target in conf.sample_targets:
        sam2bam= attach_back(workflow,
                    ShellCommand(
                        """
                        ln -s {input[sam]} {output[sam]}
                        {tool} view -q 255 -bt {param[genome]} {input[sam]} -o {output[bam]}
                        """,
                        tool="samtools",
                        input={"sam": target + "Aligned.out.sam"},
                        output={"bam":target + ".bam",
                                "sam": target + ".sam"},
                        param={"genome": conf.get(conf.get("basics", "species"), "chrom_len"),
                               },
                        name = "star sam2dam"))
        workflow.update(param=conf.items("sam2bam"))

        #From bwa/dc.py
        sam2bamnochrm = attach_back(workflow,  ## use mapping quality 1 defined by samtools official FAQ
                    ShellCommand(
                        """
                        grep -v chrM {param[chrom_bed]} > {output[nochrmbed]}
                        {tool} view -h -b -L {output[nochrmbed]} {input[bam]} > {output[nochrmbam]}
                        {tool} view -h {output[nochrmbam]}  > {output[nochrmsam]}
                        {tool} view -h {input[bam]}  > {output[usam]}
                        """,
                        tool="samtools",
                        input={"bam": target + ".bam"},
                        output={"nochrmbed": target + ".nochrM",
                                "nochrmbam": target + "_nochrM.bam",
                                "usam": target + "_u.sam", ## uniquely mapping sam for sampling
                                "nochrmsam": target + "_nochrM.sam"},
                        param={"tmp_bam": target + ".tmp.bam", "output_prefix": target,
                               "chrom_bed": conf.get(conf.get("basics", "species"), "chrom_bed"),
                               "mapq": 1,
                               "genome": conf.get(conf.get("basics", "species"), "chrom_len")},
                        name = "filtering mapping and convert")) # Use 5G memory as default
        sam2bamnochrm.update(param=conf.items("sam2bam"))

from samflow.command import ShellCommand
from samflow.workflow import attach_back

from chilin2.modules.bwa.qc import stat_bwa
from chilin2.modules.bwa.tex import tex_bwa


def bowtie(workflow, conf):   # Mapping
    """
    Use bowtie to map reads to genome,
    call __bwa_sam2bam to convert sam to bam
    :param workflow: samflow defined class
    :param conf: parsed config files
    :return: void
    """
    for target in conf.sample_targets:
        bowtie = attach_back(workflow,
                          ShellCommand(
                              "{tool} -p {param[NUM_THREADS]} -S -m 1 {param[index]} {input[fastq]} {output[sam]}",
                              tool = "bowtie",
                              input = {"fastq": target + ".fastq"},
                              output = {"sam": target + ".sam"},
                              param = {"NUM_THREADS": conf.threads,
                                       ## judge chosen species from basics section
                                       "index": conf.get_path(conf.get("basics", "species"), "genome_index")},
                              name = "bowtie aln"))
        bowtie.update(param = conf.items("bowtie"))

    _bowtie_sam2bam(workflow, conf)

    ## QC part--NOTE keeping the bwa legacy code!
    stat_bwa(workflow, conf)
    if conf.long:
        tex_bwa(workflow, conf)


def _bowtie_sam2bam(workflow, conf):  # SAM -> BAM
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
                        {tool} view -bt {param[genome]} {input[sam]} -o {output[bam]}
                        """,
                        tool="samtools",
                        input={"sam": target + ".sam"},
                        output={"bam":target + ".bam"},
                        param={"genome": conf.get(conf.get("basics", "species"), "chrom_len"),
                               },
                        name = "bowtie sam2dam"))
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

from samflow.command import ShellCommand
from samflow.workflow import attach_back

from chilin2.modules.bwa.qc import stat_bwa
from chilin2.modules.bwa.tex import tex_bwa


def bwa(workflow, conf):   # Mapping
    """
    Use BWA to map reads to genome,
    call __bwa_sam2bam to post-mapping filter
    :param workflow: samflow defined class
    :param conf: parsed config files
    :return: void
    """
    for target in conf.sample_targets:
        bwa = attach_back(workflow,
                          ShellCommand(
                              "{tool} aln -q 5 -n 0.04 -l 32 -k 2 -t {param[NUM_THREADS]} {param[index]} {input[fastq]} > {output[sai]}",
                              #"{tool} aln -q 5 -l 32 -k 2 -t {param[NUM_THREADS]} {param[index]} {input[fastq]} > {output[sai]}",
                              tool = "bwa",
                              input = {"fastq": target + ".fastq"},
                              output = {"sai": target + ".sai"},
                              param = {"NUM_THREADS": conf.threads,
                                       ## judge chosen species from basics section
                                       "index": conf.get_path(conf.get("basics", "species"), "genome_index")},
                              name = "bwa aln"))

        bwa_coordin = attach_back(workflow,
                                  ShellCommand(
                                      "{tool} samse {param[index]} {input[sai]} {input[fastq]} > {output[sam]}",
                                      tool = "bwa",
                                      input = {"fastq": target + ".fastq",
                                               "sai": target + ".sai"},
                                      output = {"sam": target + ".sam"},
                                      param = {"index": conf.get_path(conf.get("basics", "species"), "genome_index")},
                                      name = "bwa samse"))
        bwa.update(param = conf.items("bwa"))

    _bwa_sam2bam(workflow, conf)

    ## QC part
    stat_bwa(workflow, conf)
    if conf.long:
        tex_bwa(workflow, conf)


def _bwa_sam2bam(workflow, conf):  # SAM -> BAM
    """
    convert SAM to BAM and use mapping quality as cutoff
    :param workflow: samflow defined class
    :param conf: parsed config file
    :return: void
    """
    for target in conf.sample_targets:
        sam2bam = attach_back(workflow,  ## use mapping quality 1 defined by samtools official FAQ
                    ShellCommand(
                        """
                        {tool} view -bS -t {param[genome]} -q {param[mapq]} {input[sam]} > {param[tmp_bam]} && {tool} sort -m {param[max_mem]} {param[tmp_bam]} {param[output_prefix]}
                        """,
                        tool="samtools",
                        input={"sam": target + ".sam"},
                        output={"bam":target + ".bam"},
                        param={"tmp_bam": target + ".tmp.bam", "output_prefix": target,
                               "mapq": 1,
                               "genome": conf.get(conf.get("basics", "species"), "chrom_len"),
                               "max_mem": 4000000000},
                        name = "filtering mapping and convert")) # Use 5G memory as default
        sam2bam.update(param=conf.items("sam2bam"))

        ## filter chrM for chromosomoe information, generate uniquely
        ## aligned sam files (prepared for built-in sampling utility)
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
                        param={"chrom_bed": conf.get(conf.get("basics", "species"), "chrom_bed")},
                        name = "filtering chrM and convert to sam for sampling"))
        sam2bamnochrm.update(param=conf.items("sam2bam"))

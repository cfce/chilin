"""
Test whether library is mixed by other species
"""

from samflow.command import ShellCommand
from samflow.workflow import attach_back

from chilin2.modules.contamination.qc import stat_contamination
from chilin2.modules.contamination.tex import tex_contamination

def star(workflow, conf, target, output, index):   # Mapping
    """
    Use star to map reads to genome,
    """
    star = attach_back(workflow,
                      ShellCommand(
                          "{tool} --genomeDir {param[index]} --runThreadN {param[NUM_THREADS]} --readFilesIn {input[fastq]} --outFileNamePrefix {param[prefix]}",
                          tool = "STAR",
                          input = {"fastq": target + "_100k.fastq"},
                          output = {"sam": output},
                          param = {"NUM_THREADS": conf.threads,
                                   "prefix": target,
                                   ## judge chosen species from basics section
                                   "index": index},
                          name = "star aln"))
    star.update(param = conf.items("bowtie"))
    star.allow_dangling = True
    star.allow_fail = True
    return workflow


def bowtie(workflow, conf, target, output, index):   # Mapping
    """
    Use bowtie to map reads to genome,
    """
    bowtie = attach_back(workflow,
                      ShellCommand(
                          "{tool} -p {param[NUM_THREADS]} -S -m 1 {param[index]} {input[fastq]} {output[sam]}",
                          tool = "bowtie",
                          input = {"fastq": target + "_100k.fastq"},
                          output = {"sam": output},
                          param = {"NUM_THREADS": conf.threads,
                                   ## judge chosen species from basics section
                                   "index": index},
                                   name = "bowtie aln"))
    bowtie.update(param = conf.items("bowtie"))
    bowtie.allow_fail = True
    bowtie.allow_dangling = True
    return workflow

def bwa(workflow, conf, target, output, outsai, index):   # Mapping
    """
    Use BWA to map reads to genome,
    bwa before 0.6.0 support ab solid color space
    """

    if not conf.pe:
        bwa = attach_back(workflow,
                          ShellCommand(
                              "{tool} aln -q 5 -l 32 -k 2 -t {param[NUM_THREADS]} {param[index]} {input[fastq]} > {output[sai]}",
                              tool = "bwa",
                              input = {"fastq": target + "_100k.fastq"},
                              output = {"sai": outsai},
                              param = {"NUM_THREADS": conf.threads,
                                       ## judge chosen species from basics section
                                       "index": index},
                              name = "bwa aln"))
        bwa.update(param = conf.items("bwa"))
        out = attach_back(workflow,
                    ShellCommand(
                        "{tool} samse {param[index]} {input[sai]} {input[fastq]} > {output[sam]}",
                        tool = "bwa",
                        input = {"fastq": target + "_100k.fastq",
                                 "sai": outsai},
                        output = {"sam": output},
                        param = {"index": index},
                        name = "bwa samse"))
    else:
        bwa = attach_back(workflow,
                          ShellCommand(
                              "{tool} aln -q 5 -l 32 -k 2 -t {param[NUM_THREADS]} {param[index]} {input[fastq1]} > {output[sai1]} && {tool} aln -q 5 -l 32 -k 2 -t {param[NUM_THREADS]} {param[index]} {input[fastq2]} > {output[sai2]}",
                              tool = "bwa",
                              input = {"fastq1": target[0] + "_100k.fastq",
                                       "fastq2": target[1] + "_100k.fastq"},
                              output = {"sai1": outsai[0],
                                        "sai2": outsai[1]},
                              param = {"NUM_THREADS": conf.threads,
                                       ## judge chosen species from basics section
                                       "index": index},
                              name = "bwa aln"))
        bwa.update(param = conf.items("bwa"))

        out = attach_back(workflow,
                    ShellCommand(
                        "{tool} sampe {param[index]} {input[sai1]} {input[sai2]} {input[fastq1]} {input[fastq2]}> {output[sam]}",
                        tool = "bwa",
                        input = {"fastq1": target[0] + "_100k.fastq",
                                 "fastq2": target[1] + "_100k.fastq",
                                 "sai1": outsai[0],
                                 "sai2": outsai[1]},
                        output = {"sam": output},
                        param = {"index": index},
                        name = "bwa sampe"))
    bwa.allow_dangling = True
    bwa.allow_fail = True
    out.allow_dangling = True
    out.allow_fail = True


def contamination_check(workflow, conf):
    """
    bowtie mapping back to different species
    """
    if conf.items("contamination"):
        for target in conf.sample_targets:
            for species in dict(conf.items("contamination")):
                index = conf.get("contamination", species)
                if conf.mapper == "bwa":
                    output = target + species + ".sam"
                    if conf.pe:
                        outsai = [target + "pair1.sai", target + "pair2.sai"]
                        targets = [ target + "pair1", target + "pair2" ]
                    else:
                        outsai = target + ".sai"
                        targets = target
                    bwa(workflow, conf, targets, output, outsai, index)
                elif conf.mapper == "bowtie":
                    output = target + species + ".sam"
                    bowtie(workflow, conf, target, output, index)
                elif conf.mapper == "star":
                    output = target + species + "Aligned.out.sam"
                    star(workflow, conf, target, output, index)

                sam2bam = attach_back(workflow,  ## use mapping quality 1 defined by samtools official FAQ
                            ShellCommand(
                                """
                                {tool} view -bS -t {param[genome]} -q {param[mapq]} {input[sam]} > {param[tmp_bam]} && {tool} sort -m {param[max_mem]} {param[tmp_bam]} {param[output_prefix]}
                                """,
                                tool="samtools",
                                input={"sam": output},
                                output={"bam":target + species + ".bam"},
                                param={"tmp_bam": target + species + ".tmp.bam", "output_prefix": target + species,
                                       "mapq": 1,
                                       "genome": conf.get(conf.get("basics", "species"), "chrom_len"),
                                       "max_mem": 4000000000},
                                name = "filtering mapping and convert")) # Use 5G memory as default

                sam2bam.update(param=conf.items("sam2bam"))
                sam2bam.allow_dangling = True
                sam2bam.allow_fail = True

                rem = attach_back(workflow, ShellCommand(
                """
                {tool} view -Sc {input[sam]} > {output[total]}
                {tool} flagstat {input[bam]} > {output[stat]}
                """,
                tool = "samtools",
                input = {"bam": target + species + ".bam",
                         "sam": output},
                output = {"stat": target + species + "_mapped." + conf.mapper,
                          "total": target + species + "_total." + conf.mapper},
                name = "contamination calculation"))
                rem.allow_fail = True
                rem.allow_dangling = True

        ## QC part
        stat_contamination(workflow, conf)
        if conf.long:
            tex_contamination(workflow, conf)



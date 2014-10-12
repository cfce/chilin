"""
separate interface script heress
"""

import re
from samflow.command import ShellCommand
from samflow.workflow import attach_back
from chilin2.modules.config.helpers import make_link_command, sampling, json_load
import os
import json


def groom_sequencing_files(workflow, conf):  # the start of ChiLin
    """
    interface of the ChiLin input
    support fastq and fastq.gz input format,
    prepare Input with symbol links
    not allow failed, this is the first step ...
    """
    not_groomed = []
    for raw, target in conf.sample_pairs:
        if not conf.pe:
            if re.search(r"\.(fastq.gz|fq.gz)$", raw, re.I):
                attach_back(workflow, make_link_command(orig=raw, dest=target + ".fastq"))
                attach_back(workflow, sampling({"fastq": target + ".fastq"}, {"fastq_sample": target + "_100k.fastq"}, 100000, "fastq", conf))

            elif re.search(r"\.(fastq|fq)$", raw, re.I):
                attach_back(workflow, make_link_command(orig=os.path.abspath(raw), dest=target + ".fastq"))
                attach_back(workflow, sampling({"fastq": target + ".fastq"}, {"fastq_sample": target + "_100k.fastq"}, 100000, "fastq", conf))
            else:
           ##     print(raw, " is neither fastq nor bam file. Skip grooming.")
                not_groomed.append([raw, target])
        else:

            if all(map(lambda x: re.search(r"\.(fastq.gz|fq.gz)", x, re.I), raw)):

                attach_back(workflow, make_link_command(orig=raw[0], dest=target[0] + ".fastq"))
                attach_back(workflow, make_link_command(orig=raw[1], dest=target[1] + ".fastq"))

                attach_back(workflow, sampling({"fastq": target[0] + ".fastq"}, {"fastq_sample": target[0] + "_100k.fastq"}, 100000, "fastq", conf))
                attach_back(workflow, sampling({"fastq": target[1] + ".fastq"}, {"fastq_sample": target[1] + "_100k.fastq"}, 100000, "fastq", conf))

            elif all(map(lambda x: re.search(r"\.(fastq|fq)", x, re.I), raw)):
                attach_back(workflow, make_link_command(orig=os.path.abspath(raw[0]), dest=target[0] + ".fastq"))
                attach_back(workflow, make_link_command(orig=os.path.abspath(raw[1]), dest=target[1] + ".fastq"))
                attach_back(workflow, sampling({"fastq": target[0] + ".fastq"}, {"fastq_sample": target[0] + "_100k.fastq"}, 100000, "fastq", conf))
                attach_back(workflow, sampling({"fastq": target[1] + ".fastq"}, {"fastq_sample": target[1] + "_100k.fastq"}, 100000, "fastq", conf))
            else:
           ##     print(raw, " is neither fastq nor bam file. Skip grooming.")
                not_groomed.append([raw, target])

def sampling_bam(workflow, conf):   ## sampling to 4M
    """
    sampling bam files through macs2 and bedtools
    """
    for target in conf.sample_targets:
        ## sampling treat and control simultaneously
        ## sampling bam by macs2 and convert to bam by bedtools
        ## if total mapped reads < 4M, use original bam files link to *4000000.bam
        ## extract mapped reads number from json files
        ## use uniquely mapped reads sampling
        sampling_u = attach_back(workflow, sampling(target + "_u.sam", target + "_4000000.bam", 4000000, "sam", conf))
        sampling_u.allow_dangling = True
        sampling_u.allow_fail = True

        ## use encode version of 5M non chrM reads to evaluate
        if conf.frip:
            samp = attach_back(workflow, sampling(target + "_nochrM.sam", target + "_5000000_nochrM.bam", 5000000, "sam", conf))
            samp.allow_fail = True
            samp.allow_dangling = True
        else: ## default
            ## change FRiP computing with merged peaks as reference, no chrM as comparison
            samp = attach_back(workflow, sampling(target + "_nochrM.sam", target + "_4000000_nochrM.bam", 4000000, "sam", conf))
            samp.allow_fail = True
            samp.allow_dangling = True

    ## sampling merged control data to 4M control data for SPP and FRiP peaks calling
    ## change FRiP computing with merged peaks as reference, no chrM as comparison
#    if conf.control_raws:
#        attach_back(workflow, sampling(conf.prefix + "_control.bam", conf.prefix + "_control_4000000.bam", 4000000, "bam", conf))


## TODO: separate input and chip bam merge, because some of IP data may not need to be merged
def merge_bams(workflow, conf):   ## merge input and chip bam
    """
    input multiple input and multiple control to merge into one file separately
    :return:
    """
    # merge all treatments into one
    merge_bams_treat = ShellCommand(
        "{tool} merge {output[merged]} {param[bams]}",
        tool="samtools",
        input=[target + ".bam" for target in conf.treatment_targets],
        output={"merged": conf.prefix + "_treatment.bam"})
    merge_bams_treat.param = {"bams": " ".join(merge_bams_treat.input)}
    merge_bams_treat.allow_fail = True
    merge_bams_treat.allow_dangling = True

    if len(conf.treatment_targets) > 1:
        attach_back(workflow, merge_bams_treat)
    else:
        # when there's only one treatment sample, use copying instead of merging
        attach_back(workflow, make_link_command(merge_bams_treat.input[0], merge_bams_treat.output["merged"]))

    # merging step will be skipped if control sample does not exist
    # So be careful to check whether there are control samples before using `_control.bam`
    if len(conf.control_targets) > 1:
        merge_bams_control = merge_bams_treat.clone
        merge_bams_control.input = [target + ".bam" for target in conf.control_targets]
        merge_bams_control.output = {"merged": conf.prefix + "_control.bam"}
        merge_bams_control.param = {"bams": " ".join(merge_bams_control.input)}
        attach_back(workflow, merge_bams_control)
    elif len(conf.control_targets) == 1:
        attach_back(workflow, make_link_command(conf.control_targets[0] + ".bam", conf.prefix + "_control.bam"))



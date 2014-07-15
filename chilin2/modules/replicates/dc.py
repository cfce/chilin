#!/usr/bin/env python

from samflow.command import ShellCommand
from samflow.workflow import attach_back
from chilin2.modules.replicates.qc import stat_replicates


def replicates(workflow, conf):
    replicates_peaks_overlap(workflow, conf)
    replicates_bw_correlation(workflow, conf)

    ## QC parts
    stat_replicates(workflow, conf)

    

def replicates_peaks_overlap(workflow, conf):  # peaks bed from each replicate
    """
    :param workflow: class from samflow
    :param conf: external parsed config file
    :return: workflow through attach_back
    """
    for i in range(len(conf.treatment_targets)):
        for j in range(i+1, len(conf.treatment_targets)):
            replicates_overlap = attach_back(workflow,
                                             ShellCommand("{tool} -f {param[p]} -a {input[0]} -b {input[1]} | wc -l > {output}",
                                                          tool = "intersectBed",
                                                          input = [ conf.treatment_targets[i] + "_sort_peaks.narrowPeak" if conf.get("macs2", "type").lower() in ["both", "narrow"] else conf.treatment_targets[i] + "_b_sort_peaks.broadPeak",
                                                                    conf.treatment_targets[j] + "_sort_peaks.narrowPeak" if conf.get("macs2", "type").lower() in ["both", "narrow"] else conf.treatment_targets[j] + "_b_sort_peaks.broadPeak"],
                                                          output = conf.prefix + "_%s_%s.overlap" % (i, j),
                                                          param = {"p": 0.3},
                                                          name = "Replicates peaks overlap QC"))
    ## generate a barplot for meta distribution

    replicates_overlap.update(param = conf.items("replicates"))
    return workflow

def replicates_bw_correlation(workflow, conf):    ## correlation among different replicates
    """
    Use UCSC binary bigWigCorrelate to calculate reads density correlation
    collections in json files from qc
    :param workflow: samflow class
    :param conf: parsed config files
    :return: void
    """
    replicates_correlation = attach_back(workflow,
                                         ShellCommand(
                                             "{tool} {param[input_list]} > {output}",
                                             tool = "wigCorrelate",
                                             input = [ target + "_treat.bw" for target in conf.treatment_targets ],
                                             output = conf.prefix + ".cor",
                                             param = {"input_list": []},
                                             name = "correlation between bigwiggle"))
    replicates_correlation.update(param={"input_list": " ".join(replicates_correlation.input)})

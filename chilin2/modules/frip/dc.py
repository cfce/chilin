from samflow.command import ShellCommand
from samflow.workflow import attach_back

from chilin2.modules.frip.qc import stat_frip
from chilin2.modules.frip.tex import tex_frip

def FRiP(workflow, conf):  # FRiP
    """
    Fraction of Reads in Peaks regions at 4M reads level
    For example: 2 treat, 2 control
    modify: without down sampling read peaks calling, use merged peaks for comparison
    """
    ## use merged peaks for evaluation after removing chrM reads 
    for t in conf.sample_targets:
        if conf.frip: ## sampling 5M reads
            reads = t + "_5000000_nochrM.bam"
        else:
            reads = t + "_4000000_nochrM.bam"
        frip = attach_back(workflow,
                           ShellCommand("""
                                        fr=$(bedtools intersect -f {param[p]} -wa -u -abam {input[reads]} -b {input[peaks]} -bed | wc -l)
                                        total=$(samtools flagstat {input[reads]} | head -1 | cut -d" " -f1)
                                        echo $fr,$total > {output[frip]}
                                        """,
                                        tool="intersectBed",
                                        input={"reads": reads if conf.down else t+"_nochrM.bam", "peaks": conf.prefix + "_sort_peaks.narrowPeak" if conf.get("macs2", "type") in ["both", "narrow"] else conf.prefix + "_b_sort_peaks.broadPeak"},
                                        output={"frip": t + ".frip"},
                                        param={"p": "1E-9"},
                                        name="FRiP score"))
        ## in case that peaks calling on 4M reads may be very poor,
        ## no peaks generated, allow fail and dangling
        frip.allow_fail = True
        frip.allow_dangling = True
    frip.update(param=conf.items("bedtools"))
    
    ## QC part
    stat_frip(workflow, conf)
    if conf.long:
        tex_frip(workflow, conf)


def macs_4M(workflow, conf):
    """
    calculate peaks bed for FRiP 4M level
    ## in case that peaks calling on 4M reads may be very poor,
    ## no peaks generated, allow fail and dangling
    
    :param workflow:
    :param conf:
    :return:
    """
    if conf.get("tool", "macs2"):
        macs2_bin = conf.get("tool", "macs2")
    else:
        macs2_bin = "macs2"

    for target in conf.treatment_targets:
        ## DNase, H3K4, H2AZ, all acetyl marks, or TF
        if conf.get("macs2", "type").lower() in ["both", "narrow"]: ## for DNase, H3K4, H2AZ, all acetyl marks, or TF
            macs2_on_rep_narrow = attach_back(workflow,
                                              ShellCommand(
                                                  "{tool} callpeak -q {param[fdr]} --keep-dup {param[keep_dup]} --shiftsize={param[shiftsize]} --nomodel -g {param[species]} \
                                                  {param[treat_opt]} {param[control_opt]} -n {param[description]}",
                                                  tool=macs2_bin,
                                                  input={"treat": target + "_4000000.bam.bed"},
                                                  output={"peaks": target + "_4M_peaks.narrowPeak",
                                                          "peaks_xls": target + "_4M_peaks.xls"},
                                                  param={"description": target + "_4M", "keep_dup": 1, "shiftsize": 73, "species": "hs", "fdr":0.01},
                                                  name="macs2_callpeak_rep"))
            macs2_on_rep_narrow.param["treat_opt"] = " -t " + macs2_on_rep_narrow.input["treat"]
            macs2_on_rep_narrow.allow_fail = True
            macs2_on_rep_narrow.allow_dangling = True


            sort = attach_back(workflow,
                               ShellCommand(
                                   "{tool} -r -g -k 9 {input} > {output}",
                                   tool = "sort",
                                   input = target + "_4M_peaks.narrowPeak",
                                   output = target + "_4M_sort_peaks.narrowPeak"))
            sort.allow_fail = True
            sort.allow_dangling = True

            # control option is skipped if control samples does not exist
            if len(conf.control_targets) >= 1:
                macs2_on_rep_narrow.input["control"] =  conf.prefix + "_control_4000000.bam.bed"
                macs2_on_rep_narrow.param["control_opt"] = "-c " + macs2_on_rep_narrow.input["control"]

            else:
                macs2_on_rep_narrow.param["control_opt"] = ""
            macs2_on_rep_narrow.update(param=conf.items("macs2"))


        if conf.get("macs2", "type").lower() in ["both", "broad"]:  # K9, K36, K79 and K27 methylation, both for chromatin regulator, all other histone marks
            macs2_on_rep_broad = attach_back(workflow,
                                             ShellCommand(
                                                 "{tool} callpeak -q {param[fdr]} {param[treat_opt]} {param[control_opt]} --keep-dup {param[keep_dup]} --broad -g {param[species]} -n {param[description]}",
                                                 tool=macs2_bin,
                                                 input={"treat": target + "_4000000.bam.bed"},
                                                 output={"peaks": target + "_4M_peaks.broadPeak"},
                                                 param={"description": target + "_4M",
                                                        "species": "hs",
                                                        "fdr": 0.01},
                                                 name=" macs2 broad peaks"))
            macs2_on_rep_broad.param["treat_opt"] = " -t " + macs2_on_rep_broad.input["treat"]
            macs2_on_rep_broad.update(param=conf.items("macs2"))

            macs2_on_rep_broad.allow_fail = True
            macs2_on_rep_broad.allow_dangling = True


            if len(conf.control_targets) >= 1:
                macs2_on_rep_broad.input["control"] =  conf.prefix + "_control_4000000.bam.bed" ## .bam.tmp is the sampling bed files
                macs2_on_rep_broad.param["control_opt"] = "-c " + macs2_on_rep_broad.input["control"]
            else:
                macs2_on_rep_broad.param["control_opt"] = ""
            macs2_on_rep_broad.update(param=conf.items("macs2"))

            sort = attach_back(workflow,
                               ShellCommand(
                                   "{tool} -r -g -k 9 {input} > {output}",
                                   tool = "sort",
                                   input = target + "_4M_peaks.broadPeak",
                                   output = target + "_4M_sort_peaks.broadPeak",
                                   name = "sort broad peaks"))

            sort.allow_fail = True
            sort.allow_dangling = True

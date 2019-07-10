"""
Keep with a stable version of MACS
"""

import json
import sqlite3
import math

from samflow.command import ShellCommand
from samflow.workflow import attach_back

from chilin2.modules.config.helpers import JinjaTemplateCommand, template_dump, r_exec, json_load, json_dump, underline_to_space
from chilin2.modules.config.helpers import make_link_command
from chilin2.modules.macs.qc import stat_macs2, stat_macs2_on_rep


def macs2(workflow, conf):
    if conf.get("tool", "macs2"):
        macs2_bin = conf.get("tool", "macs2")
    else:
        macs2_bin = "macs2"

    format = "-f BAMPE" if conf.pe else " "

    if conf.get("macs2", "type").lower() in ["both", "narrow"]: ## for DNase, H3K4, H2AZ, all acetyl marks, or TF
        macs2_on_merged_narrow = attach_back(workflow, ShellCommand(
            """
            {tool} callpeak --SPMR -B -q {param[fdr]} --keep-dup {param[keep_dup]} -g {param[species]} {param[format]} {param[treat_opt]} {param[control_opt]} -n {param[description]} && cut -f1,2,3,4,9 {output[peaks]} > {output[bedtmp]}
            ## remove weird path characters
            cp {output[peaks]} {output[peakstmp]}
            cp {output[summits]} {output[summitstmp]}
            awk \'{{OFS="\\t";n+=1;$4="peak"n;print $0}}\' {output[peakstmp]} > {output[peaks]}
            awk \'{{OFS="\\t";n+=1;$4=n;print $1,$2,$3,"peak"$4,$5}}\' {output[bedtmp]} > {output[bed]}
            awk \'{{OFS="\\t";n+=1;$4=n;print $1,$2,$3,"peak"$4,$5}}\' {output[summitstmp]} > {output[summits]}
            """,
            tool=macs2_bin,
            input={"treat": conf.prefix + "_treatment.bam"},
            output={"peaks": conf.prefix + "_peaks.narrowPeak",
                    "peakstmp": conf.prefix + "_peaks.narrowPeak.tmp",
                    "bed": conf.prefix + "_peaks.bed",
                    "bedtmp": conf.prefix + "_peaks.bed.tmp",
                    "summits": conf.prefix + "_summits.bed",
                    "summitstmp": conf.prefix + "_summits.bed.tmp",
                    "treat_bdg": conf.prefix + "_treat_pileup.bdg",
                    "peaks_xls": conf.prefix + "_peaks.xls",
                    "control_bdg": conf.prefix + "_control_lambda.bdg"},
            param={"description": conf.prefix,
                   "keep_dup": 1,
                   "format": format,
                   "extsize": 73 * 2, # extsize=2*shiftsize
                   "fdr": 0.01,
                   "species": "hs"},
            name="macs2_callpeak_merged"))

        macs2_on_merged_narrow.param["treat_opt"] = "-t " + macs2_on_merged_narrow.input["treat"]
        # control option is skipped if control samples does not exist
        if len(conf.control_targets) >= 1:
            macs2_on_merged_narrow.input["control"] = conf.prefix + "_control.bam"
            macs2_on_merged_narrow.param["control_opt"] = "-c " + macs2_on_merged_narrow.input["control"]
        else:
            macs2_on_merged_narrow.param["control_opt"] = ""

        macs2_on_merged_narrow.update(param=conf.items("macs2"))

        sort = attach_back(workflow,
                           ShellCommand(
                               """{tool} -r -g -k 9 {input[peaks]} > {output[p]}
                               {tool} -r -g -k 5 {input[summits]} > {output[s]}
                               """,
                               tool = "sort",
                               input = {"peaks": conf.prefix + "_peaks.narrowPeak",
                                        "summits": conf.prefix + "_summits.bed"},
                               output = {"p":conf.prefix + "_sort_peaks.narrowPeak",
                                         "s":conf.prefix + "_sort_summits.bed"},
                               name = "sort peaks"))
        macs2_on_merged_narrow.allow_fail = True
        macs2_on_merged_narrow.allow_dangling = True
        sort.allow_fail = True
        sort.allow_dangling = True

    if conf.get("macs2", "type").lower() in ["both", "broad"]:  # K9, K36, K79 and K27 methylation, both for chromatin regulator, all other histone marks
        macs2_on_merged_broad = attach_back(workflow,
                                            ShellCommand(
                                                """
                                                {tool} callpeak --SPMR -B -q {param[fdr]} {param[treat_opt]} {param[control_opt]} --keep-dup {param[keep_dup]} --extsize={param[extsize]} --nomodel --broad --broad-cutoff {param[fdr]} -g {param[species]} {param[format]} -n {param[description]} && cut -f1,2,3,4,9 {output[peaks]} > {output[bedtmp]}
                                                ## remove weird path characters
                                                cp {output[peaks]} {output[peakstmp]}
                                                awk \'{{OFS="\\t";n+=1;$4="peak"n;print $0}}\' {output[peakstmp]} > {output[peaks]}
                                                awk \'{{OFS="\\t";n+=1;$4=n;print $1,$2,$3,"peak"$4,$5}}\' {output[bedtmp]} > {output[bed]}
                                                """,
                                                tool=macs2_bin,
                                                input = {"treat": conf.prefix + "_treatment.bam"},
                                                output = {"peaks": conf.prefix + "_b_peaks.broadPeak",
                                                          "peakstmp": conf.prefix + "_b_peaks.broadPeak.tmp",
                                                          "bed": conf.prefix + "_b_peaks.bed",
                                                          "bedtmp": conf.prefix + "_b_peaks.bed.tmp",
                                                          "treat_bdg": conf.prefix + "_b_treat_pileup.bdg",
                                                          "peaks_xls": conf.prefix + "_b_peaks.xls",
                                                          "control_bdg": conf.prefix + "_b_control_lambda.bdg"},
                                                param = {"description": conf.prefix + "_b",
                                                         "species": "hs",
                                                         "format": format,
                                                         "keep_dup": 1, "fdr": 0.01},
                                                name = "broad peaks calling"))
        macs2_on_merged_broad.param["treat_opt"] = " -t " + macs2_on_merged_broad.input["treat"]
        macs2_on_merged_broad.allow_fail = True
        macs2_on_merged_broad.allow_dangling = True

        if len(conf.control_targets) >= 1:
            macs2_on_merged_broad.input["control"] = conf.prefix + "_control.bam"
            macs2_on_merged_broad.param["control_opt"] = "-c " + macs2_on_merged_broad.input["control"]
        else:
            macs2_on_merged_broad.param["control_opt"] = ""
        macs2_on_merged_broad.update(param=conf.items("macs2"))

        sort = attach_back(workflow,
                    ShellCommand(
                        """{tool} -r -g -k 9 {input[peaks]} > {output[p]}
                        """,
                        tool = "sort",
                        input = {"peaks": conf.prefix + "_b_peaks.broadPeak"},
                        output = {"p": conf.prefix + "_b_sort_peaks.broadPeak"},
                        name = "sort peaks files"))
        sort.allow_dangling = True
        sort.allow_fail = True

    # For bedGraphToBigwiggle bugs, we need to remove coordinates over-border coordinates
    # As _control_lambda.bdg always exist. There are no need to check whether there are control samples.

    if conf.get("macs2", "type").lower() in ["both", "broad"]:
        cont_bdg = conf.prefix + "_b_control_lambda.bdg"
        treat_bdg = conf.prefix + "_b_treat_pileup.bdg"
    if conf.get("macs2", "type").lower() in ["both", "narrow"]:
        cont_bdg = conf.prefix + "_control_lambda.bdg"
        treat_bdg = conf.prefix + "_treat_pileup.bdg"
    import os
    bdg_trim_control = attach_back(workflow,
                                   ShellCommand(
                                       '''
                                       {tool} intersect -a {input[bdg]} -b {param[chrom_bed]} -wa -f 1.00 | bedtools sort -i - > {output}
                                       ''',
                                       tool="bedtools",
                                       input={"bdg": cont_bdg},
                                       param = {"chrom_bed": os.path.join(conf.target_dir, "chrom.bed")},
                                       output=cont_bdg+".tmp",
                                       name="bedGraph filtering control"))
    bdg_trim_control.fail = True
    bdg_trim_control.allow_dangling = True

    bdg_trim_treat = bdg_trim_control.clone
    bdg_trim_treat.input["bdg"] = treat_bdg
    bdg_trim_treat.output = treat_bdg + ".tmp"
    bdg_trim_treat.fail = True
    bdg_trim_treat.allow_dangling = True

    attach_back(workflow, bdg_trim_treat)

    bdg2bw_treat = attach_back(workflow,
                               ShellCommand(
                                   "{tool} {input[bdg]} {input[chrom_len]} {output[bw]}",
                                   tool="bedGraphToBigWig",
                                   input={"bdg": cont_bdg+".tmp",
                                          "chrom_len": conf.get_path(conf.get("basics", "species"), "chrom_len")},
                                   output={"bw": conf.prefix + "_control.bw"},
                                   name="bdg_to_bw control"))
    ## in case broad peaks failed
    bdg2bw_treat.allow_fail = True
    bdg2bw_treat.allow_dangling = True

    # prototype used here to do the similar thing on treatment file
    bdg2bw_control = bdg2bw_treat.clone
    bdg2bw_control.input["bdg"] = treat_bdg+".tmp"
    bdg2bw_control.output["bw"] = conf.prefix + "_treat.bw"

    ## in case broad peaks failed
    bdg2bw_control.allow_fail = True
    bdg2bw_control.allow_dangling = True

    attach_back(workflow, bdg2bw_control)

    stat_macs2(workflow, conf)

def macs2_rep(workflow, conf):
    # Though macs command already exists, I choose not to use prototype here
    # Because the prototype definition and usage might be far from each other, making codes not readable
    if conf.get("tool", "macs2"):
        macs2_bin = conf.get("tool", "macs2")
    else:
        macs2_bin = "macs2"

    format = " -f BAMPE " if conf.pe  else " "

    for target in conf.treatment_targets:
        ## DNase, H3K4, H2AZ, all acetyl marks, or TF
        if conf.get("macs2", "type").lower() in ["both", "narrow"]: ## for DNase, H3K4, H2AZ, all acetyl marks, or TF
            macs2_on_rep_narrow = attach_back(workflow,
                                              ShellCommand(
                                                  """
                                                  {tool} callpeak --SPMR -B -q {param[fdr]} --keep-dup {param[keep_dup]} --extsize={param[extsize]} --nomodel -g {param[species]} {param[format]} {param[treat_opt]} {param[control_opt]} -n {param[description]} && cut -f1,2,3,4,9 {output[peaks]} > {output[bedtmp]}
                                                  ## remove weird path characters
                                                  cp {output[peaks]} {output[peakstmp]}
                                                  cp {output[summits]} {output[summitstmp]}
                                                  awk \'{{OFS="\\t";n+=1;$4="peak"n;print $0}}\' {output[peakstmp]} > {output[peaks]}
                                                  awk \'{{OFS="\\t";n+=1;$4=n;print $1,$2,$3,"peak"$4,$5}}\' {output[bedtmp]} > {output[bed]}
                                                  awk \'{{OFS="\\t";n+=1;$4=n;print $1,$2,$3,"peak"$4,$5}}\' {output[summitstmp]} > {output[summits]}
                                                  """,
                                                  tool=macs2_bin,
                                                  input={"treat": target + ".bam"},
                                                  output={"peaks": target + "_peaks.narrowPeak",
                                                          "peakstmp": target + "_peaks.narrowPeak.tmp",
                                                          "summits": target + "_summits.bed",
                                                          "summitstmp": target + "_summits.bed.tmp",
                                                          "bed": target + "_peaks.bed",
                                                          "bedtmp": target + "_peaks.bed.tmp",
                                                          "treat_bdg": target + "_treat_pileup.bdg",
                                                          "peaks_xls": target + "_peaks.xls",
                                                          "control_bdg": target + "_control_lambda.bdg"},
                                                  param={"description": target, "keep_dup": 1, "extsize": 73*2, "species": "hs", "fdr":0.01, "format": format},
                                                  name="macs2_callpeak_rep"))
            macs2_on_rep_narrow.param["treat_opt"] = "-t " + macs2_on_rep_narrow.input["treat"]

            sort = attach_back(workflow,
                               ShellCommand(
                                   "{tool} -r -g -k 9 {input} > {output}",
                                   tool = "sort",
                                   input = target + "_peaks.narrowPeak",
                                   output = target + "_sort_peaks.narrowPeak"))

            # control option is skipped if control samples does not exist
            if len(conf.control_targets) >= 1:
                macs2_on_rep_narrow.input["control"] = conf.prefix + "_control.bam"
                macs2_on_rep_narrow.param["control_opt"] = "-c " + macs2_on_rep_narrow.input["control"]

            else:
                macs2_on_rep_narrow.param["control_opt"] = ""
            macs2_on_rep_narrow.update(param=conf.items("macs2"))
            macs2_on_rep_narrow.allow_dangling = True
            macs2_on_rep_narrow.allow_fail = True


        if conf.get("macs2", "type").lower() in ["both", "broad"]:  # K9, K36, K79 and K27 methylation, both for chromatin regulator, all other histone marks
            macs2_on_rep_broad = attach_back(workflow,
                                             ShellCommand(
                                                 """
                                                 {tool} callpeak --SPMR -B -q {param[fdr]} {param[treat_opt]} {param[control_opt]} --keep-dup {param[keep_dup]} --broad --broad-cutoff {param[fdr]} -g {param[species]} {param[format]} -n {param[description]} && cut -f1,2,3,4,9 {output[peaks]} > {output[bedtmp]}
                                                 ## remove weird path characters
                                                 cp {output[peaks]} {output[peakstmp]}
                                                 awk \'{{OFS="\\t";n+=1;$4="peak"n;print $0}}\' {output[peakstmp]} > {output[peaks]}
                                                 awk \'{{OFS="\\t";n+=1;$4=n;print $1,$2,$3,"peak"$4,$5}}\' {output[bedtmp]} > {output[bed]}
                                                 """,
                                                 tool=macs2_bin,
                                                 input = {"treat": target + ".bam"},
                                                 output = {"peaks": target + "_b_peaks.broadPeak",
                                                           "peakstmp": target + "_b_peaks.broadPeak.tmp",
                                                           "bed": target + "_b_peaks.bed",
                                                           "bedtmp": target + "_b_peaks.bed.tmp",
                                                           "treat_bdg": target + "_b_treat_pileup.bdg",
                                                           "peaks_xls": target + "_b_peaks.xls",
                                                           "control_bdg": target + "_b_control_lambda.bdg"},
                                                 param = {"description": target + "_b",
                                                          "species": "hs",
                                                          "format": format,
                                                          "fdr": 0.01},
                                                 name = " macs2 broad peaks"))
            macs2_on_rep_broad.param["treat_opt"] = " -t " + macs2_on_rep_broad.input["treat"]
            macs2_on_rep_broad.update(param=conf.items("macs2"))
            macs2_on_rep_broad.allow_dangling = True
            macs2_on_rep_broad.allow_fail=True

            if len(conf.control_targets) >= 1:
                macs2_on_rep_broad.input["control"] = conf.prefix + "_control.bam"
                macs2_on_rep_broad.param["control_opt"] = "-c " + macs2_on_rep_broad.input["control"]
            else:
                macs2_on_rep_broad.param["control_opt"] = ""
            ## some broad peaks cannot be called
            macs2_on_rep_broad.update(param=conf.items("macs2"))

            sort = attach_back(workflow,
                               ShellCommand(
                                   "{tool} -r -g -k 9 {input} > {output}",
                                   tool = "sort",
                                   input = target + "_b_peaks.broadPeak",
                                   output = target + "_b_sort_peaks.broadPeak",
                                   name = "sort broad peaks"))
            sort.allow_dangling = True
            sort.allow_fail=True

        ## For bedGraphToBigwiggle bugs, we need to remove coordinates outlier
        if conf.get("macs2", "type").lower() in ["both", "broad"]:
            cont_bdg = target + "_b_control_lambda.bdg"
            treat_bdg = target + "_b_treat_pileup.bdg"
        if conf.get("macs2", "type").lower() in ["both", "narrow"]:
            cont_bdg = target + "_control_lambda.bdg"
            treat_bdg = target + "_treat_pileup.bdg"
        import os
        bdg_trim_controlrep = attach_back(workflow,
                                          ShellCommand(
                                              '{tool} intersect -a {input} -b {param[chrom_bed]} -wa -f 1.00 | bedtools sort -i - > {output}',
                                              tool="bedtools",
                                              input=cont_bdg,
                                              output=cont_bdg + ".tmp",
                                              param={"chrom_bed": os.path.join(conf.target_dir, "chrom.bed")},
                                              name="bedGraph control replicate filtering"))

        bdg_trim_controlrep.allow_dangling = True
        bdg_trim_controlrep.allow_fail=True

        bdg_trim_treatrep = bdg_trim_controlrep.clone
        bdg_trim_treatrep.input = treat_bdg
        bdg_trim_treatrep.output = treat_bdg + ".tmp"

        bdg_trim_treatrep.allow_dangling = True
        bdg_trim_treatrep.allow_fail=True

        attach_back(workflow, bdg_trim_treatrep)

        bdg2bw_treatrep = attach_back(workflow,
                                      ShellCommand(
                                          "{tool} {input} {param[chrom_len]} {output}",
                                          tool="bedGraphToBigWig",
                                          input=treat_bdg+".tmp",
                                          output=target + "_treat.bw",
                                          param={"chrom_len": conf.get_path(conf.get("basics", "species"), "chrom_len")},
                                          name="bdg_to_bw treat"))
        ## in case broad peaks calling failed
        bdg2bw_treatrep.allow_dangling = True
        bdg2bw_treatrep.allow_fail=True

        # prototype used here to do the similar thing on treatment file
        bdg2bw_controlrep = bdg2bw_treatrep.clone
        bdg2bw_controlrep.input = cont_bdg + ".tmp"
        bdg2bw_controlrep.output = target + "_control.bw"
        attach_back(workflow, bdg2bw_controlrep)

        ## in case broad peaks calling failed
        bdg2bw_controlrep.allow_dangling = True
        bdg2bw_controlrep.allow_fail=True

    stat_macs2_on_rep(workflow, conf)

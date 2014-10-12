"""
This module is using 4M reads or total reads (specified by --total)
to estimate fragment size through latest macs2
"""
import json
import sqlite3
import math

from samflow.command import ShellCommand
from samflow.command import PythonCommand
from samflow.workflow import attach_back

from chilin2.modules.config.helpers import JinjaTemplateCommand, template_dump, r_exec, json_load, json_dump, underline_to_space
from chilin2.modules.config.helpers import make_link_command
from chilin2.modules.macs2_fragment.qc import stat_frag_std

def fragment(workflow, conf):
    ## this is done after FRiP
    if conf.get("tool", "macs2"):
       macs2_bin = conf.get("tool", "macs2")
    else:
        macs2_bin = "macs2"
    for target in conf.treatment_targets:
        fragment_size = attach_back(workflow, ShellCommand(
            "{tool} predictd -i {input[bam]} --rfile {param[prefix]} -g {param[species]}",
            tool = macs2_bin,
            input = {"bam": target + ".bam"},
            output = {"R": target + "_model.R"},
            param = {"prefix": target + "_model.R",
                     "species": 'hs'}))
        fragment_size.update(param = conf.items("macs2"))
        ## except too few peaks for modeling
        fragment_size.allow_fail = True
        fragment_size.allow_dangling = True

    ## extract standard deviation from MACS2 model.R,
    ## use m, p, and pileup value for standard deviation; mean fragment size is provided (choose the one with highest correlation)
    frag_qc = attach_back(workflow, PythonCommand(
        stat_frag_std,
        input = {"r": [target + "_model.R" for target in conf.treatment_targets]},
        output = {"json": conf.json_prefix + "_frag.json", "r": [ target + "_frag_sd.R" for target in conf.treatment_targets ]},
        param = {"samples": conf.treatment_bases,
                 "frag_tool": "BAMSE"},
        name = "macs2 model R script parser"))
    frag_qc.allow_fail = True
    frag_qc.allow_dangling = True

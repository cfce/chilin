import json
import os
import re
import sqlite3

from chilin2.modules.config.helpers import JinjaTemplateCommand, template_dump, r_exec, json_load, underline_to_space, count_in_million, decimal_to_latex_percent, json_dump
from samflow.command import PythonCommand
from samflow.command import ShellCommand
from samflow.workflow import attach_back
from pkg_resources import resource_filename

from chilin2.modules.phantompeak.qc import stat_phan
from chilin2.modules.phantompeak.tex import tex_phan


def Phan(workflow, conf): # NSC, RSC, Qtag
    """
    for calculating NSC, RSC score at 4M level
    http://code.google.com/p/phantompeakqualtools/
    (1) Determine strand cross-correlation peak / predominant fragment length OR print out quality measures
        Rscript run_spp.R -c=<tagAlign/BAMfile> -savp -out=<outFile>
    """
    # peaks calling by SPP needs control, for phantomqc, we do both treat and control independently

    
    for t in conf.sample_targets:
        if conf.down: ## default, this option
            ibam = t + "_4000000.bam"
        # elif conf.unsc: ## --total --unsc   
        #     ibam = t + "_rmdup.bam"
        else: ## --total
            ibam = t + ".bam"
        attach_back(workflow,
                    ShellCommand("{tool} {param[script]} -c={input[chip]} -rf -savp -out={output[spp]} -odir={param[dir]}",
                                 tool = "Rscript",
                                 input = {"chip": ibam},
                                 output = {"spp": t + ".spp", "pdf": t+"_4000000.pdf" if conf.down else t+".pdf"},
                                 param = {"script": conf.get("tool", "spp"),
                                          "dir": os.path.dirname(t + ".spp")},
                                 name = "SPP"))

    stat_phan(workflow, conf)
    if conf.long:
        tex_phan(workflow, conf)


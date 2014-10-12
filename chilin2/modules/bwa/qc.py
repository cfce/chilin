import json
import os
import re
import sqlite3

from chilin2.modules.config.helpers import JinjaTemplateCommand, template_dump, r_exec, json_load, underline_to_space, count_in_million, decimal_to_latex_percent, json_dump
from samflow.command import PythonCommand
from samflow.command import ShellCommand
from samflow.workflow import attach_back
from pkg_resources import resource_filename
from chilin2.modules.dbaccessor.caller import db_bwa


def stat_bwa(workflow, conf): ## use samtools to parse mappable reads from bwa
    """
    bam files are filtered by samtools -q 1, so mapped reads are considered to be unique
    """
    for t in conf.sample_targets:
        stat = attach_back(workflow, ShellCommand(
        """
        {tool} view -Sc {input[sam]} > {output[total]}
        {tool} flagstat {input[bam]} > {output[stat]}
        """,
        tool = "samtools",
        input = {"bam": t + ".bam",
                 "sam": t + ".sam"},
        output = {"stat": t + "_mapped.bwa",
                  "total": t + "_total.bwa"}))
        stat.allow_fail = True
        stat.allow_dangling = True
    collect = attach_back(workflow, PythonCommand(json_bwa,
        input={"bwa_mapped": [ t + "_mapped.bwa" for t in conf.sample_targets ],
               "bwa_total": [ t + "_total.bwa" for t in conf.sample_targets ]},
        output={"json": conf.json_prefix+"_map.json"},
        param={"sample":conf.sample_bases},
        name="bwa qc"))
    collect.allow_dangling = True
    collect.allow_fail = True

    if conf.long:
        long_collect = attach_back(workflow, PythonCommand(bwa_figures,
                                            input = {"dbaccessor": resource_filename("chilin2.modules.dbaccessor", "ChiLinQC.db"),
                                                     "json": conf.json_prefix + "_map.json",
                                                     "template": resource_filename("chilin2.modules.summary", "R_culmulative_plot.R")},


                                            output = {"pdf": conf.prefix + "_bwa_compare.pdf", "R": conf.prefix+"_bwa_compare.R"},
                                            param = {"sample": conf.sample_bases}))
        long_collect.allow_fail = True
        long_collect.allow_fail = True


def json_bwa(input={}, output={}, param={}): ## convert values to json files
    """
    input samtools flagstat standard output
    output json files
    kwargs for matching replicates order
    keep one value for each json for easier loading to html/pdf template
    example:
    3815725 + 0 in total (QC-passed reads + QC-failed reads)
    0 + 0 duplicates
    3815723 + 0 mapped (100.00%:-nan%)
    """
    json_dict = {"stat": {}, "input": input, "output": output, "param": param}

    for mapped, total, sam in zip(input["bwa_mapped"], input["bwa_total"], param["sample"]):
        inft = open(total, 'rU')
        infm = open(mapped, 'rU')
        json_dict["stat"][sam] = {}
        json_dict["stat"][sam]["mapped"] = int(infm.readlines()[2].split()[0])
        json_dict["stat"][sam]["total"] = int(inft.readlines()[0].strip())
        inft.close()
        infm.close()
    json_dump(json_dict)


def bwa_figures(input={"dbaccessor": "", "json": ""}, output = {"figure": ""}, param = {}):
    historyData = db_bwa(input["dbaccessor"])
    json_dict = json_load(input["json"])
    mappable_rates = [float(json_dict["stat"][i]["mapped"])/json_dict["stat"][i]['total'] for i in param["sample"]]
    mappable_rate_R = JinjaTemplateCommand(template=input["template"],
                                           param={'historic_data': historyData,
                                                  'current_data': mappable_rates,
                                                  'ids': param['sample'],
                                                  'cutoff': 0.5,
                                                  'section': 'Unique mapped rate',
                                                  "need_smooth_curve": True,
                                                  "render_dump": output["R"],
                                                  "pdf": output["pdf"]})
    template_dump(mappable_rate_R)
    r_exec(mappable_rate_R)

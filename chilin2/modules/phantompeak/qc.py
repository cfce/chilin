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

def stat_phan(workflow, conf):
    """
    collect NSC/RSC/Qtag and cross correlation figure
    """
    attach_back(workflow, PythonCommand(
        json_phan,
        input = {"spp": [t + ".spp" for t in conf.sample_targets]},
        output = {"json":conf.json_prefix + "_phan.json"},
        param = {"sample": conf.sample_bases}))

def json_phan(input = {"spp": ""}, output = {"json": ""}, param = {"sample": ""}):
    """ fragment size keep the maximus positive one
    """
    json_dict = {"stat": {}, "input": input, "output": output, "param": param}
    frag = 0
    for i, s in zip(input["spp"], param["sample"]):
        json_dict["stat"][s] = {}
        f = open(i)
        content = f.read().strip().split()
        f.close()
        json_dict["stat"][s]["NSC"] = content[8]
        json_dict["stat"][s]["RSC"] = content[9]
        for i in content[2].split(","):
            if i >= 0: 
                frag = i
                break

        json_dict["stat"][s]["frag"] = frag ## pick the most correlated ones
        json_dict["stat"][s]["Qtag"] = content[10]
    json_dump(json_dict)

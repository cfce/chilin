import json
import os
import re
import sqlite3

from chilin2.modules.config.helpers import JinjaTemplateCommand, template_dump, r_exec, json_load, underline_to_space, count_in_million, decimal_to_latex_percent, json_dump
from samflow.command import PythonCommand
from samflow.command import ShellCommand
from samflow.workflow import attach_back
from pkg_resources import resource_filename

def tex_phan(workflow, conf):
    figures = []
    for t in conf.sample_targets:
        if conf.down:
            figures.append(t + "_4000000.pdf")
        else:
            figures.append(t + ".pdf")
    attach_back(workflow,
                PythonCommand(
                    long_tex,
                    input = {"template": resource_filename("chilin2.modules.phantompeak", "phan.tex"),
                             "figure": figures},
                    output = {"latex": conf.latex_prefix + "_phan.tex"}))

    
def long_tex(input, output, param):
   latex = JinjaTemplateCommand(
       name="mapping quality",
       template=input["template"],
       param={"phan": input["figure"],
              "render_dump": output["latex"]})
   template_dump(latex)

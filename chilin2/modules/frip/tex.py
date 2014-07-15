from samflow.command import ShellCommand, PythonCommand
from samflow.workflow import attach_back
from chilin2.modules.config.helpers import JinjaTemplateCommand, template_dump, json_dump, json_load, underline_to_space, decimal_to_latex_percent
from pkg_resources import resource_filename

import os

def tex_frip(workflow, conf):
    attach_back(workflow,
        PythonCommand(
            load_latex,
            input={"json": conf.json_prefix + "_frip.json",
                   "template": resource_filename("chilin2.modules.frip", "frip.tex"),
                   },
            output={"latex": conf.latex_prefix + "_frip.tex"}))

def load_latex(input, output, param):
    json_dict = json_load(input["json"])
    frip_table = []
    stat = json_dict["stat"]
    for sample in stat: #sample is the keys in the stat dictionary
        frip_table.append([underline_to_space(str(sample)), stat[sample]["info_tag"], stat[sample]["total_tag"], decimal_to_latex_percent(stat[sample]["frip"])])
    latex = JinjaTemplateCommand(
        template=input["template"],
        param={"frip_table": frip_table,
               "render_dump": output["latex"]})
    template_dump(latex)

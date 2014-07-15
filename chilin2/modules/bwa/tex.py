from chilin2.modules.config.helpers import json_load, JinjaTemplateCommand, template_dump, underline_to_space, count_in_million, decimal_to_latex_percent, latex_start, latex_end
import os
from pkg_resources import resource_filename
from samflow.command import ShellCommand, PythonCommand
from samflow.workflow import attach_back


def tex_bwa(workflow, conf):
    attach_back(workflow,
                PythonCommand(
                    long_tex,
                    input = {"template": resource_filename("chilin2.modules.bwa", "bwa.tex"),
                             "figure": conf.prefix + "_bwa_compare.pdf"},
                    output = {"latex": conf.latex_prefix + "_map.tex"}))


def long_tex(input, output, param):
   latex = JinjaTemplateCommand(
       name="mapping quality",
       template=input["template"],
       param={"mappable_ratio_graph": input["figure"],
              "render_dump": output["latex"]})
   template_dump(latex)

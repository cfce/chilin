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
from chilin2.modules.config.helpers import template_dump, JinjaTemplateCommand, decimal_to_latex_percent, json_load

def tex_contamination(workflow, conf):
    all_species = [i for i, _ in conf.items("contamination")]
    tex = attach_back(workflow, PythonCommand(
        latex_contamination,
        input = {"template": resource_filename("chilin2.modules", "contamination/contamination.tex"),
                 "json": conf.json_prefix + "_contam.json"},
        output = {"latex": conf.latex_prefix + "_contam.tex"},
        param = {'id': conf.id, 'layout': 'c'*(len(all_species)+1)}))
    tex.allow_dangling = True
    tex.allow_fail = True

def latex_contamination(input, output, param):
    """
    input bowtie json files generated from mapping to different species
    :param input:
    :param output:
    :param param:
    :return:
    """
    json_dict = json_load(input["json"])

    if len(json_dict["stat"]["meta"]["species"]) < 1:
        section_name = ""
    else:
        section_name = "library_contamination"

    contam_values = json_dict["stat"]["value"]
    for sample in contam_values:
        for species in contam_values[sample]:
            contam_values[sample][species] = decimal_to_latex_percent(contam_values[sample][species])

    library_quality_latex = JinjaTemplateCommand(
        name="library contamination",
        template=input["template"],
        param={"section_name": section_name,
               "library_contamination": json_dict["stat"],
               'prefix_dataset_id': json_dict["param"]['id'],
               'layout': param['layout'],
               "render_dump": output["latex"]})
    template_dump(library_quality_latex)

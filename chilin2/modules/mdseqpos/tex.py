import json
import re
import os
import math
from samflow.command import ShellCommand, PythonCommand
from samflow.workflow import attach_back

from pkg_resources import resource_filename
from chilin2.modules.config.helpers import JinjaTemplateCommand, template_dump, json_load, json_dump

def tex_motif(workflow, conf):
    tex = attach_back(workflow,
                PythonCommand(
                    latex_seqpos,
                    input={"template": resource_filename("chilin2.modules", "mdseqpos/motif.tex"),"json":conf.json_prefix + "_seqpos.json"},
                    output={"tex": conf.latex_prefix + "_motif.tex"},
                    param={"id": conf.id},
                    name = "generating latex of motif info"))
    tex.allow_fail = True
    tex.allow_dangling = True

def latex_seqpos(input, output, param):
    json_dict = json_load(input["json"])
    latex = JinjaTemplateCommand(
        name = "motif finding",
        template = input["template"],
        param = {"motif_table": json_dict["stat"]["satisfied_motifs"][:5], ## use top 5 motif
                 "render_dump": output["tex"]})
    template_dump(latex)

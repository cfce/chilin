from samflow.command import ShellCommand, PythonCommand
from samflow.workflow import attach_back


from chilin2.modules.config.helpers import JinjaTemplateCommand, template_dump, json_load, json_dump

from pkg_resources import resource_filename

def latex_conservation(input, output, param):
    latex = JinjaTemplateCommand(
        name = "conservation",
        template = input["template"],
        param={"conservation_compare_graph": param["prefix"] + "_conserv_compare.pdf",
               "conservation_graph": param["prefix"] + "_conserv.pdf",
               "render_dump": output["latex"]})
    template_dump(latex)

def tex_conserv(workflow, conf):
    tex = attach_back(workflow,
        PythonCommand(
            latex_conservation,
            input={"template": resource_filename("chilin2.modules.conservation", "conservation.tex")},
            output={"latex": conf.latex_prefix + "_conserv.tex"},
            param = {"prefix": conf.prefix}))
    tex.allow_dangling = True
    tex.allow_fail = True

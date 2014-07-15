from samflow.command import ShellCommand, PythonCommand
from samflow.workflow import attach_back
from chilin2.modules.config.helpers import JinjaTemplateCommand, template_dump, json_dump, json_load, underline_to_space
from pkg_resources import resource_filename

import os

def tex_fastqc(workflow, conf):
    attach_back(workflow,
        PythonCommand(
            load_latex,
            input={"json": conf.json_prefix + "_fastqc.json",
                   "template": resource_filename("chilin2.modules.fastqc", "fastqc.tex"),
                   "pdf": conf.prefix + "_raw_sequence_qc.pdf"},
            output={"latex": conf.latex_prefix + "_fastqc.tex"}))
    #these are name, png pairings
    gccontent_graphs = [(nm.replace("_"," "),
                         os.path.join(conf.target_dir, "%s_100k_fastqc" % nm,
                                      "Images","per_sequence_gc_content.png"))\
                            for nm in conf.sample_bases]
    attach_back(workflow,
        PythonCommand(
            load_gc_latex,
            input={"template": resource_filename("chilin2.modules.fastqc", "fastqc_gc.tex"),
                   "gccontent_graphs":gccontent_graphs },
            output={"latex": conf.latex_prefix + "_fastqc_gc.tex"}))

def load_latex(input, output, param):
    json_dict = json_load(input["json"])
    fastqc_summary = []
    stat = json_dict["stat"]
    for sample in stat:
        fastqc_summary.append([underline_to_space(sample), stat[sample]["sequence_length"], stat[sample]["median"]])
    latex = JinjaTemplateCommand(
        template=input["template"],
        param={"section_name": "sequence_quality",
               "path": input["pdf"],
               "fastqc_table": fastqc_summary,
               "fastqc_graph": input["pdf"],
               'prefix_dataset_id': [ underline_to_space(i) for i in stat.keys() ],
               "render_dump": output["latex"]})
    template_dump(latex)

def load_gc_latex(input, output, param):
    latex = JinjaTemplateCommand(
        template=input["template"],
        param={"gccontent_graphs": input["gccontent_graphs"],
               "render_dump": output["latex"]})
    template_dump(latex)

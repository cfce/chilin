from chilin2.modules.config.helpers import json_load, JinjaTemplateCommand, template_dump, underline_to_space, count_in_million, decimal_to_latex_percent, latex_start, latex_end
from samflow.command import ShellCommand, PythonCommand
from samflow.workflow import attach_back
from pkg_resources import resource_filename
from chilin2.modules.summary.qc_json_summary import summary_table



def latex_environ(workflow, conf):
    """
    write out begin and end document
    including packages
    """
    attach_back(workflow,
        PythonCommand(
            latex_start,
            input={"template": resource_filename("chilin2.modules.summary", "begin.tex")},
            output={"latex": conf.latex_prefix + "_begin.tex"},
            param={"id": conf.id,
                   "version": conf.get("basics", "version"),
                   "user": conf.get('basics', 'user'),
                   "bmcard": resource_filename("chilin2.modules.summary", "bmcart.cls").rstrip('.cls')}))

    attach_back(workflow,
        PythonCommand(
            latex_end,
            input={"template": resource_filename("chilin2.modules.summary", "end.tex")},
            output={"latex": conf.latex_prefix + "_end.tex"}))


def merge_latex(workflow, conf):

## begin and end of the docs
    latex_order = [
        "_begin.tex",
        "_summary_table.tex",
    ]
    if conf.long:
        latex_order += ["_fastqc.tex",
                        "_fastqc_gc.tex",
                        "_map.tex",
                        "_conserv.tex",
                        # "_macs2.latex", "_macs2_on_sample.latex",
                        # "_phan.tex",
                        "_motif.tex",
                        "_contam.tex",
                        "_frip.tex",
                        ]
    latex_order.append("_end.tex")

    latex_list = [conf.latex_prefix + i for i in latex_order]
    merge_cmd = attach_back(workflow,
                            ShellCommand(
                                "cat {param[tex]} > {output}",
                                output=conf.prefix + ".tex"))
    merge_cmd.allow_fail = True
    merge_cmd.param = {"tex": " ".join(latex_list)}


def latex_summary_table(input, output, param):
    summary = summary_table(param["conf"])

    latex = JinjaTemplateCommand(
        template=input["template"],
        param={"summary_table": summary,
               "layout": param["layout"],
               "render_dump": output["latex"]})
    template_dump(latex)


def summary_table_latex(workflow, conf):
    summary_tab = attach_back(workflow,
         PythonCommand(
             latex_summary_table,
             input={"template": resource_filename("chilin2.modules.summary", "summary_table.tex")},
             output={"latex": conf.latex_prefix + "_summary_table.tex"},
             param={"conf": conf,
                    "layout": "l"+"c"*(1+len(conf.sample_bases))}))
    summary_tab.allow_fail = True
    summary_tab.allow_dangling = True


def render_pdf(workflow, conf, long = True):
    latex_environ(workflow, conf)
    summary_table_latex(workflow, conf)
    merge_latex(workflow, conf)
    render = attach_back(workflow,
                ShellCommand(
                    # Somehow the pdflatex has to be invoked twice..
                    "{tool} -output-directory {output[dir]} -jobname={param[name]} {input} \
                    && {tool} -output-directory {output[dir]} -jobname={param[name]} {input}",
                    tool="pdflatex",
                    input=conf.prefix + ".tex",
                    # output[pdf] should use "conf.prefix" to have the absolute path
                    output={"dir": conf.target_dir, "pdf": conf.prefix + ".pdf"},
                    # param[name] should use "conf.id" to avoid using absolute path
                    param={"name": conf.id},
                    name="report"))

    render.allow_fail = True

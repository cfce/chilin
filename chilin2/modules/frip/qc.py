from samflow.command import PythonCommand
from samflow.workflow import attach_back

from chilin2.modules.config.helpers import JinjaTemplateCommand, template_dump, r_exec, json_load, underline_to_space, count_in_million, decimal_to_latex_percent, json_dump


def stat_frip(workflow, conf):    # collect frip score
    """
    collect FRiP informative tag number and effective peaks number
    """
    stat = attach_back(workflow, PythonCommand(
        json_frip,
        input={"frip": [t+".frip" for t in conf.sample_targets]},
        output={"json": conf.json_prefix+"_frip.json"},
        param={"samples": conf.sample_bases}))
    stat.allow_fail = True
    stat.allow_dangling = True


def json_frip(input={}, output={}, param={}):    # convert to json
    """
    input is *.frip
    output is conf.json_prefix + "_frip.json"
    param for matching samples
    """
    json_dict = {"stat": {}, "input": input, "output": output, "param": param}
    for i,s in zip(input["frip"], param["samples"]):
        inf = open(i).read().strip().split(",")
        json_dict["stat"][s] = {}
        json_dict["stat"][s]["info_tag"] = int(inf[0])
        json_dict["stat"][s]["total_tag"] = int(inf[1])
        json_dict["stat"][s]["frip"] = float(int(inf[0]))/int(inf[1])
    json_dump(json_dict)

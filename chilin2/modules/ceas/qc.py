import json
import re

from samflow.command import ShellCommand
from samflow.command import PythonCommand
from samflow.workflow import attach_back
from chilin2.modules.config.helpers import JinjaTemplateCommand, template_dump, json_load, json_dump


def stat_ceas(workflow, conf, has_dhs, has_velcro):    # collect meta gene distribution info
    """ Describe peaks' distribution
    ###########################################################################
    DEPRECATED!!!!--see stat_bedAnnotate below
    ###########################################################################
    """
    attach_back(workflow, PythonCommand(
        json_meta,
        input={"meta": conf.prefix + ".meta",
               "top_peaks": 5000},
        output={"json": conf.json_prefix + "_meta.json"},
        param={"id": conf.id},
        name="DHS summary"))

    if has_dhs:
        attach_back(workflow, PythonCommand(
            json_dhs,
            input={"dhs": conf.prefix + ".dhs",
                   "top_peaks": 5000},
            output={"json": conf.json_prefix + "_dhs.json"},
            name="DHS summary"))

    if has_velcro:
        attach_back(workflow, PythonCommand(
            json_velcro,
            input={"velcro": conf.prefix + ".velcro",
                   "top_peaks": 5000},
            output={"json": conf.json_prefix + "_velcro.json"},
            name="Velcro summary"))

def stat_bedAnnotate(workflow, conf, has_dhs, has_velcro):
    """ Describe peaks' distribution
    # collect meta gene distribution info
    """
    attach_back(workflow, PythonCommand(
        json_meta2,
        input={"meta": conf.prefix + ".meta"},
        output={"json": conf.json_prefix + "_meta.json"},
        param={"id": conf.id},
        name="bedAnnotate summary"))

    if has_dhs:
        attach_back(workflow, PythonCommand(
            json_dhs,
            input={"dhs": conf.prefix + ".dhs",
                   "top_peaks": 5000},
            output={"json": conf.json_prefix + "_dhs.json"},
            name="DHS summary"))

    if has_velcro:
        attach_back(workflow, PythonCommand(
            json_velcro,
            input={"velcro": conf.prefix + ".velcro",
                   "top_peaks": 5000},
            output={"json": conf.json_prefix + "_velcro.json"},
            name="Velcro summary"))


def json_meta(input={}, output={}, param={}):
    """
    ###########################################################################
    DEPRECATED!!! see json_meta2
    ###########################################################################
    generate json of promoter, intergenic and exon
    only one output either from merged one or the best one
    overlap percentage info
    """
    f = open(input["meta"])
    content = f.read().strip().split(",")
    exon = content[0]
    intron = content[1]
    inter = content[2]
    promoter = content[3]
    json_dict = {"input": input, "stat": {}, "output": output, "param": param}
    json_dict["stat"]["exon"] = float(exon)
    json_dict["stat"]["intron"] = float(intron)
    json_dict["stat"]["promoter"] = float(promoter)
    json_dict["stat"]["inter"] = float(inter)
    f.close()
    json_dump(json_dict)

def json_meta2(input={}, output={}, param={}):
    """
    generate json of genomic distribution (given by the bedAnnotate output)
    ***THE key difference between json_meta and this fn is that bedAnnotate
    conveniently outputs the distribution as a dictionary of peak counts
    """
    f = open(input["meta"])
    #f = something like: {'Intron': 68017, 'Exon': 7659, 'Intergenic': 73090, 'Promoter': 11229}
    content = eval(f.read())
    total = 0 
    for k in content.keys():
        total += content[k]
    json_dict = {"input": input, "stat": {}, "output": output, "param": param}
    json_dict["stat"]["exon"] = content['Exon']/float(total)
    json_dict["stat"]["intron"] = content['Intron']/float(total)
    json_dict["stat"]["promoter"] = content['Promoter']/float(total)
    json_dict["stat"]["inter"] = content['Intergenic']/float(total)
    f.close()
    json_dump(json_dict)
    

def json_velcro(input={}, output={}, param={}):
    result_dict = {"stat": {}, "input": input, "output": output, "param": param}
    result_dict["stat"] = 1-float(open(input["velcro"]).read().strip())
    json_dump(result_dict)

def json_dhs(input={"top_peaks": "", "dhs_peaks": ""}, output={"json": ""}, param={}):
    result_dict = {"stat": {}, "input": input, "output": output, "param": param}
#    result_dict["stat"] = float(open(input["dhs"]).read().strip())
    content = open(input["dhs"]).read().strip().split(",")
    result_dict["stat"]["overlap"] = int(content[1])
    result_dict["stat"]["number"] = int(content[0])
    json_dump(result_dict)

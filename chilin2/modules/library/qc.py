from samflow.command import PythonCommand
from samflow.workflow import attach_back
from chilin2.modules.config.helpers import json_dump


def stat_pbc(workflow, conf): # collect pbc value
    """
    statistics collected from *.pbc
    """
    attach_back(workflow, PythonCommand(
        json_pbc,
        input = {"pbc": [t + ".pbc" for t in conf.sample_targets]},
        output = {"json": conf.json_prefix + "_pbc.json"},
        param = {"samples":conf.sample_bases}))


def json_pbc(input={}, output={}, param={}): # convert to json format
    """
    input is the target + ".pbc"
    output is the json files conf.json_prefix + "_pbc.json"
    param for matching samples order
    """
    json_dict = {"stat": {}, "input": input, "output": output, "param": param}

    for i, s in zip(input["pbc"], param["samples"]):
        inl = open(i).readlines()[0].strip().split()
        json_dict["stat"][s] = {}
        json_dict["stat"][s]["N1"] = int(inl[0])
        json_dict["stat"][s]["Nd"] = int(inl[1])
        json_dict["stat"][s]["PBC"] = round(float(inl[2]), 3)

    json_dump(json_dict)

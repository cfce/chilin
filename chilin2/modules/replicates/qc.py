from chilin2.modules.config.helpers import JinjaTemplateCommand, template_dump, json_dump, json_load
import json
import os
import re
import sqlite3

from samflow.command import PythonCommand
from samflow.command import ShellCommand
from samflow.workflow import attach_back
from pkg_resources import resource_filename
from chilin2.modules.dbaccessor.caller import db_bwa


def stat_replicates(workflow, conf):  ## replicates peaks and bigwiggle
    """
    input:wigCorrelate of multiple replicates results
          replicates peaks overlap number(percentage: 0.3)
    output: *replicates.json
    """
    stat = attach_back(workflow,
                PythonCommand(
                    json_reps,
                    input = {"cor": conf.prefix+".cor",
                             "overlap": [conf.prefix + "_%s_%s.overlap" % (i, j) for i in range(len(conf.treatment_targets)) for j in range(i+1, len(conf.treatment_targets))]},
                    output = {"json": conf.json_prefix + "_rep.json"},
                    param = {"param": conf.id}))
    stat.allow_fail = True
    stat.allow_dangling = True

def json_reps(input, output, param):
    json_dict = {"stat": {}, "input": input, "output": output, "param": param}
    json_dict['stat']['cor'] = [ float(i.strip().split()[2]) for i in open(input['cor']).readlines() ]
    json_dict["stat"]['overlap'] = [ float(open(i).read().strip()) for i in input['overlap'] ]
    json_dump(json_dict)

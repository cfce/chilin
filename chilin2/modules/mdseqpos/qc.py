import re
import json
import os
import math
from samflow.command import ShellCommand, PythonCommand
from samflow.workflow import attach_back

from pkg_resources import resource_filename
from chilin2.modules.config.helpers import JinjaTemplateCommand, template_dump, json_load, json_dump


def stat_motif(workflow, conf):
    attach_back(workflow,
                PythonCommand(
                    stat_seqpos,
                    input={"seqpos":conf.prefix + "_seqpos/" + "motif_list.json"},
                    output={"json": conf.json_prefix + "_seqpos.json"},
                    param={"prefix": conf.prefix + "_seqpos/seqLogo/", "z_score_cutoff": -1},
                    name = "collect motif info"))


def stat_seqpos(input = {"template": "", "seqpos": ""}, output={"latex_section": ""}, param = {"prefix": "", "z_score_cutoff":-15}):
    """parse mdsepose html file"""
    z_score_cutoff = param["z_score_cutoff"]
    seqpos_html_content = open(input['seqpos']).readlines()
    mdseqpos_result = []

    ## parse motif list json file
    for m in seqpos_html_content:
        mdseqpos_result.append(json.loads(m.strip()))
    satisfied_motif_list = []

    for a_motif in mdseqpos_result:
        if a_motif['seqpos_results']['zscore'] == 'None':
            a_motif['seqpos_results']['zscore'] = 65535
        if a_motif['factors'] == None:
            a_motif['factors'] = ['denovo']
        satisfied_motif_list.append(a_motif)

    satisfied_motif_list.sort(key=lambda x:x['seqpos_results']['zscore'])
    satisfied_count = 0
    top_motifs = []
    for a_motif in satisfied_motif_list:

        if a_motif['id'].find('observed')>0:
            continue
        if satisfied_count == 10:
            break

        # z_score is a negative score, the smaller, the better
        if a_motif['seqpos_results']['zscore'] < z_score_cutoff :
            satisfied_count += 1
            top_motifs.append(a_motif)

    ## choose first 5 motifs to fit into latex document

    for n, _ in enumerate(top_motifs):
        top_motifs[n]["logoImg"] = param["prefix"] + top_motifs[n]['id'] + ".png"

    result_dict = {"stat": {}, "input": input, "output": output, "param": param}
    result_dict["stat"]["satisfied_motifs"] = top_motifs
    json_dump(result_dict)


import json
from samflow.command import ShellCommand, PythonCommand
from samflow.workflow import attach_back

import re
import random
import os
import json
from chilin2.modules.config.helpers import JinjaTemplateCommand, template_dump, json_load, underline_to_space, decimal_to_latex_percent, json_dump

def stat_contamination(workflow, conf):
    all_species = [i for i, _ in conf.items("contamination")]
    summ = []
    for target in conf.sample_targets:
        summ.append([ (target + species + "_mapped." + conf.mapper, target + species + "_total." + conf.mapper) for species in all_species ])

    attach_back(workflow,
                PythonCommand(json_contamination,
                              input={"summaries": summ},
                              output={"json": conf.json_prefix + "_contam.json"},
                              param={"samples": conf.sample_bases,
                                     "id": conf.id,
                                     "species": all_species},
                              name = "stat contamination"))

## summary of library contamination
def json_contamination(input = {"summaries": [[]]}, output = {"json": ""}, param = {"samples": "", "species": "", "id": ""}):
    library_contamination = {}
    library_contamination["meta"] = {"sample": param["id"], "species": param["species"]}
    library_contamination["value"] = {}
    for a_summary, s in zip(input["summaries"], map(underline_to_space, param["samples"])):
        ## each bowtie_summary has several species information
        library_contamination["value"][s] = {}
        for i, j in zip(a_summary, param["species"]):
            ## species 1, species2, species3
            mapped = int(open(i[0]).readlines()[2].strip().split()[0])
            total = int(open(i[1]).read().strip())
            library_contamination["value"][s][j] = float(mapped)/total

    json_dict = {"stat": {}, "input": input, "output": output, "param": param}
    json_dict["stat"] = library_contamination
    json_dump(json_dict)

def bowtie_summary_parse(input=""): ## extract bowtie standard error
    """
    This is used for library contamination only
    """
    summary_content = open(input).read()
    print(summary_content)
    print("_" * 100)
    total_reads = int(re.findall(r"# reads processed: (.*)\n", summary_content)[0])

    # WARN: this is `unique mappable reads` as we set `-m 1`
    mappable_reads = int(re.findall(r"# reads with at least one reported alignment: (.*) \(", summary_content)[0])

    mappable_rate = float(mappable_reads) / float(total_reads)

    return {"total_reads": total_reads, "mappable_reads": mappable_reads, "mappable_rate": mappable_rate}


## Mapper changes, so background bowtie mapping statitics are obsolete
#def stat_bowtie(input={"bowtie_summaries": [], "dbaccessor":"", "template": ""},
#                output={"json": "", "R": "", "pdf": ""},
#                param={"sams": [], }):
#    """ sams = [{'name1':a, 'total1': 5...}, {'name2':c, 'total2': 3...}...] **args """
#
#    # unique location is in text_macs2_summary part
#    json_dict = {"stat": {}, "input": input, "output": output, "param": param}
#
#    db = sqlite3.connect(input["dbaccessor"]).cursor()
#    db.execute("select map_ratio from mapping")
#    historyData = [str(i[0]) for i in (db.fetchall())]
#    bowtie_summaries = {"total_reads": [],
#                        "mappable_reads": [],
#                        "mappable_rate": []}
#
#    for summary, sam in zip(input["bowtie_summaries"], param["sams"]):
#        json_dict["stat"][sam] = _bowtie_summary_parse(summary)
#        json_dict["stat"][sam]["cutoff"] = 5000000 # mappable reads
#        json_dict["stat"][sam]["judge"] = "Pass" if json_dict["stat"][sam]["mappable_reads"] >= 5000000 else "Fail"
#
#    mappable_rates = [json_dict["stat"][i]["mappable_rate"] for i in json_dict["stat"]]
#
#    mappable_rate_R = JinjaTemplateCommand(template=input["template"],
#        param={'historic_data': historyData,
#               'current_data': mappable_rates,
#               'ids': [ os.path.basename(i) for i in param["sams"]],
#               'cutoff': 0.5,
#               'section': 'Unique mapped rate',
#               "need_smooth_curve": True,
#               "render_dump": output["R"],
#               "pdf": output["pdf"], })
#    template_dump(mappable_rate_R)
#    r_exec(mappable_rate_R)
#
#    with open(output["json"], "w") as f:
#        json.dump(json_dict, f, indent=4)

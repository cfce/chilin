import re
from samflow.command import ShellCommand, PythonCommand
from samflow.workflow import attach_back
from chilin2.modules.config.helpers import JinjaTemplateCommand, template_dump, json_dump, json_load, underline_to_space, r_exec
from chilin2.modules.dbaccessor.caller import db_fastqc
from pkg_resources import resource_filename


def stat_fastqc(workflow, conf):  # collect raw reads quality and GC contents
    """
    long: generate long pages or not
    """
    attach_back(workflow,
                PythonCommand(
                    json_fastqc,
                    input={"fastqc_summaries": [target + "_100k_fastqc/fastqc_data.txt" for target in conf.sample_targets]},
                    output={"json": conf.json_prefix + "_fastqc.json"},
                    param={"ids": conf.sample_bases,
                           "id": conf.id},
                    name = "collect fastqc results"))
    
    if conf.long:  ## prepare long document images and tex
        attach_back(workflow,
        PythonCommand(fastqc_detailed_figure,
                      input = {"dbaccessor": resource_filename("chilin2.modules.dbaccessor", "ChiLinQC.db"),
                               "template": resource_filename("chilin2.modules.summary", "R_culmulative_plot.R"), 
                               "json": conf.json_prefix + "_fastqc.json"},
                      output = {"R": conf.prefix + "_raw_sequence_qc.R",
                                "pdf": conf.prefix + "_raw_sequence_qc.pdf"},
                      param={"ids": conf.sample_bases}))
        

def json_fastqc(input={"fastqc_summaries": []},
                output={"R": "", "json": "", "pdf": ""},
                param={"ids": [], "id": ""}):
    json_dict = {"stat": {}, "input": input, "output": output, "param": param}
    stat = json_dict["stat"]
    for a_summary, a_id in zip(input["fastqc_summaries"], param["ids"]):
        parsed = _fastqc_parse(input=a_summary)
        stat[a_id] = {}
        stat[a_id]["median"] = parsed["median"]
        stat[a_id]["sequence_length"] = parsed["sequence_length"]
    json_dump(json_dict)

def fastqc_detailed_figure(input = {"template": "", "dbaccessor": "", "json": ""},output = {"R":"", "pdf":""},
                           param={"ids": ""}):
    quality_medians = []
    stat = json_load(input["json"])['stat']
    for s in param['ids']:
        quality_medians.append(stat[s]['median'])

    ## from dbaccessor module
    historic_data = db_fastqc(input["dbaccessor"])
    # The table of fastqc_summary that will be used for rendering
    # Col 1: sample ID
    # Col 2: sequence length
    # Col 3: median of sequence quality
    fastqc_dist_r = JinjaTemplateCommand(
        template=input["template"],
        param={'historic_data': historic_data,
               'current_data': quality_medians,
               'ids': [underline_to_space(i) for i in param["ids"]],
               'cutoff': 25,
               "need_smooth_curve": True,
               "section": "FastQC score distribution",
               "pdf": output["pdf"],
               "render_dump": output["R"]})
    template_dump(fastqc_dist_r)
    r_exec(fastqc_dist_r)

def _fastqc_parse(input, output=None, param=None):
    data = open(input).readlines()
    sequence_length = 0
    quality_dict = {}
    in_seq_quality_section = False
    for line in data:
        if re.search(r"^Sequence length", line):
            assert sequence_length == 0
            sequence_length = int(re.findall(r"^Sequence length\t(\d+)", line)[0])
        elif re.search(r"^>>Per sequence quality", line):
            assert not in_seq_quality_section
            in_seq_quality_section = True
            continue
        if re.search(r"^>>END_MODULE", line) and in_seq_quality_section:
            in_seq_quality_section = False

        if (not line.startswith("#")) and in_seq_quality_section:
            sequence_quality = re.findall("^(\w+)\t(\w+)", line)[0]
            quality_dict[sequence_quality[0]] = float(sequence_quality[1])
    total = sum(quality_dict.values())
    n = 0
    for item in sorted(quality_dict.items(), key=lambda e: e[0], reverse=True):
        n = n + item[1]
        if n / total > 0.5:
            median = int(item[0])
            break
    return {"sequence_length": sequence_length,
            "median": median}



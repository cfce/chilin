import os
import sys
from chilin2.modules.config.helpers import JinjaTemplateCommand, template_dump, json_dump
from samflow.command import PythonCommand
from samflow.workflow import attach_back

def stat_macs2(workflow, conf):   # collect peaks
    """
    merged peaks and replicates peaks
    high confident peaks
    duplicates level
    """
    xls = conf.prefix + "_peaks.xls" if conf.get("macs2", "type") in ["both", "narrow"] else conf.prefix + "_b_peaks.xls"

    stat = attach_back(workflow, PythonCommand(
        json_macs2,
        input = {"macs2_peaks_xls": xls},
        output = {"json":  conf.json_prefix + "_macs2.json"},
        param = {"id": conf.id}))
    stat.allow_fail = True
    stat.allow_dangling = True


def stat_macs2_on_rep(workflow, conf):
    if conf.get("macs2", "type") in ["both", "narrow"]:
        xls = [ t + "_peaks.xls" for t in conf.treatment_targets ]
    else:
        xls = [ t + "_b_peaks.xls" for t in conf.treatment_targets ]
    if len(conf.treatment_targets) > 1:
        stat = attach_back(workflow, PythonCommand(
            json_macs2_on_reps,
            input = {"all_peak_xls": xls},
            output = {"json":  conf.json_prefix + "_macs2_rep.json"},
            param = {"samples": conf.treatment_bases}))
        stat.allow_fail = True
        stat.allow_dangling = True

def json_macs2(input={"macs2_peaks_xls": ""}, output={"json": ""}, param={"id": ""}):
    """
    input macs2 _peaks.xls
    output conf.json_prefix + "_macs2.json"
    """
    json_dict = {"stat": {}, "input": input, "output": output, "param": param}
    if os.path.exists(input['macs2_peaks_xls']): ## in case only broad peaks would break down sometimes, narrowPeak very seldom no peaks
        json_dict["stat"] = _peaks_parse(input["macs2_peaks_xls"])
        json_dump(json_dict)

def json_macs2_on_reps(input={"all_peak_xls": []}, output={"json": ""}, param={"samples": []}):
    """
    collect replicates macs2 info to json files
    compared to merged one, collect redundant ratio with --keep-dup 1 option
    """
    json_dict = {"stat": {}, "input": input, "output": output, "param": param}
    parsed = []
    for i in input["all_peak_xls"]:
        if os.path.exists(i): ## in case only broad peaks would break down sometimes, narrowPeak very seldom no peaks
            parsed.append(_peaks_parse(i))

    if all(map(os.path.exists, input['all_peak_xls'])):
        for sample, stat in zip(param["samples"], parsed):
            json_dict["stat"][sample] = stat
        json_dump(json_dict)

def _peaks_parse(input):
    total = 0
    fc20n = 0
    fc10n = 0
    peaks_info = {}
    with open(input) as peaks_xls:
        for line in peaks_xls:
            if line.startswith('# tags after filtering in treatment'):
                # tags after filtering in treatment: 13438948
                peaks_info["treat_unic"] = int(line.strip().split(":")[1])
            if line.startswith("# Redundant rate in treatment: "):
                peaks_info["treat_unic_ratio"] = 1 - float(line.split(':')[1])

            if line.startswith('# tags after filtering in control'):
                peaks_info["control_unic"] = int(line.strip().split(":")[1])
            if line.startswith("# Redundant rate in control: "):
                peaks_info["control_unic_ratio"] = 1 - float(line.split(':')[1])

            if line.startswith('# d'):
                peaks_info["distance"] = int(line.strip().split("=")[1])
            if line.strip() != "" and not line.startswith("#") and not line.startswith("chr\t"):
                l = line.strip().split("\t")
                total += 1
                ## column 7th denotes fold change value
                fc = float(l[7])
                if fc >= 20:
                    fc20n += 1
                if fc >= 10:
                    fc10n += 1
            if line.startswith("# qvalue cutoff"):
                q_value_cutoff = float(line.split('=')[1])
            if line.startswith("# d"): # parse shift-size, # d =
                shift_size = int(line.strip().split("=")[1])/2

    peaks_info["totalpeak"] = total
    peaks_info["peaksge20"] = fc20n
    peaks_info["peaksge10"] = fc10n
    if peaks_info["totalpeak"] >= 200:
        peaks_info["peaksge20ratio"] = peaks_info["peaksge20"] / peaks_info["totalpeak"]
        peaks_info["peaksge10ratio"] = peaks_info["peaksge10"] / peaks_info["totalpeak"]
    elif 0 < peaks_info["totalpeak"] < 200:
        peaks_info["peaksge20ratio"] = peaks_info["peaksge20"] / peaks_info["totalpeak"]
        peaks_info["peaksge10ratio"] = peaks_info["peaksge10"] / peaks_info["totalpeak"]
        print >> sys.stderr, "Warning: peaks so few for motif scan"
    elif peaks_info["totalpeak"] == 0 :
        peaks_info["peaksge20ratio"] = 0
        peaks_info["peaksge10ratio"] = 0
        print >> sys.stderr, "Warning: 0 peaks"

    peaks_info["qvalue"] = q_value_cutoff
    peaks_info["shiftsize"] = shift_size
    return peaks_info


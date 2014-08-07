#!/usr/bin/env python

"""
Time: 
Thu Mar  6 14:30:19 EST 2014

Description:
    get enrichment ratio on union DHS, exons and promoters

python version too slow
bigWigAverageBed only report average values for each regions
Use CoverageBed with input BAM files to get raw counts in meta regions and ratio over total reads
"""
from samflow.command import ShellCommand
from samflow.workflow import attach_back
from chilin2.modules.config.helpers import json_dump
from samflow.command import PythonCommand


def read_enrichment_on_meta(workflow, conf):
    """ total reads enrichment in exon, promoter and union DHS regions
    """
    try:
        has_dhs = conf.get(conf.get("basics", "species"), "dhs")
    except:
        has_dhs = ""
    import os
    for t in conf.sample_targets:
        attach_back(workflow, ShellCommand(
            ## use bash modules/ceas/meta_info.sh to get latest annotations
            ## Qian: coverage should not consider the strand information for promoters and exons, 5 column for reads count
            ## because of no strand information, we plus the 4th column
            """
            {tool} -abam {input[bam]} -b {param[exon]} -counts | awk \'{{n+=$4}}END{{print n}}\' - > {output[exon]}
            {tool} -abam {input[bam]} -b {param[promoter]} -counts | awk \'{{n+=$4}}END{{print n}}\' - > {output[promoter]}
            """,
            tool = "coverageBed",
            input = {"bam": t + "_4000000.bam" if conf.down else t + ".bam"},
            output = {"exon":t+".enrich.exon",
                      "promoter": t+".enrich.promoter"},
            param = {"promoter": os.path.join(conf.target_dir, "gene.bed_promoter"), "exon": os.path.join(conf.target_dir, "gene.bed_exon")}))

        if has_dhs:
            ## Qian: dhs has no strand problem, use 4th column as reads count
            attach_back(workflow, ShellCommand(
            ## 4 column reads count
            """
            {tool}  -abam {input[bam]} -b {param[dhs]} -counts | awk \'{{n+=$4}}END{{print n}}\' - > {output[dhs]}
            """,
            tool = "coverageBed",
            input = {"bam": t + "_4000000.bam" if conf.down else t + ".bam",
                    "dhs":conf.get_path(conf.get("basics", "species"), "dhs")},
            output = {"dhs": t+".enrich.dhs"}, param = {"dhs": conf.get_path(conf.get("basics", "species"), "dhs")},
            ))

    attach_back(workflow, PythonCommand(
        enrich_in_meta,
        input = {"exon":[ t+".enrich.exon" for t in conf.sample_targets ],
                 "promoter": [ t+".enrich.promoter" for t in conf.sample_targets ],
                 "mapped": [ t+"_mapped.bwa" for t in conf.sample_targets ]},
        output = {"json": conf.json_prefix + "_enrich_meta.json"},
        param = {"samples": conf.sample_bases, "id":conf.id, "has_dhs":has_dhs,
                 "dhs": [ t+".enrich.dhs" for t in conf.sample_targets ]}))

def enrich_in_meta(input = {'exon':'','dhs':'','promoter':'', "mapped": ""}, output = {"json": ""}, param = {'id':"", 'samples':""}):
    """ enrichment in meta regions
    """
    json_dict = {"stat": {}, "input": input, "output": output, "param":param}


    for n, s in enumerate(param['samples']):

        mapped = float(open(input["mapped"][n]).readlines()[2].split()[0])
        json_dict['stat'][s] = {}
        json_dict['stat'][s]['exon'] = float(open(input['exon'][n]).read().strip())/mapped
        json_dict['stat'][s]['promoter'] = float(open(input['promoter'][n]).read().strip())/mapped
        if param['has_dhs']:
            json_dict['stat'][s]['dhs'] = float(open(param['dhs'][n]).read().strip())/mapped
        else:
            json_dict['stat'][s]['dhs'] = 0
    json_dump(json_dict)


#!/usr/bin/env python

"""
Time:
Thu Mar  6 14:30:19 EST 2014

Description:
    get enrichment ratio on union DHS, exons and promoters

import from frip module using intersectBed
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
        enrich = attach_back(workflow, ShellCommand(
            """
            exon=$(bedtools intersect -f {param[p]} -wa -u -abam {input[bam]} -b {param[exon]} -bed | wc -l)
            promoter=$(bedtools intersect -f {param[p]} -wa -u -abam {input[bam]} -b {param[promoter]} -bed | wc -l)
            total=$(samtools flagstat {input[bam]} | head -1 | cut -d" " -f1)
            echo $exon,$promoter,$total > {output[meta]}
            """,
            tool = "coverageBed",
            input = {"bam": t + "_4000000.bam" if conf.down else t + ".bam"},
            output = {"meta":t+".enrich.meta"},
            param = {"promoter": os.path.join(conf.target_dir, "gene.bed_promoter"), 
                     "p": "1E-9",
                     "exon": os.path.join(conf.target_dir, "gene.bed_exon")}))
        enrich.allow_dangling = True
        enrich.allow_fail = True

        if has_dhs:
            dhs = attach_back(workflow, ShellCommand(
            """
            dhs=$(bedtools intersect -f {param[p]} -wa -u -abam {input[bam]} -b {param[dhs]} -bed | wc -l)
            total=$(samtools flagstat {input[bam]} | head -1 | cut -d" " -f1)
            echo $dhs,$total > {output[dhs]}
            """,
            tool = "coverageBed",
            input = {"bam": t + "_4000000.bam" if conf.down else t + ".bam",
                    "dhs":conf.get_path(conf.get("basics", "species"), "dhs")},
            output = {"dhs": t+".enrich.dhs"}, param = {"p": "1E-9","dhs": conf.get_path(conf.get("basics", "species"), "dhs")},
            ))
            dhs.allow_fail = True
            dhs.allow_dangling = True

    em = attach_back(workflow, PythonCommand(
        enrich_in_meta,
        input = {"meta":[ t+".enrich.meta" for t in conf.sample_targets ],
                 "mapped": [ t+"_mapped.bwa" for t in conf.sample_targets ]}, ## use 4M reads for down sampling ones, and all reads instead
        output = {"json": conf.json_prefix + "_enrich_meta.json"},
        param = {"samples": conf.sample_bases, "id":conf.id, "has_dhs":has_dhs,
                 "down": conf.down,
                 "dhs": [ t+".enrich.dhs" for t in conf.sample_targets ]}))
    em.allow_fail = True
    em.allow_dangling = True

def enrich_in_meta(input = {'meta':'', 'mapped':''}, output = {"json": ""}, param = {'dhs': '', 'down': '', 'has_dhs':'', 'id':"", 'samples':""}):
    """ enrichment in meta regions
    """
    json_dict = {"stat": {}, "input": input, "output": output, "param":param}
    for n, s in enumerate(param['samples']):
        ## total mapped reads

        mapped = float(open(input["mapped"][n]).readlines()[2].split()[0])
        json_dict['stat'][s] = {}

        meta = open(input['meta'][n]).read().strip().split(",")
        meta = map(float, meta)
        if not param["down"]:
            json_dict['stat'][s]['exon'] = meta[0]/mapped
            json_dict['stat'][s]['promoter'] = meta[1]/mapped ## use all mapped reads
        else:
            json_dict['stat'][s]['exon'] = meta[0]/meta[2]
            json_dict['stat'][s]['promoter'] = meta[1]/meta[2] ## use 4M reads

        if param['has_dhs']:
            dhs = open(param["dhs"][n]).read().strip().split(",")
            dhs = map(float, dhs)
            if not param["down"]:
                json_dict['stat'][s]['dhs'] = dhs[0]/mapped
            else:
                json_dict['stat'][s]['dhs'] = dhs[0]/dhs[1]

    json_dump(json_dict)


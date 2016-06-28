from chilin2.modules.config.helpers import json_load, JinjaTemplateCommand, template_dump, underline_to_space, count_in_million, decimal_to_latex_percent, latex_start, latex_end
import os

from samflow.workflow import attach_back
from pkg_resources import resource_filename


def _stat(json_path):
    return json_load(json_path)["stat"]

def summary_table(conf):
    table = []
    exist = os.path.exists
    pre = conf.json_prefix
    samples = conf.sample_bases
    print(conf.control_targets)

    ## summary table header
    if len(conf.control_targets) > 0:
        table.append(["Metrics ",
                      "&".join(['treat' + str(i) for i in range(1, len(conf.treatment_targets)+1)]),
                      "&".join(['control' + str(i) for i in range(1, len(conf.control_targets)+1)])])
    else:
        table.append(["Metrics ",
                      "&".join(['treat' + str(i) for i in range(1, len(conf.treatment_targets)+1)])])

    js = pre + "_fastqc.json"
    ## not skip fastqc and do exist json file
    if exist(js):
        stat = _stat(js)
        fastqc = ["FastQC"]
        for s in samples:
            fastqc_s = stat[s]['median']
            if  fastqc_s >= 25:
                fastqc.append("\color{blue} %s" % fastqc_s)
            else:
                fastqc.append("\color{red} %s" % fastqc_s)
        table.append(fastqc)

    ## not skip fastqc and do exist json file
    mapped = None
    js = pre + "_map.json"
    if exist(js):
        stat = _stat(js)
        if conf.pe:
            labels = '(PE)'
        else:
            labels = ''

        table.append(["Original total reads %s" % labels] + [ count_in_million(stat[s]["total"]) for s in samples ])

        mapping = ["Unique mapped reads %s" % labels]
        for s in samples:
            ratio = float(stat[s]["mapped"])/stat[s]["total"]
            if stat[s]["mapped"] >= 4e6 and ratio >= 0.6:
                color = 'blue'
            else:
                color = 'red'

            mapping.append("\color{%s} %s (%s)" % (color, count_in_million(stat[s]["mapped"]), decimal_to_latex_percent(ratio)))

        table.append(mapping)
        mapped = stat[s]['mapped']

    js = pre + "_pbc.json"
    samples = conf.sample_bases
    def NRF(Nd, map_reads): ## encode cutoff 0.8
        if conf.down:
            if float(map_reads) > 4e6:
                return float(Nd)/4e6
            else:
                return float(Nd)/map_reads
        else:
            return float(Nd)/map_reads

    if exist(js):
        stat = _stat(js)
        if mapped:
            if conf.down:
                ## distinct location => unique locations(with reads >= 1)
                nrf = ["Unique locations of 4M reads"] 
                n1 = ["Locations with only 1 read from 4M reads, number (ratio)"]
                for s in samples:
                    nrf_s = NRF(stat[s]['Nd'],float(mapped))
                    if nrf_s >= 0.7:
                        labels = 'blue'
                    else:
                        labels = 'red'
                    nrf.append("\color{%s} %s (%s)" % (labels, count_in_million(stat[s]["Nd"]), decimal_to_latex_percent(nrf_s)))
                    ## N1 locations
                    n1.append("%s (%s)" % (count_in_million(stat[s]["N1"]), decimal_to_latex_percent(NRF(stat[s]['N1'],float(mapped)))))
                table.append(nrf)
                table.append(n1)
            else:
                ## distinct location => unique locations(with reads >= 1)
                table.append(["Unique locations of total reads"] + [ "%s (%s)" % (count_in_million(stat[s]["Nd"]), decimal_to_latex_percent(NRF(stat[s]['Nd'],float(mapped)))) for s in samples ])
                ## N1 locations
                table.append(["Locations with only 1 read from total reads, number (ratio)"] + [ "%s (%s)" % (count_in_million(stat[s]["N1"]), decimal_to_latex_percent(NRF(stat[s]['N1'],float(mapped)))) for s in samples ])
        else:
            ## output N1 locations(locations with only one read)
            table.append(
                ["Unique locations of 4M reads"] + [ "%s" % count_in_million(stat[s]["N1"]) for s in samples ])

        ## PBC1 = N1/Nd
        if conf.down:
            pbcs = ["PBC of 4M reads"]
            for s in samples:
                pbc = stat[s]["PBC"]
                if pbc >= 0.8:
                    color = 'blue'
                else:
                    color = 'red'
                pbcs.append("\color{%s} %s" % (color, decimal_to_latex_percent(pbc)))
            table.append(pbcs)
        else:
            table.append(
                ["PBC of total reads"] + [ decimal_to_latex_percent(stat[s]["PBC"]) for s in samples ])

    # js = pre + "_phan.json"
    # if exist(js):
    #     stat = _stat(js)
    #     if conf.down:
    #         table.append(
    #             ["NSC of 4M reads"] + [ str(round(float(stat[s]["NSC"]), 3)) for s in samples ])
    #         table.append(
    #             ["RSC of 4M reads"] + [ str(round(float(stat[s]["RSC"]), 3)) for s in samples ])
    #         table.append(
    #             ["Qtag of 4M reads"] + [ str(stat[s]["Qtag"]) for s in samples ])
    #         table.append(
    #             ["Fragment size of 4M reads"] + [ str(stat[s]["frag"]) for s in samples ])
    #     else:
    #         table.append(
    #             ["NSC of Total reads"] + [ str(round(float(stat[s]["NSC"]), 3)) for s in samples ])
    #         table.append(
    #             ["RSC of Total reads"] + [ str(round(float(stat[s]["RSC"]), 3)) for s in samples ])
    #         table.append(
    #             ["Qtag of Total reads"] + [ str(stat[s]["Qtag"]) for s in samples ])
    #         ## fragment size based on all reads
    #         table.append(
    #             ["Fragment size of Total reads"] + [ str(stat[s]["frag"]) for s in samples ])

    js = pre + "_frag.json"
    samples = conf.treatment_bases
    if exist(js):
        stat = _stat(js)
        if conf.down:
            table.append(["Fragment size of 4M reads"] + [ str(stat[s]) for s in samples ])
        else:
            table.append(["Fragment size of total reads"] + [ str(stat[s]) for s in samples ])

    js = pre + "_enrich_meta.json"
    samples = conf.sample_bases
    if exist(js):
        stat = _stat(js)
        if conf.down:
            table.append(["DHS/Promoter/Exon ratio of 4M reads"] + [ str(decimal_to_latex_percent(stat[s]['dhs']))+"/"+str(decimal_to_latex_percent(stat[s]['promoter']))+"/"+str(decimal_to_latex_percent(stat[s]['exon']))  for s in samples ])
        else:
            table.append(["DHS/Promoter/Exon ratio of total reads"] + [ str(decimal_to_latex_percent(stat[s]['dhs']))+"/"+str(decimal_to_latex_percent(stat[s]['promoter']))+"/"+str(decimal_to_latex_percent(stat[s]['exon']))  for s in samples ])

    js = pre + "_frip.json"
    samples = conf.sample_bases
    if exist(js):
        stat = _stat(js)
        if conf.down:
            if conf.frip:
                label = '5M'
            else:
                label = '4M'
            frips = ["FRiP of %s non-chrM reads" % label]
            for s in samples:
                frip = stat[s]['frip']
                if frip >= 0.01:
                    color = 'blue'
                else:
                    color = 'red'
                frips.append("\color{%s} %s" % (color, str(decimal_to_latex_percent(frip))))
        else:
            frips = ["FRiP of total non-chrM reads"] + [ str(decimal_to_latex_percent(stat[s]['frip'])) for s in samples ]
        table.append(frips)

    ## IP layer
    samples = conf.treatment_bases

    js = pre + "_macs2_rep.json"
    if exist(js):
        stat = _stat(js)
        table.append(
            ["Replicates total peaks"] +
            [stat[s]["totalpeak"] for s in samples])
        table.append(
            ["Replicates 10 fold  confident peaks"] +
            [stat[s]["peaksge10"] for s in samples])
        table.append(
            ["Replicates 20 fold  confident peaks"] +
            [stat[s]["peaksge20"] for s in samples])

    js = pre + "_rep.json"
    if exist(js):
        stat = _stat(js)

        if stat["cor"][0] >= 0.6:
            color = 'blue'
        else:
            color = 'red'
        table.append(["Replicates reads correlation/replicates peaks overlap"] + \
                     ['\multicolumn{%s}{c}{%s}' % (str(len(conf.sample_bases)), '\color{%s} ' % (color) + "Correlation " + '\t'.join(map(lambda x: str(round(x,2)), stat["cor"])) + '/Overlap ' + '\t'.join(map(lambda x:str(int(x)), stat["overlap"]))) ])

    ## Pool layer
    js = pre + "_macs2.json"
    if exist(js):
        stat = _stat(js)
        table.append(
            ["Merged total/10 fold/20 fold peaks"] + \
            ['\multicolumn{%s}{c}{%s}' % (str(len(conf.sample_bases)), "/".join(map(str,[stat["totalpeak"], stat["peaksge10"], stat["peaksge20"]])))])

    js = pre + "_velcro.json"
    if exist(js):
        stat = _stat(js)
        table.append(
            ["Top peaks not overlap with blacklist regions ratio"] + ['\multicolumn{%s}{c}{\parbox[c][][c]{2.2cm}{%s}}' % (str(len(conf.sample_bases)), decimal_to_latex_percent(stat))])

    js = pre + "_dhs.json"
    if exist(js):
        stat = _stat(js)
        ratio = float(stat["overlap"])/stat["number"]
        if ratio >= 0.7:
            color = 'blue'
        else:
            color = 'red'

        table.append(
            ["Top peaks overlap with union DHS number (ratio)"] + ['\multicolumn{%s}{c}{\parbox[c][][c]{2.2cm}{\color{%s} %s (%s)}}' % (str(len(conf.sample_bases)), color, stat["overlap"], decimal_to_latex_percent(float(stat["overlap"])/stat["number"]))])

    js = pre + "_meta.json"
    if exist(js):
        stat = _stat(js)
        table.append(["Exon/Intron/Intergenic/Promoter ratio of peak summits"] + ['\multicolumn{%s}{c}{\parbox[c][][c]{2.2cm}{%s}}' % (str(len(conf.sample_bases)), "/".join(map(decimal_to_latex_percent, [stat["exon"], stat["intron"], stat["inter"], stat["promoter"]])))])

    js = pre + "_conserv.json"
    if exist(js):
        stat = _stat(js)
        table.append(["Top peaks conservation plot"] + ['\multicolumn{%s}{c}{\parbox[c][][c]{2in}{\includegraphics[scale=0.25]{%s}}}\\' % (str(len(conf.sample_bases)), conf.prefix + "_conserv.pdf")])

    js = pre + "_seqpos.json"
    if exist(js):
        stat = _stat(js)['satisfied_motifs']
        ## multicolumn of latex motif row
        try:
            table.append(["Top peaks motif analysis"] + ['\multicolumn{%s}{c}{\parbox[c][][c]{2in}{%s}}' % (str(len(conf.sample_bases)), ' , '.join(stat[0]['factors']) + " "+ str(round(stat[0]['seqpos_results']['zscore'], 1)))])
        except IndexError:
            table.append(["Top peaks motif analysis"] + ['\multicolumn{%s}{c}{\parbox[c][][c]{2.2cm}{%s}}' % (str(len(conf.sample_bases)), "No significant motif found")])
            pass

    return table

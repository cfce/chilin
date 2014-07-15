from samflow.command import ShellCommand, PythonCommand
from samflow.workflow import attach_back

from pkg_resources import resource_filename
from chilin2.modules.config.helpers import JinjaTemplateCommand, template_dump, json_load, json_dump
from chilin2.modules.mdseqpos.qc import stat_motif
from chilin2.modules.mdseqpos.tex import tex_motif


def seqpos(workflow, conf):
    """

    :param workflow:
    :param conf:
    :return:
    """
    get_top_peaks = attach_back(workflow,
                                ShellCommand(
                                    "{tool} -n {param[peaks]} {input} | cut -f 1,2,3,4,9 > {output}",
                                    tool="head",
                                    input=conf.prefix + "_sort_summits.bed" if conf.get("macs2", "type") in ["both", "narrow"] else conf.prefix + "_b_sort_peaks.broadPeak",
                                    output=conf.prefix + "_peaks_top_motif.bed",
                                    param={"peaks": 5000},
                                    name="top summits for mdseqpos"))
    get_top_peaks.update(param=conf.items("seqpos"))

    if conf.get("basics", "species"):
        species = conf.get("basics", "species")
    else:
        species = 'hg19'

    if conf.get("tool", "mdseqpos"):
        mdseqpos_bin = conf.get("tool", "mdseqpos")
    else:
        mdseqpos_bin = "MDSeqPos.py"

    mdseqpos = attach_back(workflow,
                           ShellCommand(
            "{tool} -d  -w 600  -p 0.001  -m cistrome.xml -O {output[result_dir]} {input} {param[species]} ",
            tool=mdseqpos_bin,
            input=conf.prefix + "_peaks_top_motif.bed",
            output={"result_dir": conf.prefix + "_seqpos", "seqpos": conf.prefix + "_seqpos/" + "motif_list.json"},
            param={"species": species}, name="motif finding"))

    # ERR: note these values are hard-coded!
    # need to make it more universal for all species
    mdseqpos.update(param=conf.items("seqpos"))

    ## QC part
    stat_motif(workflow, conf)
    if conf.long:
        tex_motif(workflow, conf)

## regulatory potential

from samflow.command import ShellCommand
from samflow.workflow import attach_back
from pkg_resources import resource_filename

def reg_potential(workflow, conf):
    """
    """
    if conf.get('macs2', 'type') in ["both", "narrow"]:
        get_fold5_peaks = attach_back(workflow, 
                                      ShellCommand(
                                        """
                                        {tool} \'($1 != "chr" && $1 !="#" && $8>=5)\' {input} | awk \'($1,$2,$3,$10,$9)\' -f 1,2,3,10,9 > {output}
                                        """,
                                        tool = 'awk',
                                        input = conf.prefix + '_peaks.xls',
                                        output = conf.prefix + '_5foldpeaks.bed'
                                     ))
        reg = attach_back(workflow,
                          ShellCommand(
                              "{tool} -t {input[peaks]} -g {param[geneTable]} -n {param[prefix]} -d {param[dist]}",
                              tool = "RegPotential.py",
                              input = {"peaks": conf.prefix + '_5foldpeaks.bed'},
                              output = {"potential": conf.prefix + "_gene_score_5fold.txt"},
                              param = {"geneTable": conf.get_path(conf.get("basics", "species"), "geneTable"),
                                       "tool": resource_filename("chilin2.modules", "regulatory/RegPotential.py"),
                                       "prefix": conf.prefix + "_gene_score_5fold.txt",
                                       "dist": 100000},
                              name = "Regulatory Potential"))
        reg.update(param=conf.items("reg_potential"))
        reg.allow_fail = True
        reg.allow_dangling = True

    get_top_peaks = attach_back(workflow,
                                ShellCommand(
                                    "{tool} -n {param[peaks]} {input} | cut -f 1,2,3,4,9> {output}",
                                    tool="head",
                                    input=conf.prefix + "_sort_peaks.narrowPeak" if conf.get("macs2", "type") in ["both", "narrow"] else conf.prefix + "_b_sort_peaks.broadPeak",
                                    output=conf.prefix + "_peaks_top_reg.bed",
                                    param={"peaks": 10000},
                                    name="top summits for regpotential"))
    get_top_peaks.update(param=conf.items("reg_potential"))
    get_top_peaks.allow_fail = True
    get_top_peaks.allow_dangling = True

    reg = attach_back(workflow,
                      ShellCommand(
                          "{tool} -t {input[peaks]} -g {param[geneTable]} -n {param[prefix]} -d {param[dist]}",
                          tool = "RegPotential.py",
                          input = {"peaks": conf.prefix + "_peaks_top_reg.bed"},
                          output = {"potential": conf.prefix + "_gene_score.txt"},
                          param = {"geneTable": conf.get_path(conf.get("basics", "species"), "geneTable"),
                                   "tool": resource_filename("chilin2.modules", "regulatory/RegPotential.py"),
                                   "prefix": conf.prefix + "_gene_score.txt",
                                   "dist": 100000},
                          name = "Regulatory Potential"))
    reg.update(param=conf.items("reg_potential"))
    reg.allow_fail = True
    reg.allow_dangling = True

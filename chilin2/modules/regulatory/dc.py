

## regulatory potential

from samflow.command import ShellCommand
from samflow.workflow import attach_back
from pkg_resources import resource_filename

def reg_potential(workflow, conf):
    """

    """

    ## exclude some peaks regions, such as velcro

    ##
    get_top_peaks = attach_back(workflow,
                                ShellCommand(
                                    "{tool} -n {param[peaks]} {input} | cut -f 1,2,3,4,9> {output}",
                                    tool="head",
                                    input=conf.prefix + "_sort_peaks.narrowPeak" if conf.get("macs2", "type") in ["both", "narrow"] else conf.prefix + "_b_sort_peaks.broadPeak",
                                    output=conf.prefix + "_peaks_top_reg.bed",
                                    param={"peaks": 10000},
                                    name="top summits for regpotential"))
    get_top_peaks.update(param=conf.items("reg_potential"))

    reg = attach_back(workflow,
                      ShellCommand(
                          "{tool} -t {input[peaks]} -g {param[table]} -n {param[prefix]} -d {param[dist]}",
                          tool = "RegPotential.py",
                          input = {"peaks": conf.prefix + "_peaks_top_reg.bed"},
                          output = {"potential": conf.prefix + "_gene_score.txt"},
                          param = {"table": conf.get_path(conf.get("basics", "species"), "regpotential"),
                                   "tool": resource_filename("chilin2.modules", "regulatory/RegPotential.py"),
                                   "prefix": conf.prefix,
                                   "dist": 10000},
                          name = "Regulatory Potential"))
    reg.update(param=conf.items("reg_potential"))

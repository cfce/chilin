from samflow.command import ShellCommand
from samflow.workflow import attach_back
from chilin2.modules.library.qc import stat_pbc


def PBC(workflow, conf):  # PBC1
    """
    Introduce ENCODE II library complexity assessment methods
    N1 / Nd, N1 is the location with exact one read, Nd is distinct location number
    :param workflow: samflow class
    :param conf: parsed config
    :return: void
    """
    for t in conf.sample_targets:
        pbc1 = attach_back(workflow,
                           ShellCommand(
                               """
                               bamToBed -i {input[bam]} | {tool} \'{{l[$1"\\t"$2"\\t"$3"\\t"$6]+=1}} END {{for(i in l) print l[i]}}\' \\
                                 | awk \'{{n[$1]+=1}} END {{for (i in n) print i"\\t"n[i]}}\'  \\
                                 | sort -k1n -  > {output[hist]}
                               awk '{{
                               if (NR==1) {{N1=$2}}
                               Nd+=$2
                               }} END {{print N1,Nd,N1/Nd}}' {output[hist]} > {output[pbc]}
                               """,
                               tool = "awk",
                               input = {"bam":t + "_4000000.bam" if conf.down else t + ".bam"},
                               output = {"pbc": t + ".pbc",
                                         "hist": t + ".hist"},
                               name = "PBC"))
        pbc1.allow_fail = True
        pbc1.allow_dangling = True

    ## QC part
    stat_pbc(workflow, conf)

from samflow.command import ShellCommand, PythonCommand
from samflow.workflow import attach_back
from chilin2.modules.fastqc.qc import stat_fastqc
from chilin2.modules.fastqc.tex import tex_fastqc

def fastqc(workflow, conf):
    """
    fastqc to extract gc contents(not yet) and median sequence quality
    :param workflow:
    :param conf:
    :return:
    """
    for raw, target in conf.sample_pairs:
        if conf.pe:
            fastqc_run = attach_back(workflow,
                                     ShellCommand(
                                         "{tool} {input} --extract -t {param[threads]} -o {output[target_dir]}",
                                         ## only check one pair
                                         input = target[0] + "_100k.fastq",
                                         output = {"target_dir": conf.target_dir,
                                                   "fastqc_summary": target[0] + "_100k_fastqc/fastqc_data.txt"},
                                         tool = "fastqc",
                                         param = {"threads": conf.threads},
                                         name = "fastqc"))
        else:
            fastqc_run = attach_back(workflow,
                                     ShellCommand(
                                         "{tool} {input} --extract -t {param[threads]} -o {output[target_dir]}",
                                         input=target + "_100k.fastq",
                                         output={"target_dir": conf.target_dir,
                                                 "fastqc_summary": target + "_100k_fastqc/fastqc_data.txt"},
                                         tool="fastqc",
                                         param={"threads": conf.threads},
                                         name = "fastqc"))
            fastqc_run.update(param=conf.items("fastqc"))
        fastqc.allow_fail = True
        fastqc.allow_dangling = True

    ## QC part of chilin
    ## use conf property conf.long = True
    stat_fastqc(workflow, conf)
    if conf.long:
        tex_fastqc(workflow, conf)

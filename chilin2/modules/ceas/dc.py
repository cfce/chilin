"""
DC version of ceas
"""
from samflow.command import ShellCommand, PythonCommand
from samflow.workflow import attach_back
from chilin2.modules.ceas.qc import stat_ceas, stat_bedAnnotate


def bedtools_ceas(workflow, conf):    # use bedtools
    """
    ###########################################################################
    DEPRECATED!! Please use bedAnnotate below
    ###########################################################################

    statistics of meta gene distribution,
    1bp overlap as cutoff
    :return:
    """
    summits = conf.prefix + "_sort_summits.bed" if conf.get("macs2", "type") in ["both", "narrow"] else conf.prefix + "_b_sort_peaks.broadPeak"
    ceas = attach_back(workflow, ShellCommand(
                       """
                       n=$(head -n {param[peaks]} {input[summits]} | wc -l)
                       head -n {param[peaks]} {input[summits]} > {input[summits]}.tmp
                       if [ {param[type]} == summit ]
                       then
                           exon=$({tool} intersect -f {param[p]} -wa -u -a {input[summits]}.tmp -b {param[exon]} | wc -l)
                           intron=$({tool} intersect -f {param[p]} -wa -u -a {input[summits]}.tmp -b {param[intron]} | wc -l)
                           intergenic=$({tool} intersect -f {param[p]} -wa -u -a {input[summits]}.tmp -b {param[intergenic]} | wc -l)
                           promoter=$({tool} intersect -f {param[p]} -wa -u -a {input[summits]}.tmp -b {param[promotor]} | wc -l)
                           exon=$(echo \"scale=5;$exon/$n\" | bc)
                           intron=$(echo \"scale=5;$intron/$n\" | bc)
                           promoter=$(echo \"scale=5;$promoter/$n\" | bc)
## we assume the complementary set of exon, intron and promoter(2kb) to be intergenic
                           intergenic=$(echo \"scale=5;$intergenic/$n\" | bc)
                       else
                           exon=$({tool} intersect -f {param[p]} -wo -a {input[summits]}.tmp -b {param[exon]} | awk '{{n+=$7}}END{{print n}}')
                           intron=$({tool} intersect -f {param[p]} -wo -a {input[summits]}.tmp -b {param[intron]} | awk '{{n+=$7}}END{{print n}}')
                           intergenic=$({tool} intersect -f {param[p]} -wo -a {input[summits]}.tmp -b {param[intergenic]} |awk '{{n+=$7}}END{{print n}}')
                           promoter=$({tool} intersect -f {param[p]} -wo -a {input[summits]}.tmp -b {param[promotor]} |awk '{{n+=$7}}END{{print n}}')
                           n=$(echo \"scale=2;$exon+$intron+$intergenic+$promoter\" | bc)
                           exon=$(echo \"scale=5;$exon/$n\" | bc)
                           intron=$(echo \"scale=5;$intron/$n\" | bc)
                           intergenic=$(echo \"scale=5;$intergenic/$n\" | bc)
                           promoter=$(echo \"scale=5;$promoter/$n\" | bc)
                       fi
                       echo $exon,$intron,$intergenic,$promoter > {output[meta]}
                       """,
                       tool = "bedtools",
                       input = {"summits": summits},
                       output = {"meta": conf.prefix + ".meta"},
                       param = {"p": "1E-9",## 0.3, 30% too strict
                                "peaks": 5000,
                                "type": "summit" if conf.get("macs2", "type") in ["both", "narrow"] else "base",
                                "exon": conf.get_path(conf.get("basics", "species"), "ceas_exon"),
                                "intron": conf.get_path(conf.get("basics", "species"), "ceas_intron"),
                                "intergenic": conf.get_path(conf.get("basics", "species"), "ceas_intergenic"),
                                "promotor": conf.get_path(conf.get("basics", "species"), "ceas_promotor")},
                       name = "bedtools ceas"))

    ## DHS, velcro, meta gene statistics
    ## DHS, velcro judgements
    try:
        has_velcro = conf.get(conf.get("basics", "species"), "velcro")
        has_dhs = conf.get(conf.get("basics", "species"), "dhs")
    except:
        has_velcro = ""
        has_dhs = ""
    if has_dhs:
        DHS(workflow, conf)
    if has_velcro:
        velcro(workflow, conf)
    stat_ceas(workflow, conf, has_dhs, has_velcro)

#def stat_bedAnnotate(workflow, conf, has_dhs, has_velcro):
#    """collect the stats for 
def bedAnnotate_ceas(workflow, conf):
    """
    Calls bedAnnotate to get the genome distribution of the summits
    """
    import os
    summits = conf.prefix + "_sort_summits.bed" if conf.get("macs2", "type") in ["both", "narrow"] else conf.prefix + "_b_sort_peaks.broadPeak"
    ceas = attach_back(workflow, ShellCommand(
            """{tool} -g {param[geneTable]} -b {input} -e {output[exon]} -t {output[gene]}> {output[meta]}
            meta_info.sh {output[gene]} {output[exon]} 2000 {param[chrominfo]}
            """,
            tool="bedAnnotate.py",
            input=summits,
            output={"meta":conf.prefix + ".meta", "gene":os.path.join(conf.target_dir, "gene.bed"), "exon": os.path.join(conf.target_dir, "exon.bed"),
                    "promoter": os.path.join(conf.target_dir, "gene.bed_promoter"), "exon": os.path.join(conf.target_dir, "gene.bed_exon")},
            param={"geneTable": conf.get_path(conf.get("basics", "species"), "geneTable"),
                   "chrominfo": conf.get_path(conf.get("basics", "species"), "chrom_len")},
            name="bedAnnotate (ceas)"))
    try:
        has_velcro = conf.get(conf.get("basics", "species"), "velcro")
        has_dhs = conf.get(conf.get("basics", "species"), "dhs")
    except:
        has_velcro = ""
        has_dhs = ""
    if has_dhs:
        DHS(workflow, conf)
    if has_velcro:
        velcro(workflow, conf)
    stat_bedAnnotate(workflow, conf, has_dhs, has_velcro)


def DHS(workflow, conf):   # DHS overlap percentage
    """
    get peaks overlapping percentage with union DHS
    :param workflow: uniform pipeline workflow from samflow
    :param conf: parsed config files
    :return: workflow
    """
    peaks = conf.prefix + "_sort_peaks.narrowPeak" if conf.get("macs2", "type") in ["both", "narrow"] else conf.prefix + "_b_sort_peaks.broadPeak"
    DHS = attach_back(workflow,
                      ShellCommand("""
                                   n=$(head -n {param[p]} {input[MACS2_bed]} | wc -l)
                                   dhs=$(head -n {param[p]} {input[MACS2_bed]} | {tool} -wa -u -a - -b {input[DHS_peaks_bed]}|wc -l)
                                   ##dhs=$(echo \"scale=5;$dhs/$n\" | bc)
                                   echo $n,$dhs > {output}
                                   """,
                                   tool="intersectBed",
                                   input={"MACS2_bed": peaks,
                                          "DHS_peaks_bed": conf.get(conf.get("basics", "species"), "dhs")},
                                   output=conf.prefix + ".dhs",
                                   param={"p": 5000},
                                   name = "intersect DHS"))


def velcro(workflow, conf):
    attach_back(workflow,
                ShellCommand(
                    """
                    n=$(head -n {param[p]} {input[MACS2_bed]} | wc -l)
                    velcro=$(head -n {param[p]} {input[MACS2_bed]} | {tool} -wa -u -a - -b {input[velcro_peaks_bed]} | wc -l)
                    velcro=$(echo \"scale=5;$velcro/$n\" | bc)
                    echo $velcro > {output}
                    """,
                    tool="intersectBed",
                    input={"MACS2_bed": conf.prefix + "_sort_peaks.narrowPeak" if conf.get("macs2", "type") in ["both", "narrow"]
                            else conf.prefix + "_b_sort_peaks.broadPeak",
                           "velcro_peaks_bed": conf.get(conf.get("basics", "species"), "velcro")},
                    output=conf.prefix + ".velcro",
                    param={"p": 5000},
                    name = "velcro overlap"))


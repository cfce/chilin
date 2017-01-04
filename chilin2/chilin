#!/usr/bin/env python
"""
a ChIP-seq pipeline for benchmarking quality control of the experiment
"""

import sys
import os
import argparse
from string import Template
import time

from chilin2.modules.config.helpers import write_conf
from samflow.command import ShellCommand
from samflow.workflow import Workflow, attach_back
from pkg_resources import resource_filename

from chilin2.modules.config.config import ChiLinConfig
from ConfigParser import SafeConfigParser

from chilin2.modules.bwa.dc import bwa
from chilin2.modules.bowtie.dc import bowtie
from chilin2.modules.star import star
from chilin2.modules.interface.dc import groom_sequencing_files
from chilin2.modules.ceas.dc import bedAnnotate_ceas  # bedtools_ceas
from chilin2.modules.conservation.dc import conservation
from chilin2.modules.contamination.dc import contamination_check
from chilin2.modules.fastqc.dc import fastqc
from chilin2.modules.frip.dc import FRiP
from chilin2.modules.library.dc import PBC
from chilin2.modules.macs.dc import macs2_rep, macs2
from chilin2.modules.macs2_fragment.dc import fragment
from chilin2.modules.mdseqpos.dc import seqpos
from chilin2.modules.enrichment.dc import read_enrichment_on_meta
# from chilin2.modules.phantompeak.dc import Phan
from chilin2.modules.replicates.dc import replicates
from chilin2.modules.regulatory.dc import reg_potential
from chilin2.modules.interface.dc import merge_bams
from chilin2.modules.interface.dc import sampling_bam
from chilin2.modules.summary.qc_summary_table import render_pdf


class FriendlyArgumentParser(argparse.ArgumentParser):
    """
    Override argparse to show full length help information
    """

    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(1)


def parse_args(args=None):
    """
    parse input arguments
    """
    parser = FriendlyArgumentParser(description=__doc__)
    parser.add_argument('-V', '--version', action='version',
                        version="%(prog)s (code version 2.0, db version 1.0)")
    sub_parsers = parser.add_subparsers(help="sub-command help", dest="sub_command")

    parser_run = sub_parsers.add_parser("run", help="run pipeline using a config file",
                                        description="ChiLin-run: Run ChiLin pipeline using a config file")
    parser_run.add_argument("-c", "--config", dest="config", required=True,
                            help="specify the config file to use")

    ## add gen parser
    parser_gen = sub_parsers.add_parser("gen", help="generate config files")
    parser_gen.add_argument("-o", "--output", dest="out", help="name of config file, default: chilin.conf",
                            default="chilin.conf")

    ## add simple mode for chilin
    # generate a species list
    cf = SafeConfigParser()
    cf.read(resource_filename('chilin2.modules', 'config/chilin.conf'))
    sections_ls = ['basics', 'tool', 'bwa', 'macs2', 'reg', 'conservation',
                   'seqpos', 'qc', 'contamination']
    species_ls = [s for s in cf.sections() if s not in sections_ls]

    parser_simple = sub_parsers.add_parser("simple",
                                           help="run simple modes, should prepare species conf files when installing")
    parser_simple.add_argument("-t", "--treat", dest="treat", required=True,
                               help="treatment data, comma(,) separate replicates, format(fastq,bam) fixed by extension")
    parser_simple.add_argument("-c", "--cont", dest='cont', default='',
                               help="control data, comma(,) separate replicates , format(fastq,bam) fixed by extension")
    parser_simple.add_argument("-s", "--species", dest="species", help="species--from chilin.conf: %s" % species_ls,
                               required=True)
    parser_simple.add_argument("-u", "--user", dest="user", default="anonymous",
                               help="your name, DEFAULT: anonymous")

    parser_simple.add_argument("-r", "--rtype", dest="rtype", help="input factor type for conservation", default="tf",
                               choices=["tf", "histone", "dnase"])
    parser_simple.add_argument("-p", "--type", dest="ptype",
                               help="input peaks type for macs2 peaks calling, narrow, broad or both", default="both",
                               choices=["narrow", "both", "broad"])
    parser_simple.add_argument("-o", "--output", dest="out", help="output directory name", default="output")
    parser_simple.add_argument("-i", "--id", dest="id",
                               help="name for your chilin run, suggestion: use your factor name, e.g. ER or MYC",
                               default="run")

    parser_batch = sub_parsers.add_parser("batch", help="run pipeline for multiple datasets")
    parser_batch.add_argument("-b", "--batch-config", dest="batch_config", required=True, help="batch file")

    ## batch mode, conf run mode, simple run mode
    for p in (parser_run, parser_batch, parser_simple):
        p.add_argument("-v", "--verbose-level", dest="verbose_level", type=int,
                       default=2)
        p.add_argument("--short", dest="long", action="store_false",
                       default=True,
                       help="Flag to output a short report--default: False")
        p.add_argument("--from", dest="start_step", default=0, type=int,
                       help="Only step after this number will be processed")
        p.add_argument("--to", dest="end_step", default=100, type=int,
                       help="Only step before this number will be processed ")
        p.add_argument("--skip", dest="skip_step", default="",
                       help="Steps to skip, use comma as seperator")
        p.add_argument("--mapper", dest="mapper", default="bwa",
                       choices=["bwa", "bowtie", "star"],
                       help="choose your mapper")  ## should be consistent with the genome index you specified in the gene
        ## new option for use total reads or down sampling 4M
        p.add_argument("--total", dest="down", action="store_false",
                       default=True,
                       help="default: 4M, down sampling 4M or use total reads, only support PBC, NSC/RSC/Qtag, FRiP")
        p.add_argument("--frip", dest="frip", action="store_true",
                       default=False, help="default: 4M for FRiP, decide using 5M or 4M for FRiP calculation")
        #        p.add_argument("--unsc", dest="unsc", action="store_true",
        #                       default=False, help="default: False, this option can be used with --total, decide using total non redundant uniquely mapped reads NSC or total uniquely mapped reads NSC.")

        p.add_argument("--dry-run", dest="dry_run", action="store_true",
                       default=False)
        p.add_argument("--not-allow-dangling", dest="allow_dangling", action="store_false",
                       default=True)
        p.add_argument("--not-allow-fail", dest="allow_fail", action="store_false",
                       default=True)
        p.add_argument("--dont_resume", dest="resume", action="store_false",
                       default=True)
        p.add_argument("--dont_remove", dest="clean", action="store_false",
                       default=True)
        p.add_argument("--threads", dest="threads", default=4, type=int,
                       help="threads # to use when aligning (default: 4)")
        p.add_argument("--pe", dest="pe", action="store_true", default=False,
                       help="default single end mode, turn on for pair end sequencing")

    return parser.parse_args(args), parser


class StepChecker:
    """
    For step control and selection
    """

    def __init__(self, start, end, skips):
        self.start = start
        self.end = end
        self.skips = skips

    def need_run(self, step_id):
        if step_id < self.start:
            return False
        if step_id > self.end:
            return False
        if step_id in self.skips:
            return False
        return True


class ChiLinBuilder:
    def __init__(self, workflow, conf):
        self.workflow = workflow
        self.conf = conf
        self.LaTex_fragments = []
        self.plain_fragments = []
        self.finished = set()

    def build(self, prep_func, tag=None):
        prep_func(self.workflow, self.conf)
        if tag:
            self.finished.add(tag)

    def attach_back(self, command):
        attach_back(self.workflow, command)


## for clean up work directory
def clean_up(workflow, conf):
    """
    package all the necessary results and delete temporary files
    """
    preserved = []
    not_moved = []
    added_bed = []

    for target in conf.sample_targets:
        preserved.append(target + ".bam")
    for target in conf.treatment_targets:
        preserved.append(target + "_peaks.bed")
        not_moved.append(target + "_peaks.bed")
        preserved.append(target + "_b_peaks.bed")
        not_moved.append(target + "_b_peaks.bed")
        preserved.append(target + "_treat.bw")
        not_moved.append(target + "_treat.bw")
        preserved.append(target + "_control.bw")
        not_moved.append(target + "_control.bw")
        preserved.append(target + "_sort_peaks.narrowPeak")
        not_moved.append(target + "_sort_peaks.narrowPeak")
        added_bed.append(target + "_sort_peaks.narrowPeak")
        # remove bedgraph
        #preserved.append(target + "_treat_pileup.bdg")
        #preserved.append(target + "_control_lambda.bdg")
        preserved.append(target + "_100k_fastqc")
        preserved.append(target + "_peaks.xls")
        not_moved.append(target + "_peaks.xls")
        preserved.append(target + "_sort_summits.bed")
        not_moved.append(target + "_sort_summits.bed")
        preserved.append(target + "_b_peaks.xls")
        not_moved.append(target + "_b_peaks.xls")
        preserved.append(target + "_b_sort_peaks.broadPeak")
        not_moved.append(target + "_b_sort_peaks.broadPeak")
        added_bed.append(target + "_b_sort_peaks.broadPeak")
    if len(conf.treatment_targets) >= 2:
        preserved.append(conf.prefix + "_treatment.bam")
    if len(conf.control_targets) >= 2:
        preserved.append(conf.prefix + "_control.bam")
    preserved.append(conf.prefix + "_treat.bw")
    not_moved.append(conf.prefix + "_treat.bw")
    preserved.append(conf.prefix + "_control.bw")
    not_moved.append(conf.prefix + "_control.bw")
    preserved.append(conf.prefix + "_conserv_img.png")
    preserved.append(conf.prefix + "_conserv_compare.pdf")
    preserved.append(conf.prefix + "_sort_peaks.narrowPeak")
    not_moved.append(conf.prefix + "_sort_peaks.narrowPeak")
    added_bed.append(conf.prefix + "_sort_peaks.narrowPeak")
    preserved.append(conf.prefix + "_peaks.xls")
    not_moved.append(conf.prefix + "_peaks.xls")
    preserved.append(conf.prefix + "_peaks.bed")
    not_moved.append(conf.prefix + "_peaks.bed")
    preserved.append(conf.prefix + "_sort_summits.bed")
    not_moved.append(conf.prefix + "_sort_summits.bed")
    preserved.append(conf.prefix + "_b_peaks.xls")
    not_moved.append(conf.prefix + "_b_peaks.xls")
    preserved.append(conf.prefix + "_b_sort_peaks.broadPeak")
    not_moved.append(conf.prefix + "_b_sort_peaks.broadPeak")
    added_bed.append(conf.prefix + "_b_sort_peaks.broadPeak")
    preserved.append(conf.prefix + "_seqpos")
    preserved.append(conf.prefix + "_gene_score.txt")
    preserved.append(conf.prefix + "_gene_score_5fold.txt")
    preserved.append(conf.prefix + "_conserv.txt")

    # remove bedgraph
    #preserved.append(conf.prefix + "_b_treat_pileup.bdg")
    #preserved.append(conf.prefix + "_b_control_lambda.bdg")
    #preserved.append(conf.prefix + "_treat_pileup.bdg")
    #preserved.append(conf.prefix + "_control_lambda.bdg")
    preserved.append(conf.prefix + ".tex")
    preserved.append(os.path.join(conf.target_dir, "json"))
    preserved.append(conf.prefix + ".pdf")
    not_moved.append(conf.prefix + ".pdf")
    d_pattern = map(lambda x: os.path.join(conf.target_dir, x), os.listdir(conf.target_dir))

    attach_back(workflow, ShellCommand(
        "if [ ! -d '{output}' ]; then mkdir -p {output}; fi",
        output=os.path.join(conf.target_dir, "attic")))

    for df in d_pattern:
        if not df in preserved:
            deleted = attach_back(workflow,
                                  ShellCommand("rm -rf {param[delete]}",
                                               param={"delete": df}))
            deleted.allow_fail = True

    for pf in preserved:
        if not pf in not_moved:
            mv = attach_back(workflow,
                             ShellCommand("mv {param[move]} {param[attic]} >/dev/null 2>&1",
                                          param={"move": pf, "attic": os.path.join(conf.target_dir, "attic")}))
            mv.allow_fail = True

    for af in added_bed:
        add_bed_suffix = attach_back(workflow,
                                     ShellCommand("mv {param[file1]} {param[file2]} > /dev/null 2>&1",
                                                  param={"file1": af, "file2": af + ".bed"}))
        add_bed_suffix.allow_fail = True


def create_cleanup(args, conf):
    """creates a cleanup workflow--invoked AFTER a successful run"""
    workflow = Workflow(name="ChiLin2 Cleanup")
    workflow.set_option(
        verbose_level=args.verbose_level,
        dry_run_mode=args.dry_run,
        resume=args.resume,
        allow_dangling=args.allow_dangling)
    bld = ChiLinBuilder(workflow, conf)
    bld.build(clean_up)
    return workflow


def create_workflow(args, conf):
    """
    args is the sys.argv
    conf is the config files
    :param args: sub_command(sub_parser)
    :param conf: input config files
    :return: workflow
    """
    workflow = Workflow(name="ChiLin Main")
    bld = ChiLinBuilder(workflow, conf)

    if args.skip_step:
        skipped_steps = [int(i) for i in args.skip_step.split(",")]
    else:
        skipped_steps = []

    step_checker = StepChecker(args.start_step, args.end_step, skipped_steps)

    ## flexible options for replicates and species annotations
    have_treat_reps = len(conf.treatment_pairs) >= 2
    has_conservation = conf.get(conf.get("basics", "species"), "conservation")
    has_motifdb = conf.get("seqpos", "db")

    need_run = step_checker.need_run

    bld.attach_back(ShellCommand(
        "if [ ! -d '{output}' ]; then mkdir -p {output}; fi",
        output=conf.target_dir))

    bld.build(write_conf)
    ## DC prepare sampling reads and raw reads fastqs, support Fastq currently
    bld.build(groom_sequencing_files)

    # -------------------------
    ## need input fastq
    # A. reads qc
    if need_run(1):
        ## fastqc on 100k reads
        bld.build(fastqc)
        bld.build(contamination_check)

    if need_run(2):  ## this step can not be skipped, provide bam files
        if conf.mapper == "bwa":
            bld.build(bwa)
        if conf.mapper == "bowtie":
            bld.build(bowtie)
        if conf.mapper == "star":
            bld.build(star)
        ## merge bams
        bld.build(merge_bams)
    #----------------------------
    # B. ChIP qc
    if need_run(3):
        bld.build(sampling_bam)  ## provide sam and bam
        bld.build(fragment)

    if need_run(4):
        bld.build(sampling_bam)  ## provide sam and bam
        bld.build(PBC)

    if need_run(5):  ## if no replicates, skip automatically
        if have_treat_reps:
            bld.build(macs2_rep)
            bld.build(replicates)
    if need_run(6):  ## this step can not skipped if run following steps
        bld.build(macs2)

        ## following all needs merged peaks or summits
        if need_run(7):
            ## merged peaks FRiP
            bld.build(FRiP)  ## sample and IP qc
            # --------------------
            # C. annotation qc
        if need_run(8):
            bld.build(bedAnnotate_ceas)
            ## generate gene.bed, exon.bed for read_enrichment analysis
            if need_run(9):
                bld.build(sampling_bam)  ## provide sam and bam
                bld.build(read_enrichment_on_meta)
                ## following steps dependent on step 6th
        if need_run(10):
            if has_conservation:
                ## wig from ucsc and convert to bigwiggle
                ## make this more common
                bld.build(conservation)
        if need_run(11):
            bld.build(reg_potential)
        if need_run(12):
            if has_motifdb:
                bld.build(seqpos)
        if need_run(13):
            bld.build(render_pdf)

    return workflow


def main(args=None):
    args, parser = parse_args(args)

    if args.sub_command == "run":
        conf = ChiLinConfig(args.config, args)
        # long document or only summary
        conf.long = args.long
        conf.frip = args.frip
        conf.down = args.down
        conf.mapper = args.mapper

        ## edit macs2 section species part to
        ## your specified species effective genome size
        # handle threads
        conf.threads = args.threads
        workflow = create_workflow(args, conf)

        workflow.set_option(
            verbose_level=args.verbose_level,
            dry_run_mode=args.dry_run,
            resume=args.resume,
            logger=conf.log,
            allow_dangling=args.allow_dangling,
            allow_fail=args.allow_fail)

        if workflow.invoke():  #invoke returns True on success
            if args.clean:  #clean_up
                wkflw = create_cleanup(args, conf)
                wkflw.invoke()

    if args.sub_command == "gen":
        fin = open(resource_filename("chilin2.modules", "config/chilin.conf"))
        t = Template(fin.read())
        fin.close()

        # overwrite basics
        ls = ['USER', 'ID', 'TIME', 'SPECIES', 'FACTOR', 'TREAT', 'CONT', 'OUTPUT']
        d = dict([(fld, '') for fld in ls])
        d['TIME'] = time.strftime("%Y-%m-%d")
        fout = open(args.out, "w")
        fout.write(t.safe_substitute(d))
        fout.close()

    if args.sub_command == "simple":
        conf = ChiLinConfig(resource_filename("chilin2.modules", "config/chilin.conf"), args)  ## pass command parameter
        conf.frip = args.frip
        conf.down = args.down
        conf.long = args.long
        conf.mapper = args.mapper

        conf.root_dir = os.path.dirname(args.treat.split(",")[0])  ## input and output should be in the same directory
        ## long document or only summary
        #handle threads
        conf.threads = args.threads
        conf.set("basics", "treat", args.treat)
        conf.set("basics", "time", time.strftime("%Y-%m-%d"))
        if args.cont:
            conf.set("basics", "cont", args.cont)
        else:
            conf.set("basics", "cont", "")
        conf.set("basics", "species", args.species)

        if args.species in ["hg18", "hg19", "hg38"]:
            conf.set("macs2", "species", "hs")
        if args.species in ["mm8", "mm9", "mm10"]:
            conf.set("macs2", "species", "mm")

        conf.set("seqpos", "species", args.species)

        conf.set("conservation", "type", args.rtype)

        conf.set("basics", "id", args.id)
        conf.set("basics", "output", os.path.abspath(args.out))
        conf.set("basics", "factor", args.rtype)
        conf.set("basics", "user", args.user)
        conf.set("macs2", "type", args.ptype)

        ## uniform extension for different types of factors
        if args.rtype.lower() in ['tf', 'histone']:
            conf.set("macs2", "extsize", "146")  # for histone and tf, use 146 extsize(73 shiftsize)
        elif args.rtype.lower() == "dnase":
            conf.set("macs2", "extsize", "100")  # for dnase, use 100 extsize(50 shiftsize)
        else:
            print >> sys.stderr, "No such factor\n"

        workflow = create_workflow(args, conf)
        workflow.set_option(
            verbose_level=args.verbose_level,
            dry_run_mode=args.dry_run,
            resume=args.resume,
            logger=conf.log,
            allow_dangling=args.allow_dangling,
            allow_fail=args.allow_fail)

        if workflow.invoke():  #invoke returns True on success
            if args.clean:  #clean_up
                wkflw = create_cleanup(args, conf)
                wkflw.invoke()

    if args.sub_command == "batch":
        with open(args.batch_config) as batch_file:
            for a_conf in batch_file:
                a_conf = a_conf.strip()
                conf = ChiLinConfig(a_conf, args)
                conf.frip = args.frip
                conf.down = args.down
                conf.mapper = args.mapper
                conf.long = args.long
                # conf.unsc = args.unsc
                conf.threads = args.threads
                workflow = create_workflow(args, conf)
                workflow.set_option(
                    verbose_level=args.verbose_level,
                    dry_run_mode=args.dry_run,
                    resume=args.resume,
                    allow_dangling=args.allow_dangling,
                    allow_fail=args.allow_fail)

                ## one by one clean up
                if workflow.invoke():  #invoke returns True on success
                    if args.clean:  #clean_up
                        wkflw = create_cleanup(args, conf)
                        wkflw.invoke()


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr("User interrupt:) \n")

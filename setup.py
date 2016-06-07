"""
Time-stamp: <2014-11-13 09:38:31 qqin>
install ChiLin and related program and dependent data
"""
from pkg_resources import resource_filename
import os
from distribute_setup import use_setuptools
import platform
import argparse
from argparse import RawDescriptionHelpFormatter
from pprint import pprint
import sys
import subprocess

use_setuptools()

from setuptools import setup
from ConfigParser import SafeConfigParser

def which(program):
    import os

    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = os.path.expanduser(path.strip('"'))
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return 'Missing'

def module_list():
    dependency = [
    ['--------', '------'],
    ['software', 'status'],
    ['--------', '------'],
    ['seqtk', which('seqtk')], ## default bwa
    ['fastqc', which('fastqc')],
    ['bwa', which('bwa')], ## default bwa
    ['samtools', which('samtools')],
    ['macs2', which('macs2')],
    ['bedtools', which('bedtools')],
    ['bedClip', which('bedClip')],
    ['bedGraphToBigWig', which('bedGraphToBigWig')],
    ['wigCorrelate', which('wigCorrelate')],
    ['wigToBigWig', which('wigToBigWig')],
    ]
    try:
        from bx.bbi.bigwig_file import BigWigFile
        dependency.append(['bx-python', BigWigFile])
    except:
        dependency.append(['bx-python', 'Missing'])

    dependency+=[['MDSeqPos.py', which('MDSeqPos.py')]]
    dependency += [
    ['--------', '------'],
    ['data', 'status'],
    ['--------', '------']]
    cf = SafeConfigParser()
    cf.read(resource_filename('chilin2.modules', 'config/chilin.conf'))
    sections_ls = ['basics', 'tool', 'bwa', 'macs2', 'reg', 'conservation',
                   'seqpos', 'qc']
    species_ls = [s for s in cf.sections() if s not in sections_ls]

    for sp in species_ls:
        options = cf.options(sp)
        for op in options:
            if os.path.exists(cf.get(sp, op)) or os.path.exists(cf.get(sp, op)+".amb"):
               status = "ok at " + '  ' + cf.get(sp, op)
            else:
               status = "missing"
            dependency += [[sp+'  '+op, status]]
            
    for d, s in dependency:
        pprint('{: ^25}: {: ^50}'.format(d, s))

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

def generate_defaults():
    """Tries to generate the default run configuration files from the
    paths specified in chilin.conf.  The idea is that an admin, who is
    installing ChiLin system-wide could define the defaults which will allow
    users to avoid generating/editing conf files!!
    """

    cf = SafeConfigParser()
    cf.add_section('basics')
    ls = ['user', 'id', 'time', 'species', 'factor', 'treat', 'cont', 'output']
    for fld in ls:
        cf.set('basics', fld, '$'+fld.upper())
    #SET the chilin version number
    cf.set('basics', "version", _CHILIN_VERSION)

    #read in the chilin.conf file--look first for developer defaults
    if os.path.exists('chilin.conf.filled'):
        cf.read('chilin.conf.filled')
    else:
        cf.read('chilin.conf')

    #write the template file!
    f = open(os.path.join('chilin2','modules','config','chilin.conf'),'w')
    cf.write(f)
    f.close()

def default_species_section():
    """add default settings in chilin.conf.filled for hg38, hg19, mm9 and mm10
    """
    cf = SafeConfigParser()
    #read in the chilin.conf file--look first for developer defaults
    if os.path.exists('chilin.conf'):
        cf.read('chilin.conf')
    cf.add_section('basics')
    ls = ['user', 'id', 'time', 'species', 'factor', 'treat', 'cont', 'output']
    for fld in ls:
        cf.set('basics', fld, '$'+fld.upper())
    #SET the chilin version number
    cf.set('basics', "version", _CHILIN_VERSION)

    current_dir = os.getcwd()
    cf.set("contamination", "mycoplasma", os.path.join(current_dir, "db", "mycoplasma.fasta"))

    for sp in ['hg38', 'hg19', 'mm9', 'mm10']:
        cf.add_section(sp)
        cf.set(sp, "genome_index", os.path.join(current_dir, "db", sp, sp+".fa"))
        cf.set(sp, "genome_dir", os.path.join(current_dir, "db", sp))
        cf.set(sp, "genetable", os.path.join(current_dir, "db", sp, sp+".refGene"))
        cf.set(sp, "chrom_len", os.path.join(current_dir, "db", sp, "chromInfo_"+sp+".txt"))
        cf.set(sp, "dhs", os.path.join(current_dir, "db", sp, "DHS_"+sp+".bed"))
        if sp in ["mm9", "hg19"]:
           cf.set(sp, "conservation", os.path.join(current_dir, "db", sp, "placental"))
           cf.set(sp, "velcro", os.path.join(current_dir, "db", sp, sp+"-blacklist.bed"))
        elif sp in ["mm10", "hg38"]:
           cf.set(sp, "conservation", os.path.join(current_dir, "db", sp, "phastcon.bw"))
           cf.set(sp, "velcro", "")


    #write the template file!
    f = open('chilin.conf.filled','w')
    cf.write(f)
    f.close()
    return

def clean():
    """ clean up 
    """
    import os
    clean_up = os.system("""
cd software
cd mdseqpos && python setup.py clean --all && cd ..
for i in bedtools-2.17.0  bowtie  bwa samtools-0.1.19  seqtk
do
  cd $i && make clean && cd ..
done
cd ..
rm -rf build chilin2.egg-info chilin.iml dist distribute-0.6.49-py2.6.egg distribute-0.6.49-py2.7.egg distribute-0.6.49.tar.gz flycheck_setup.py TAGS
rm -rf chilin_env
""")
    return

def install():
    setup(
        name='chilin2',
        #LEN: please set the version number above
        version=_CHILIN_VERSION,
        packages=['chilin2', 'chilin2.modules', 'chilin2.modules.bwa',
                  'chilin2.modules.bowtie', 'chilin2.modules.ceas',
                  'chilin2.modules.config', 'chilin2.modules.conservation',
                  'chilin2.modules.enrichment', 'chilin2.modules.macs2_fragment',
                  'chilin2.modules.contamination', 'chilin2.modules.dbaccessor',
                  'chilin2.modules.fastqc', 'chilin2.modules.frip', 'chilin2.modules.interface',
                  'chilin2.modules.library', 'chilin2.modules.macs', 'chilin2.modules.mdseqpos',
                  ## 'chilin2.modules.phantompeak',
                  'chilin2.modules.regulatory',
                  'chilin2.modules.star',
                  'chilin2.modules.replicates', 'chilin2.modules.summary', 'chilin2.modules.washU',
                  'samflow'],

        package_dir = {'chilin2': 'chilin2',
                       'samflow': 'samflow'},
        package_data = {"chilin2.modules" : ["config/*.conf", "dbaccessor/ChiLinQC.db",
                                             "dbaccessor/*.txt", "summary/*tex", "summary/*R",
                                             "mdseqpos/*tex",
                                             "phantompeak/*tex", "bwa/*tex", "conservation/*tex", "contamination/*tex",
                                             "summary/*cls", "fastqc/*tex", "frip/*tex",
                                             "summary/CFCE_Logo_Final.jpg"]},
        url='http://cistrome.org/chilin/',
        license='MIT',
        platforms='linux/unix',
        author='Hanfei Sun, Shenglin Mei, Qian Qin, Len Taing',
        author_email='1410771@tongji.edu.cn',
        description=read("README.md"),
        scripts = ["chilin2/chilin", "chilin2/modules/conservation/conservation_plot.py",
		   "chilin2/modules/conservation/conservation_onebw_plot.py",
                   "chilin2/modules/ceas/bedAnnotate.py",
                   "chilin2/modules/ceas/meta_info.sh",
                   "chilin2/modules/interface/sampling_pe_sam.py",
                   "chilin2/modules/interface/filter_pe_sam_faster.py",
                   "chilin2/modules/regulatory/RegPotential.py"],
        install_requires=['jinja2','argparse','macs2','numpy','cython', 'bx-python'])


def install_full():
    if not os.path.exists("chilin_env"):
	    setup_env = subprocess.call(
	    """
	    python virtualenv.py -ppython2.7  chilin_env
	    python virtualenv.py -ppython2.7  chilin_env --relocatable
	    #python virtualenv.py -ppython2.7 --system-site-packages --distribute chilin_env
	    #python virtualenv.py -ppython2.7 --system-site-packages --distribute chilin_env --relocatable
	    """, shell=True)

    execfile("chilin_env/bin/activate_this.py", dict(__file__="chilin_env/bin/activate_this.py"))
    setup_env = subprocess.call(
        """
        source chilin_env/bin/activate
        cd software
        if [ ! -s chilin_env/bin/bedtools ] 
        then
		cd FastQC && chmod 755 * && cp -rf * ../../chilin_env/bin && cd ..
		cd bwa && make clean && make && cp bwa ../../chilin_env/bin && cd ..
		cd seqtk && make clean && make && cp seqtk ../../chilin_env/bin && cd ..
		cd samtools* && make clean && make && cp samtools ../../chilin_env/bin && cd ..
		cd bedtools* && make clean && make && cp bin/* ../../chilin_env/bin && cd ..
        fi
        """, shell=True)


    if platform.system() == 'Linux':
        setup_env = subprocess.check_call(
            """
            cp chilin2/chilin chilin_env/bin/
            . chilin_env/bin/activate
            cd software
            chmod 755 ucsc/linux/*
            cp -f ucsc/linux/* ../chilin_env/bin
            chmod -R 755 ../chilin_env/bin
            """, shell=True, executable='/bin/bash')
    else:
        setup_env = subprocess.check_call(
            """
            . chilin_env/bin/activate
            cd software
            chmod 755 ucsc/mac/*
            cp -f ucsc/mac/* ../chilin_env/bin
            chmod -R 755 ../chilin_env/bin
            """, shell=True, executable='/bin/bash')

    #execfile("chilin_env/bin/activate_this.py", dict(__file__="chilin_env/bin/activate_this.py"))
    #install()

    setup_env = subprocess.check_call(
        """
        source chilin_env/bin/activate
        python setup.py install
        cd software
        cd bx-python && python setup.py install && cd ..
        cd mdseqpos && python setup.py install && cd ..
        which R && R -e "source('http://bioconductor.org/biocLite.R');biocLite('seqLogo');library(seqLogo)"
        """, shell=True, executable='/bin/bash')


class FriendlyArgumentParser(argparse.ArgumentParser):
    """
    Override argparse to show full length help information
    """
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(1)

_CHILIN_VERSION="2.0.0"
def main():
    '''
    main setup script
    :return:
    '''
    #READ in chilin.conf and generate the default confs.
    p = FriendlyArgumentParser(description=__doc__, formatter_class=RawDescriptionHelpFormatter, usage='setup.py [list|install]')

    p.add_argument("install", nargs="?", help="install chilin")
    p.add_argument("--list", '-l', default=False, action="store_true", help="list chilin dependent software")
    p.add_argument("--full", '-f', default=False, action="store_true", help="install dependent software")
    args = p.parse_args()

    if not (args.list or args.install):
       p.print_help()
       sys.exit(1)
    if args.list:
        if os.path.exists("chilin_env/bin/activate_this.py"):
            execfile("chilin_env/bin/activate_this.py", dict(__file__="chilin_env/bin/activate_this.py"))
        module_list()
        sys.exit(0)

    if not args.install in ["install", "develop", "clean"]:
        pprint("maybe typo, use `install` to install or `develop` to develop or `-l` to list or clean to clean up")
        p.print_help()
        sys.exit(1)

    if args.install:
        if args.install == "clean":
            clean()
            return
        if os.path.exists('chilin.conf'):
            generate_defaults()
            if args.full:
                default_species_section()
                install_full()
            else:
                install()

if __name__ == '__main__':
    main()


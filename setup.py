"""
Time-stamp: <2014-01-16 10:46:09 qqin>
install ChiLin2 and related program and dependent data
"""
import os
from distribute_setup import use_setuptools
use_setuptools()
from setuptools import setup
from ConfigParser import SafeConfigParser


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



_CHILIN_VERSION="2.0.0"
def main():
    '''
    main setup script
    :return:
    '''
    #READ in chilin.conf and generate the default confs.
    if os.path.exists('chilin.conf'):
        generate_defaults()

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
                  'chilin2.modules.phantompeak', 'chilin2.modules.regulatory',
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
        author='Hanfei Sun, Shenglin Mei, Qian Qin, Len',
        author_email='qianqind@gmail.com',
        description=read("README.md"),
        scripts = ["chilin2/ChiLin2.py", "chilin2/modules/conservation/conservation_plot.py",
                   "chilin2/modules/ceas/bedAnnotate.py",
                   "chilin2/modules/interface/sampling_pe_sam.py",
                   "chilin2/modules/regulatory/RegPotential.py"],
        install_requires=['jinja2','argparse'])


if __name__ == '__main__':
    main()

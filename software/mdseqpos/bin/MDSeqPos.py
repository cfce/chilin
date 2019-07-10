#!/usr/bin/env python

"""
SYNOPSIS: This tool searches for known and new (denovo) motifs that are in
regions specified by a bed file.
"""

import os
import sys
import optparse
import time
import random
import shutil
import json
import numpy
import math
from datetime import datetime

from jinja2 import Template

import mdseqpos
from mdseqpos.motif import MotifList
from mdseqpos.chipregions import ChipRegions
import mdseqpos.settings as settings
import mdseqpos.bayesian_motif_comp as bmc
import mdseqpos.pwmclus_motif_comp as pmc

_DEBUG = False

#REMOVE _DEBUG as a parameter
def read_known_motifs(motif_dbs, _DEBUG = False):
    """Given a list of xml file names, this function tries to load the motifs
    in those databases
    """
    DATA_DIR = os.path.join(settings.DEPLOY_DIR, 'database')

    known_motifs = MotifList()
    for db in motif_dbs:
        if _DEBUG: print "loading (time): %s (%s)" % (db, time.ctime())
        tmp = MotifList()
        tmp.from_xml_file(os.path.join(DATA_DIR, db))
        known_motifs.extend(tmp)
        if _DEBUG: print "load Complete (time): %s (%s)" % (db, time.ctime())
    return known_motifs

def seqpos_cluster_known_motifs(known_motifs, chip_regions, width, cutoff):
    """Cluster the known motifs based on their similarity score,then one motif
    could represent the whole cluster
    """
    CLUSTER = os.path.join(settings.DEPLOY_DIR, 'database', 'cistrome.cluster')
    cluster_motifs = {} #cluster_id to motifs
    id2motif = {}
    for m in known_motifs:
        id2motif[m.id] = m
    for line in open(CLUSTER,'r').readlines():
        mid, cid = line.strip().split('\t')
        if mid in id2motif:
            if cid in cluster_motifs:
                cluster_motifs[cid].append(id2motif[mid])
            else:
                cluster_motifs[cid] = [id2motif[mid]]
    
    fitered_motifs = []
    for motifs in cluster_motifs.values():
        m0 = motifs[0]
        m0.seqpos(chip_regions, width)
        if m0.seqpos_results['pvalue'] <= cutoff:
            fitered_motifs.append(m0)
            for m in motifs[1:]:
                m.seqpos(chip_regions, width)
                fitered_motifs.append(m)
    return MotifList(fitered_motifs)

def sample(p):
    """Given an array of probabilities, which sum to 1.0, randomly choose a 'bin',
    e.g. if p = [0.25, 0.75], sample returns 1 75% of the time, and 0 25%;
    NOTE: p in this program represents a row in the pssm, so its length is 4"""
    r = random.random()
    i = 0
    while r > p[i]:
        r -= p[i]
        i += 1
    return i

def pssm_to_fasta(pssm, fastafilename, n=1000):
    """Generate a fastafile that approximates the base frequences of the pssm"""
    base = ('A', 'C', 'G', 'T')
    fastafile = open(fastafilename, 'w')
    for seqnum in xrange(n):
        seqstr = "".join(base[sample(p)] for p in pssm)
        print >> fastafile, ">%d" % seqnum
        print >> fastafile, seqstr
    fastafile.close()

def reverse_pssm(pssm):
    """Return a reversed pssm"""
    pssm_rev = [t[::-1] for t in pssm[::-1]]
    return pssm_rev

#NOTE: this fn is obsolete!
def make_logo(motif, width, height, img_dir):
    """Generate motif logo, place it in the img_dir and return its file name.
    """
    # create directory for holding motif logos if directory does not
    # already exist
    if not os.path.exists(img_dir):
        os.makedirs(img_dir)
    # generate temporary FASTA file which approximates the motif PSSM; input for seqlogo
    fasta_file_name = os.path.join(img_dir, 'temp.fa')
    pssm_to_fasta(motif.pssm, fasta_file_name)
    logo_file_name = os.path.join(img_dir, motif.id)
    # run the seqlogo command
    command = "%s -f %s -o %s -w %d -h %d -F PNG -c -n -Y" % \
              (os.path.join(settings.DEPLOY_DIR, 'weblogo', 'seqlogo'),
               fasta_file_name, logo_file_name, width, height)
    os.system(command)
    # delete the temporary FASTA file
    os.remove(fasta_file_name)
                    
    # return the filename of the motif logo, strip out the _OUTPUT_DIR
    logo_file_name += '.png'
    return logo_file_name

def save_to_html(output_dir, motifList, motifDists):
    """Saves the list of motifs as an html file named 'mdseqpos_out.html'.
    The motif logos associated with mdseqpos_out.html will be stored
    under output_dir/img i.e. results/img.
    """
    # OBSOLETE: moving to goLogo to generate the motif logos
    # #make results and results/img
    # if not os.path.exists(output_dir):
    #     os.makedirs(output_dir)
    # if not os.path.exists(os.path.join(output_dir, 'img')):
    #     os.makedirs(os.path.join(output_dir, 'img'))

    # #make the motif logos, put them in results/img
    # for (i, m) in enumerate(motifList):
    #     #NOTE: logoImg isn't part of the class definition, but we add it on 
    #     #to the instanc in hopes that the Motif.to_json will know how to 
    #     #serialize it.
    #     tmp = make_logo(motifList[i], 192 / 20.0, 120 / 20.0,
    #                     os.path.join(output_dir, 'img'))
    #     #need to strip out the output_dir from the path, e.g. 
    #     #results/img/denovo0.png --> img/denovo.png
    #     motifList[i].logoImg = tmp.replace(output_dir+"/",'')

    #create the output file--
    mdseqpos_out_file = open(os.path.join(output_dir,"mdseqpos_out.html"), 'w')
    motifList_json = motifList.to_json() 
    motifDists_json = json.dumps(motifDists)
    t = loader.get_template('mdseqpos_out.html')
    c = Context({'motifList': motifList_json, 'motifDists':motifDists_json, 
                 'now':datetime.now()}) 
    mdseqpos_out_file.write(t.render(c))

    #copy over supporting files--i.e. everything in the template/static dir
    rootdir = os.path.join(settings.DEPLOY_DIR, 'template', 'static','')
    for (path, dirs, files) in os.walk(rootdir):
        for f in files:
            subdir = path.replace(rootdir, '')
            if subdir:
                #check to make sure that the subdir exists
                if not os.path.exists(os.path.join(output_dir, subdir)):
                    os.makedirs(os.path.join(output_dir, subdir))

            shutil.copyfile(os.path.join(path, f), 
                            os.path.join(output_dir, subdir, f))

    #END save_to_html

def save_to_html_plain(output_dir, motifList, distCutoff):   
    # > dir(motifList[0])
    # ['_ATTRIBUTES', '__doc__', '__init__', '__module__', '__str__', '_results_fields', 
    # '_validpssm', 'antisense', 'dbd', 'entrezs', 'equals', 'factors', 'from_dict', 
    # 'from_flat_file', 'from_xml', 'fullname', 'getText', 'getcutoff()', 'getmeanposition()', 
    # 'getnumhits()', 'getpvalue()', 'getzscore()', 'id', 'numseqs', 'pmid', 'pssm', 
    # 'pssm_from_xml', 'pssm_to_xml', 'refseqs', 'seqpos', 'seqpos_results', 'seqpos_stat', 
    # 'setpssm', 'source', 'sourcefile', 'species', 'status', 'symbols', 'synonyms', 
    # 'to_json', 'to_xml']
    #
    # > motifList[0].seqpos_results.keys()
    # ['plot', 'zscore', 'seq', 'cutoff', 'pssm', 'meanposition', 'numhits', 'pvalue']

    # fix pssm, and change pvalue to -10log(pval, e)
    mlist = []
    for motif in motifList:
        for row in motif.seqpos_results['pssm']:
            rsum = sum([float(t) for t in row])
            for t in range(len(row)):
                row[t] = round(float(row[t]) / rsum, 3)
                if row[t] < 0.001:
                    row[t] = 0.001
            imax = row.index(max(row))
            row[imax] -= sum(row) - 1
        motif.seqpos_results['pvalue'] = round(-10*math.log(motif.getpvalue(), math.e), 3)
        motif.seqpos_results['pssm_rev'] = reverse_pssm(motif.seqpos_results['pssm'])
        mlist.append(motif)
    
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    
    #create seqLogo by run Rscript.
    seqLogoFolder = os.path.join(output_dir, 'seqLogo')
    if not os.path.exists(seqLogoFolder):
        os.mkdir(seqLogoFolder)
    
    #draw seqLogo to png
    rfile = os.path.join(seqLogoFolder, 'draw_seqLogo.r')
    rscript = open(rfile, 'w')
    rscript.write('setwd("%s")\n' %os.path.abspath(seqLogoFolder))
    rscript.write('library("seqLogo")\n')
    for motif in mlist:
        # create png for both seqlogo and reversed seqlogo.
        for mid, pssm in ((motif.id, motif.seqpos_results['pssm']), (motif.id + '_rev', motif.seqpos_results['pssm_rev'])):
            t1 = ['c(%s)' % ','.join([str(t2) for t2 in t]) for t in pssm]
            t2 = 'data<-cbind(%s)\n' % ','.join(t1)
            rscript.write('png("%s.png", width=660, height=300)\n' %mid)
            rscript.write(t2)
            rscript.write('seqLogo(as.matrix(data))\n')
            rscript.write('dev.off()\n\n')
    rscript.close()
    cmd = 'Rscript %s' %rfile
    os.system(cmd)
        
    #create each motif's single page
    motifFolder = os.path.join(output_dir, 'motif')
    if not os.path.exists(motifFolder):
        os.mkdir(motifFolder)
    pssmFolder = os.path.join(output_dir, 'pssm')
    if not os.path.exists(pssmFolder):
        os.mkdir(pssmFolder)
    for motif in mlist:
        render_to_file('single_motif.html', {'motif': motif}, os.path.join(motifFolder, motif.id + '.html'))
        render_to_file('single_pssm.html', {'motif': motif}, os.path.join(pssmFolder, motif.id + '.html'))
    
    #create hit seq pages
    #hitseqFolder = os.path.join(output_dir, 'hitseq')
    #if not os.path.exists(hitseqFolder):
    #    os.mkdir(hitseqFolder)
    #for arg in args:
    #    outf = open(os.path.join(hitseqFolder, arg['id']+'.txt'), 'w')
    #    outf.writelines([t+'\n' for t in arg['hitseq']])
    #    outf.close()
    
    #collapse motifs
    flat_clusters = pmc.motif_hcluster2(mlist, distCutoff)
    flat_clusters.sort(key = lambda x: min([t.getzscore() for t in x]))
    
    m_collapse = []
    for c in flat_clusters:
        c.sort(key=lambda x:x.getzscore())
        m_collapse.append(c.pop(0))
        
        m_collapse[-1].similar_motifs = [m_collapse[-1]] #also put self into similarity_motifs list.
        for m in c:
            similarity_score = pmc.similarity(m.seqpos_results['pssm'], m_collapse[-1].seqpos_results['pssm'])[0]
            m.similarity_score = round(similarity_score, 3)
            m_collapse[-1].similar_motifs.append(m)

    for i in range(len(m_collapse)):
        motif = m_collapse[i]
        motif.collapse_num = len(motif.similar_motifs)
        motif.class_id = i + 1
    
    #create table page and home page.
    render_to_file('table.html', {'motifs': m_collapse}, os.path.join(output_dir, 'table.html'))
    render_to_file('mdseqpos_index.html', {}, os.path.join(output_dir, 'mdseqpos_index.html'))

def render_to_file(template_html, render_dict, filen):
    template_d = os.path.join(settings.DEPLOY_DIR, 'template')
    template_f = open(os.path.join(template_d, template_html))
    template = Template(template_f.read())
    outf = open(filen, 'w')
    outf.write(template.render(render_dict))
    outf.close()
    template_f.close()

def calc_motif_dist(motifList):
    """Given a list of motifs, returns a dictionary of the distances
    for each motif pair, e.g. {Motif1:Motif2:(dist, offset, sense/antisense)}
    """
    ret = {}
    for m1 in motifList:
        for m2 in motifList:
            if m1.id != m2.id:
                #check if the score is already calculated
                if ("%s:%s" % (m1.id,m2.id)) not in ret and \
                        ("%s:%s" % (m2.id,m1.id)) not in ret:
                    #must send in pssms as numpy arrays
                    dist = bmc.BLiC_score(numpy.array(m1.pssm, float),
                                          numpy.array(m2.pssm, float))
                    ret["%s:%s" % (m1.id,m2.id)] = dist
    return ret
    
def calc_motif_dist_pcc(motifList):
    """Given a list of motifs, returns a dictionary of the distances
    for each motif pair, e.g. {Motif1:Motif2:(dist, offset, sense/antisense)}
    """
    ret = {}
    for m1 in motifList:
        for m2 in motifList:
            if m1.id != m2.id:
                #check if the score is already calculated
                if ("%s:%s" % (m1.id,m2.id)) not in ret and \
                        ("%s:%s" % (m2.id,m1.id)) not in ret:
                    #must send in pssms as numpy arrays
                    dist = pmc.similarity(m1.pssm, m2.pssm)
                    #dist = bmc.BLiC_score(numpy.array(m1.pssm, float),
                    #                      numpy.array(m2.pssm, float))
                    ret["%s:%s" % (m1.id,m2.id)] = dist
    return ret
            
def main():
    #ALWAYS PRINT OUT VERSION INFO: 
    print mdseqpos.__version__
    print 'Library path:', mdseqpos.__file__
    print 
    
    USAGE = """USAGE: MDSeqPos.py [options] BEDFILE GENOME
    
    BEDFILE - regions file
    GENOME  - assembly which the regions pertain to, e.g. 'hg18', 'mm9', etc.
              as defined in BUILD_DICT in lib/settings.py"""
              
    parser = optparse.OptionParser(usage=USAGE)
    parser.add_option('-g', '--genome-dir', dest="genome_dir", default=None,
                      help="Path to the genome assembly dir")
    parser.add_option('-d', '--denovo', default=False, action="store_true",
                      help="flag to run denovo motif search (default: False)")
    parser.add_option('-m', '--known-motifs', default=None,
                      help="comma separated list of known motifs dbs to use in the motif search, e.g. -m pbm.xml,transfac.xml")
    parser.add_option('-n', '--new-motifs', default='denovo.xml',
                      help="name of the output XML file which stores new motifs found during adenovo search, e.g. -n foo.xml (default: denovo.xml)")
    parser.add_option('-p', '--pval', default=0.001,
                      help="pvalue cutoff for motif significance, (default: 0.001)")
    parser.add_option('-s', '--species-list', default=None, 
                      help="name of species to filter the results with--if multuple species, comma-separate them, e.g. hs,mm,dm")
    parser.add_option('-w', '--width', default=600,
                      help="width of the region to be scanned for motifs; depends on resoution of assay, (default: 600 basepairs)")
    parser.add_option('-v', '--verbose', default=False, action="store_true",
                      help="flag to print debugging info (default: False)")
    parser.add_option('--hcluster', type = 'float', default=0.8,
                      help="The similarity cutoff for hierarchical clustering of the result, (default: 0.8, The higher, the more groups, 0 ~ 1)")
    parser.add_option('-c', '--cluster', default=False, action="store_true",
                      help="This option only for know-motifs cistrome.xml, If you want to use pre-clustered database to accelerate seqpos. default (not set): False")
    parser.add_option('--maxmotif', default=0,
                      help="maximum number of motifs to report, (default: 0, i.e. no max)")
    parser.add_option('-O', '--output-directory', default="results", 
                      help="output directory name (default: results)")
    
    #parse the command line options
    (opts, args) = parser.parse_args(sys.argv)
    if len(args) < 3: 
        parser.print_help()
        sys.exit(-1)
    bedfile_name = args[1]
    genome = args[2]
    _DEBUG = opts.verbose
    output_dir = opts.output_directory
    print opts.genome_dir

    #READ in the regions that are specified in the BED file
    print "read regions start time: %s" % time.ctime()
    #HERE we should rely on a standard package to read in bed files; stub it
    chip_regions = ChipRegions(bedfile_name, genome, genome_dir=opts.genome_dir)
    print "read regions end time: %s" % time.ctime()

    #LOAD the motifs (both known and denovo)
    known_motifs, new_motifs = None, None
    if opts.known_motifs:
        motif_dbs = [x.strip() for x in opts.known_motifs.split(',')]
        known_motifs = read_known_motifs(motif_dbs, _DEBUG)

    if opts.denovo:
        print "starting denovo search...(time: %s)" % time.ctime()
        new_motifs = chip_regions.mdmodule(width=int(opts.width))
        print "completed denovo search...(time: %s)" % time.ctime()
        new_motifs.save_to_xml(os.path.join(output_dir, opts.new_motifs))
        
    #Run seqpos stats on all_motifs
    print "starting seqpos stats...(time: %s)" % time.ctime()
    if new_motifs:
        for m in new_motifs:
            m.seqpos(chip_regions, width=int(opts.width))
    if opts.cluster and opts.known_motifs == 'cistrome.xml': #only for cistrome.xml to use cistrome.cluster
        known_motifs = seqpos_cluster_known_motifs(known_motifs, chip_regions, int(opts.width), float(opts.pval))
    elif known_motifs:
        for m in known_motifs:
            m.seqpos(chip_regions, width=int(opts.width))
    print "completed seqpos stats...(time: %s)" % time.ctime()

    #Combine both known and new motifs
    all_motifs = None
    if known_motifs and new_motifs:
        all_motifs = MotifList(known_motifs + new_motifs)
    elif known_motifs:
        all_motifs = known_motifs
    else:
        all_motifs = new_motifs

    #CULL the results to see only the relevant results, and output
    sig_motifs = all_motifs.cull(float(opts.pval), int(opts.maxmotif))
    
    #filter by species?
    if opts.species_list:
        species_list = opts.species_list.split(',')
        sig_motifs = sig_motifs.filterBySpecies(species_list)

    #dists = calc_motif_dist_pcc(sig_motifs)
    #save_to_html(output_dir, sig_motifs, dists)
    save_to_html_plain(output_dir, sig_motifs, opts.hcluster)

    json_list = [t.to_json() for t in sig_motifs]
    jsonf = open(os.path.join(output_dir, 'motif_list.json'),'w')
    for js in json_list:
        jsonf.write(js +'\n')
    jsonf.close()

if __name__ == '__main__':
    main()

#! /usr/bin/env python

"""A command-line script for scanning ChIP regions for one motif."""

import sys, os

from optparse import OptionParser

import mdseqpos
from mdseqpos.chipregions import ChipRegions
from mdseqpos.motif import Motif, MotifList

import mdseqpos.settings as settings
#from mdseqpos.motif import Motif

USAGE = """Usage: MotifScan.py [options] BEDFILE GENOME MOTIFID

Arguments:             
  BEDFILE              name of an existing BED file containing ChIP regions,
                       to be scanned for the given motif
  GENOME               abbreviation identifying the genome build that is
                       used in the BED file (must be 'hg18', 'mm3', or ... )
  MOTIFID              string ID of the motif to be used for scanning the
                       ChIP regions"""

NEW_MOTIFS_FILE = 'NULL'
NEW_MOTIFS_DIR = None #'results/'
KNOWN_MOTIFS_FILE = 'NULL' #'known.xml'
KNOWN_MOTIFS_DIR = os.path.join(settings.DEPLOY_DIR, "database")
BED_FILE = 'motif_scan_out.bed'
BED_DIR = '.'
FASTA_FILE = 'motif_scan_out.fsa'
FASTA_DIR = '.'

if __name__ == '__main__':
    #ALWAYS PRINT OUT VERSION INFO: 
    print mdseqpos.__version__

    # parse command line arguments
    # 
    # Required Arguments:
    # 1. chip_regions_path -- input BED file containing user's ChIP regions
    # 2. genome -- genome identifier (e.g. 'hg18')
    # 3. motif_id -- ID of the motif
    # 
    # Optional Arguments:
    # 1. known_motifs_file, known_motifs_dir -- input XML file where known motifs are stored
    # 2. new_motifs_file, new_motifs_dir -- input XML file where the user's motif is stored
    # 3. bed_file, bed_dir -- output BED file
    # 4. fasta_file, fasta_dir -- output FASTA file
    # 
    parser = OptionParser(usage=USAGE)
    parser.add_option('-g', '--genome-dir', dest="genome_dir", default=None,
                      help="Path to the genome assembly dir")
    parser.add_option('-m', '--known-motifs-file', default=KNOWN_MOTIFS_FILE,
                      help="name of input XML file containing known motifs")
    parser.add_option('-M', '--known-motifs-dir', default=KNOWN_MOTIFS_DIR,
                      help="directory of input XML file containing known motifs")
    parser.add_option('-n', '--new-motifs-file', default=NEW_MOTIFS_FILE,
                      help="name of input XML file containing new motifs discovered by \
                      de novo scan of ChIP regions")
    parser.add_option('-N', '--new-motifs-dir', default=NEW_MOTIFS_DIR,
                      help="directory of input XML file containing new motifs")
    parser.add_option('-b', '--bed-file', default=BED_FILE,
                      help="output BED file containing the genomic regions that are MotifScan\
                      hits")
    parser.add_option('-B', '--bed-dir', default=BED_DIR,
                      help="directory of the output BED file")
    parser.add_option('-f', '--fasta-file', default=FASTA_FILE,
                      help="output FASTA file containing the genomic sequences that are\
                      MotifScan hits")
    parser.add_option('-F', '--fasta-dir', default=FASTA_DIR,
                      help="directory of the output FASTA file")
    parser.add_option('-p', '--pssm-file', default=None,
                      help="input a pssm file that represents the motif")
    (opts, args) = parser.parse_args(sys.argv)
    
    # required arguments
    chip_regions_path = args[1]
    genome = args[2]
    motif_id = args[3]
    
    # set known motifs file path
    if opts.known_motifs_file == 'NULL':
        known_motifs_path = None
    else:
        known_motifs_path = os.path.join(opts.known_motifs_dir, opts.known_motifs_file)
        
    # set new motifs file path
    if opts.new_motifs_file == 'NULL':
        new_motifs_path = None
    else:
        if opts.new_motifs_dir is None:
            new_motifs_path = opts.new_motifs_file
        else:
            new_motifs_path = os.path.join(opts.new_motifs_dir, opts.new_motifs_file)
            
    # set BED file path
    bed_path = os.path.join(opts.bed_dir, opts.bed_file)
        
    # set FASTA file path
    fasta_path = os.path.join(opts.fasta_dir, opts.fasta_file)
    
    # retrieve known motifs
    motifs = MotifList()
    if known_motifs_path is not None:
        motifs.from_xml_file(known_motifs_path)
    
    # retrieve new motifs if new motifs file is specified by user
    if new_motifs_path is not None:
        new_motifs = MotifList()
        new_motifs.from_xml_file(new_motifs_path)
        #motifs.append(new_motifs)
        motifs += new_motifs
    
    # scan all motifs for the desired motif ID
    if opts.pssm_file is None:
        for motif in motifs:
            if motif.id == motif_id:
                desired_motif = motif
                break
    else: #use a pssm as the flatfile
        desired_motif = Motif.from_flat_file(opts.pssm_file)
    
    # retrieve ChIP regions
    chip_regions = ChipRegions(chip_regions_path, genome, opts.genome_dir)
    
    # scan ChIP regions for motif
    hits = chip_regions.motifscan(desired_motif)
    
    # save hits as BED file and/or FASTA file
    if bed_path is not None:
        bed_file = open(bed_path, 'w')
        bed_file.write(hits.to_bed())
        bed_file.close()
    if fasta_path is not None:
        fasta_file = open(fasta_path, 'w')
        fasta_file.write(hits.to_fasta())
        fasta_file.close()

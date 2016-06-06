"""
SYNOPSIS: A module for ChIP regions.

Adapted from Cliff Meyer's cistrome.py module.
"""

import os
import sys
from operator import itemgetter
from copy import deepcopy

#NOTE: when seqscan is modified to take pssms that aren't numpy arrays,
#then i can remove this!
import numpy

try:
    from _MDmod import MDmod
    import _seq
except:
    from mdseqpos._MDmod import MDmod
    import mdseqpos._seq as _seq

import motif
import settings
import count

#NOTE: this should maybe be moved to lenlib.core.regions - for maximum reuse
class ChipRegions:
    """A class for ChIP regions.

    Adapted from Cliff Meyer's ChIP_region class."""
    def __init__(self, bedfilename=None, genome=None, genome_dir=None):
        """Initialize ChIP region.

        Adapted from Cliff Meyer's ChIP_region.__init__() method."""
        self.chrom = []
        self.chromStart = []
        self.chromEnd = []
        self.sequence = []
        #note: strand & hits info is only used/loaded when motifscan is called
        self.strand = []
        self.hitscore = []
        self.genome = None
        self.genome_dir = None
        if bedfilename is not None:
            self.read(bedfilename)
        if genome is not None:
            self.genome = genome
        if genome_dir is not None:
            self.genome_dir = genome_dir
        print self.genome_dir
        self.preprocessed_regions = False

    def preprocess(self, width, margin):
        """Preprocess the chip regions, you only need to do this once (this
        is also destructive) so that is why we set the preprocessed_regions
        flag to True.
        """
        self.preprocessed_regions = True;
        self.unique()
        self.trim(width + margin)

    def read(self, bedfilename):
        """Read ChIP regions from a BED file.

        This function only reads the first three columns of BED files:
          (1) Chromosome Name (e.g., 'chr7')
          (2) Start Position (e.g., '41567')
          (3) End Position (e.g., '41699')

        Adapted from Cliff Meyer's ChIP_region.ReadBed() method."""
        bedfile = open(bedfilename, 'r')
        for line in bedfile:
            if line.isspace():
                continue
            if line[:5] == 'track':
                continue
            row = line.split()
            self.chrom.append(row[0])
            self.chromStart.append(int(row[1]))
            self.chromEnd.append(int(row[2]))

    def len_read_sequence(self, repeatMasked=False):
        """Read in region sequences, len style!"""
        build_dir = self.genome_dir if self.genome_dir else os.path.join(settings.ASSEMBLY_DIR,
                                 settings.BUILD_DICT[self.genome])
        #SORT the reads so that we're more efficient about our File I/O
        tuples = zip(range(len(self.chrom)), self.chrom, self.start, self.end)
        #sort by chr and then start
        sorted_tuples = sorted(tuples, key=itemgetter(1,2))

        unsorted_seq = []
        last_file, fpchr = None, None
        for (order, chr, start, end) in sorted_tuples:
            if repeatMasked:
                chr_file = os.path.join(build_dir, chr+'.fa')
            else:
                pass

            #only open the files when we move to a different chr
            if chr_file != last_file:
                try:
                    if fpchr: fpchr.close()
                    fpchr = open(chr_file, 'r')
                    last_file = chr_file
                except:
                    print "Could not open file:", chr_file
                    unsorted_seq.append((order, None))
                    continue
            else: #continue using same fpchr
                #SCREW it, movie time!
                pass

    def ext_read_sequences(self, repeatMAsked=False):
        """Uses Jim Kent's twoBitToFa to read in the sequences associated
        with the regions.  twoBitToFa is included in the mdseqpos package,
        under 'tools/twoBitToFa'.

        The location of the two bit files are defined in settings, i.e.
        mm8's twoBit is in settings.TWOBIT_DIR/mm8.2bit
        """
        twoBit_path = os.path.join(settings.DEPLOY_DIR, 'tools', 'twoBitToFa')
        assembly_path = os.path.join(settings.TWOBIT_DIR, self.genome+".2bit")
        #twoBitToFa is fussy about bed file formats, so we have to
        #make a temporary bed to ensure its compliance
        f = open('temp.bed','w')
        f.write(self.to_4col_bed())
        f.close()

        cmd = "%s %s delete_this.fa -bed=%s" % (twoBit_path, assembly_path,
                                                'temp.bed')
        #print "CMD is %s" % cmd
        (status, log) = commands.getstatusoutput(cmd)
        #print "Status is %s, Log is:\n%s\n" % (status, log)

        #read in the sequences from the fasta file:
        f = open('delete_this.fa')
        try:
            self.sequence = []
            tmp = ""
            for line in f:
                #print "#%s#" % line.strip()
                if line.strip() == ">IGNORE_THIS": #new sequence
                    self.sequence.append(tmp.upper())
                    tmp = ""
                else:
                    tmp += line.strip()
            if tmp: self.sequence.append(tmp.upper())
        finally:
            f.close()

        #NOTE: the first sequence is always empty--so we drop it
        if self.sequence and len(self.sequence) > 0:
            self.sequence = self.sequence[1:]

        #delete temporary files
        if os.path.exists('delete_this.fa'): os.remove('delete_this.fa')
        if os.path.exists('temp.bed'): os.remove('temp.bed')


    def read_sequence(self, repeatMasked=False):
        """Retrieve sequences for ChIP regions.

        Adapted from Cliff Meyer's ChIP_region.getSequence() method."""
        #BUILD_DICT has been moved to settings.py
        build_dir = self.genome_dir if self.genome_dir else os.path.join(settings.ASSEMBLY_DIR,
                                 settings.BUILD_DICT[ self.genome ])
        #print build_dir
        self.sequence = []
        for i,x in enumerate(self.chrom):
            chr   = x
            left  = self.chromStart[i]
            right = self.chromEnd[i]
            if repeatMasked == True:
                chr_file = os.path.join(build_dir, chr+'.fa')
                try:
                    fpchr = open(chr_file, 'r')
                    #print chr_file + "aa"
                except:
                    continue
            else:
                subdir = os.listdir(build_dir)
                if "raw" in subdir:
                    chr_file = os.path.join(build_dir, 'raw', chr + ".fa")
                    #print chr_file + "bb"
                else:
                    chr_file = os.path.join(build_dir, 'rawgenome', chr+".fa")
                    #print chr_file + "cc"
                try:
                    fpchr = open(chr_file, 'r')
                except:
                    print "Could not open file: %s--ignoring:%s,%s,%s" % (chr_file, chr, left, right)
                    #NOTE: Should be ignored--not set to None
                    #self.sequence.append(None)
                    continue
            fpchr.readline()
            here = fpchr.tell()
            line = fpchr.readline()
            line = line.rstrip()
            line_len = len(line)
            left  -= 1
            right -= 1
            start = (left  / line_len) * (line_len + 1) + (left  % line_len)
            stop  = (right / line_len) * (line_len + 1) + (right % line_len)

            try:
                fpchr.seek(start + here - 25)
            except:
                print 'WARNING: requesting invalid chr location before chromoso\
me start', start, here, left, right, line_len, x
                fpchr.close()
                continue
            try:
                fpchr.seek(stop + here + 25)
            except:
                print 'WARNING: requesting invalid chr location over chromosome\
 end', stop, here, left, right, line_len, x
                fpchr.close()
                continue

            fpchr.seek(start + here)

            seq = fpchr.read(stop - start)
            seq = seq.replace('\n', '')
            self.sequence.append(seq.upper())
            fpchr.close()

    def to_fasta(self):
        """Return ChIP regions as a string in FASTA format.

        Adapted from Cliff Meyer's print_fasta() function."""
        fasta_out = ''
        for i,chrom in enumerate(self.chrom):
            fasta_out += ">"+str(self.chromStart[i])+"\n"+self.sequence[i]+"\n"
        return fasta_out

    def to_bed(self):
        """Return ChIP regions as a string in BED format.

        Adapted from Cliff Meyer's __str__() method."""
        bed_out = ''
        for i,chrom in enumerate(self.chrom):
            bed_out += chrom+"\t"+str(self.chromStart[i])+"\t"+str(self.chromEnd[i])+"\t"+self.sequence[i]+"\t"+str(self.hitscore[i])+"\t"+self.strand[i]+"\n"
        return bed_out

    def trim(self, length=600):
        """Trim all ChIP regions to a given length.

        The same number of basepairs (give or take 1 basepair) will be
          trimmed from either end of each ChIP region.

        Adapted from Cliff Meyer's ChIP_region.TrimBed() method."""
        half_length = length / 2.0
        for i in xrange(len(self.chrom)):
            midpoint = (self.chromStart[i] + self.chromEnd[i]) / 2.0
            self.chromStart[i] = int(midpoint - half_length)
            self.chromEnd[i] = int(midpoint + half_length)

    def unique(self):
        """Remove duplicate ChIP regions."""
        pass

    def copy(self):
        """Return a deep copy of self."""
        return deepcopy(self)

    def mdmodule(self, motif_widths=(7, 10, 13, 16, 19), width=600):
        """Run a de novo motif scan of ChIP regions.

        De novo motifs will be returned as a MotifList object.

        Adapted from Cliff Meyer's ChIP_region.MDscan() method.
        Directly calls Xiaole Shirley Liu's _MDmod program."""
        motifs = motif.MotifList()
        # prepare input sequences
        input = self.copy()
        input.unique()
        input.trim(int(width/2))
        input.read_sequence(True)

        #OBSOLETE: DELETE!!
        # prepare background sequences
        background = self.copy()
        background.unique()
        background.trim(int(width))
        background.read_sequence(True)

        # scan for motifs
        for motif_width in motif_widths:
            raw_pssms = MDmod(i=input.sequence, b=background.sequence, w=motif_width, t=50, s=50, n=100, r=10)
            for raw_pssm in raw_pssms:
                motif_id = 'denovo%d' % len(motifs)
                tmp = motif.Motif()
                tmp.id = motif_id
                tmp.pssm = raw_pssm
                motifs.append(tmp)
        return motifs

    #SHOULD make this a static method??--need to pass in chipregions.??
    def motifscan(self, motif):
        """Scan sequences with a motif.

        Motif must be provided as a Motif object.
        Hits will be returned as a ChipRegions object.

        Directly calls [Somebody]'s _seq program.
        """
        SENSE = 0
        self.read_sequence(True)

        # LEN: BINOCH UPGRADE
        bgseqprob_mat = count.count(self.sequence)
        markov_order = 2
        prob_option = _seq.CUTOFF_OPTION

        #GRR: THIS IS THE ONLY numpy dependency, and it is b/c seqscan
        #expects self.pssm to be a numpy array!
        pssm = numpy.array(motif.pssm, float)
        s_idx, start, end, orient, score = \
               _seq.seqscan(self.sequence, pssm, bgseqprob_mat,
                            markov_order, prob_option)
        hits = ChipRegions(genome=self.genome,genome_dir=self.genome_dir)
        hits.genome = self.genome

        for i,idx in enumerate(s_idx):
            #sorry for the i vs idx confusion, but they're really different!
            hits.chrom.append(self.chrom[idx])
            hits.chromStart.append(int(self.chromStart[idx] + int(start[i])))
            hits.chromEnd.append(int(self.chromStart[idx] + int(end[i])))
            if int(orient[i]) == SENSE:
                hits.strand.append("+")
            else:
                hits.strand.append("-")
            #hits.hitscore.append(score[idx])--not sure if i should use i/idx
            hits.hitscore.append(score[i])

        # LEN: BINOCH UPGRADE END

        hits.read_sequence(True)

        return hits

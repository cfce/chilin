#!/usr/bin/env python

"""
Time-stamp: <2014-11-14 02:34:09 qqin>
-d: input one merged Phastcon conservation bw file
    e.g. rsync -avz --progress rsync://hgdownload.soe.ucsc.edu/goldenPath/hg38/phastCons7way/hg38.phastCons7way.bw 
"""
# ------------------------------------
# python modules
# ------------------------------------

import os
import sys
import re
import logging
import subprocess
import math
from optparse import OptionParser

try:
    from bx.bbi.bigwig_file import BigWigFile
except:
    sys.stderr.write("Need bx-python!")
    sys.exit()

# ------------------------------------
# constants
# ------------------------------------
logging.basicConfig(level=20,
                    format='%(levelname)-5s @ %(asctime)s: %(message)s ',
                    datefmt='%a, %d %b %Y %H:%M:%S',
                    stream=sys.stderr,
                    filemode="w"
                    )
#bigWigSummary = 'bigWigSummary'

# ------------------------------------
# Misc functions
# ------------------------------------

error  = logging.critical		# function alias
warn   = logging.warning
debug  = logging.debug
info   = logging.info

class PeakIO:
    """IO for peak region information.

    This class can hold peak information from MAT/MA2C/MACS, and
    provide some extra functions:

    1. filtering functions, filter_pvalue/score/fdr/fold are functions
    to filter the peaks according to pvalue/score/fdr/fold range given
    by the user.

    2. overlapping functions

    """
    def __init__ (self, comment=""):
        """Initialization function.

        comment: you can add any comments to the peakIO object like
        whether or not it is from a ChIP-chip or a ChIP-seq
        experiments.
        """
        self.peaks = {}
        self.comment = comment

    def dup (self):
        """return a duplicate peakI.
        """
        r = PeakIO(comment=self.comment)
        peaks = self.peaks
        new_peaks = {}
        chrs = peaks.keys()
        chrs.sort()
        for chrom in chrs:
            new_peaks[chrom]=[]
            new_peaks[chrom].extend(peaks[chrom])
        r.peaks = new_peaks
        return r


    def add (self, chromosome, start, end, summit=None,
             score=None, total_p=None,
             pvalue=None, fold_enrichment=None, fdr=None):
        """Use this function to add items to PeakIO object.

        items: (peak start,peak end, peak length, peak summit, peak
        score, number of tags/probes in peak region, peak pvalue, peak
        fold_enrichment, fdr) <-- tuple type

        Parameters:
        1. chromosome
        2. start
        3. end
        4. summit: the highest position for the peak region
        5. score:  the score for peak region
        6. total_p: total points in peak region. For ChIP-seq, it's
        how many tags in the region; for ChIP-chip, it's the number
        of probes.
        7. pvalue: -10*log(10,p-value) for peak region
        8. fold_enrichment: fold enrichment for the region
        9. fdr: False Discovery Rate for the region
        """
        if not self.peaks.has_key(chromosome):
            self.peaks[chromosome]=[]
        self.peaks[chromosome].append((start,end,end-start,summit,
                                       score,total_p,
                                       pvalue,fold_enrichment,fdr))

    def filter_pvalue (self, pvalue_cut_low, pvalue_cut_up=None ):
        """Filter peaks in a given pvalue range.

        Note, pvalue is actually -10*log(10,pvalue)

        If pvalue_cut_low and pvalue_cut_up is assigned, the peaks with pvalue in [pvalue_cut_low,pvalue_cut_up).
        """
        peaks = self.peaks
        new_peaks = {}
        chrs = peaks.keys()
        chrs.sort()
        if pvalue_cut_up:
            for chrom in chrs:
                new_peaks[chrom]=[p for p in peaks[chrom] if p[6] >= pvalue_cut_low and p[6]<pvalue_cut_up]
                if not new_peaks[chrom]: del new_peaks[chrom]
        else:
            for chrom in chrs:
                new_peaks[chrom]=[p for p in peaks[chrom] if p[6] >= pvalue_cut_low]
                if not new_peaks[chrom]: del new_peaks[chrom]
        self.peaks = new_peaks

    def filter_score (self, score_low, score_up=None ):
        """Filter peaks in a given score range.

        If score_low and score_up is assigned, the peaks with score in [score_low,score_up).

        """
        peaks = self.peaks
        new_peaks = {}
        chrs = peaks.keys()
        chrs.sort()
        if score_up:
            for chrom in chrs:
                new_peaks[chrom]=[p for p in peaks[chrom] if p[4] >= score_low and p[4]<score_up]
                if not new_peaks[chrom]: del new_peaks[chrom]
        else:
            for chrom in chrs:
                new_peaks[chrom]=[p for p in peaks[chrom] if p[4] >= score_low]
                if not new_peaks[chrom]: del new_peaks[chrom]
        self.peaks = new_peaks

    def filter_fold (self, fold_low, fold_up=None ):
        """Filter peaks in a given fold enrichment range.

        If fold_low and fold_up is assigned, the peaks with fold in [fold_low,fold_up)

        """
        peaks = self.peaks
        new_peaks = {}
        chrs = peaks.keys()
        chrs.sort()
        if fold_up:
            for chrom in chrs:
                new_peaks[chrom]=[p for p in peaks[chrom] if p[7] >= fold_low and p[7]<fold_up]
                if not new_peaks[chrom]: del new_peaks[chrom]
        else:
            for chrom in chrs:
                new_peaks[chrom]=[p for p in peaks[chrom] if p[7] >= fold_low]
                if not new_peaks[chrom]: del new_peaks[chrom]
        self.peaks = new_peaks

    def filter_fdr (self, fdr_up, fdr_low=None ):
        """Filter peaks in a given FDR range.

        If fdr_low and fdr_up is assigned, the peaks with fold in (fdr_low,fdr_up]. Otherwise, return the peaks with FDR lower or equal to fdr_up.

        """
        peaks = self.peaks
        new_peaks = {}
        chrs = peaks.keys()
        chrs.sort()
        if fdr_low:
            for chrom in chrs:
                new_peaks[chrom]=[p for p in peaks[chrom] if p[8] > fdr_low and p[8]<=fdr_up]
                if not new_peaks[chrom]: del new_peaks[chrom]
        else:
            for chrom in chrs:
                new_peaks[chrom]=[p for p in peaks[chrom] if p[8] <= fdr_up]
                if not new_peaks[chrom]: del new_peaks[chrom]
        self.peaks = new_peaks

    def sort (self):
        peaks = self.peaks
        chrs = peaks.keys()
        chrs.sort()
        for chrom in chrs:
            peaks[chrom].sort(lambda x,y: cmp(x[0],y[0]))


    def get_chr_names (self):
        """Return all the chromosome names stored.

        """
        l = set(self.peaks.keys())
        return l


def cmp(a, b):
    return (a > b) - (a < b)


def parse_BED (fhd):
    """Parse a tab-delimited bed file

    Return a PeakIO object containing peak regions.
    """
    import subprocess
    peaks = PeakIO()
    n=0
    for thisline in fhd:
        n+=1
        if n>5000: ## use top 5000 peak
            break
        thisline = thisline.rstrip()
        if not thisline: continue #return ("blank",None,None)
        if thisline.startswith("#"): continue #return ("comment line",None,None) # comment line is skipped
        if thisline.startswith("track"): continue
        if thisline.startswith("browser"): continue
        thisfields = thisline.split()
        startpos = max(0,int(thisfields[1]))

        peaks.add(thisfields[0],startpos,int(thisfields[2]),1,1,1,1,1,1)
    return peaks


# ------------------------------------
# Main function
# ------------------------------------
def main():
    usage = "usage: %prog <-d path> [options] <bed files> ..."
    description = "Draw conservation plot for many bed files."

    optparser = OptionParser(version="%prog 0.1",description=description,usage=usage,add_help_option=False)
    optparser.add_option('-H','--height', dest='height',type='int',default=10, help="height of plot")
    optparser.add_option('-W','--width',dest='width',type='int',default=10, help="width of plot")
    optparser.add_option('-w',dest='w',type='int',default=1000, help="window width centered at middle of bed regions,default: 1000")
    optparser.add_option('-t','--title',dest='title',help="title of the figure. Default: 'Average Phastcons around the Center of Sites'",default= 'Average Phastcons around the Center of Sites')
    optparser.add_option('-d','--phasdb',dest='phasdb',help= 'The directory to store phastcons scores in the server')
    optparser.add_option('-o','--outimg',dest='outimg',help= 'output image file prefix')
    optparser.add_option("-l","--bed-label",dest="bedlabel",type="string",action="append",
                         help="the BED file labels in the figure. No space is allowed. This option should be used same times as -w option, and please input them in the same order as BED files. default: will use the BED file filename as labels.")
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")

    (options,bedfiles) = optparser.parse_args()
    options.pf_res = options.w / 100 # get 100 points to plot
    options.w = options.pf_res * 100 # trim

    bedfiles = map(os.path.abspath,bedfiles)
    bedfilenames = map(os.path.basename,bedfiles)

    bedfilenum = len(bedfiles)

    if bedfilenum < 1 or not options.phasdb:
        optparser.print_help()
        sys.exit(1)

    if options.bedlabel and len(options.bedlabel) == bedfilenum:
        bedlabel = options.bedlabel
    else:                               # or use the filename
        bedlabel = map(lambda x:os.path.basename(x),bedfiles)

    if options.height < 10:
        error("Height can not be lower than 10!")
        sys.exit(1)
    if options.width < 10:
        error("Width can not be smaller than 10!")
        sys.exit(1)

    # check the files
    for f in bedfiles:
        if not os.path.isfile(f):
            error("%s is not valid!" % f)
            sys.exit(1)

    # check phastcons db
    if not os.path.isfile(options.phasdb):
        error("%s is not valid!" % options.phasdb)
        sys.exit(1)

    if not options.phasdb:
        error("%s has no valid phastcons db bw files!" % options.phasdb)
        sys.exit(1)

    info("number of bed files: %d" % bedfilenum)

    avgValues = []
    
    # for each bed file
    for f in bedfiles:
        info("extract phastcons scores using %s" % f)
        scores = extract_phastcons(f, options.phasdb, options.w, options.pf_res)
        avgValues.append(scores)
    if options.w == 4000:
        ## 100 points for 4000, 40bp resolution
        print("\t".join([str(avgValues[0][i]) for i in [12,25,38,50,62,75,88]]))
    elif options.w == 400:
        print("\t".join([str(avgValues[0][i]) for i in [45,48,50,52,55]]))
    olddir = os.path.dirname(options.outimg)
    makeBmpFile(avgValues,olddir, options.outimg ,options.height,options.width,options.w,options.pf_res,options.title,bedlabel)

def extract_phastcons ( bedfile, phasdb, width, pf_res ):
    """Extract phastcons scores from a bed file.

    Return the average scores
    """
    info("read bed file...")
    bfhd = open(bedfile)
    bed = parse_BED(bfhd)

    # calculate the middle point of bed regions then extend left and right by 1/2 width
    bchrs = bed.peaks.keys()
    bchrs.sort()

    sumscores = []
    bw = BigWigFile(open(phasdb, 'rb'))
    for chrom in bchrs:
        info("processing chromosome: %s" %chrom)
        pchrom = bed.peaks[chrom]
        for i in range(len(pchrom)):
            mid = int((pchrom[i][0]+pchrom[i][1])/2)
            left = int(mid - width/2)
            right = int(mid + width/2)

            if left < 0:
                left = 0
                right = width
            summarize = bw.summarize(chrom, left, right, width/pf_res)
            if not summarize:
                continue
            dat = summarize.sum_data / summarize.valid_count
            sumscores.append(dat)

    ## a list with each element is a list of conservation score at the same coordinate
    sumscores = map(list, zip(*sumscores))

    ## exclude na
    sumscores = [[t2 for t2 in t if not math.isnan(t2)] for t in sumscores]
    try:
        conscores = [sum(t)/len(t) for t in sumscores]
    except ZeroDivisionError:
        conscores = [0] * (width/pf_res)

    return conscores

def makeBmpFile(avgValues, wd, outimg, h,w, width, pf_res, title, bedlabel):

    #creating R file in which to write the rscript which defines the correlation plot
    #create and save the file in the current working directory

    ## outimg should be id/prefix, that is, conf.prefix
    fileName = os.path.join(wd, os.path.basename(outimg))
    rFile = open(fileName+'.R','w')
    bmpname = fileName+'.pdf'
    rscript = 'sink(file=file("/dev/null", "w"), type="message")\n'
    rscript += 'sink(file=file("/dev/null", "w"), type="output")\n'
    rscript += 'pdf("%s",height=%d,width=%d)\n' %(bmpname,h,w)
    xInfo = range(int(-width/2),int(width/2), pf_res)
    rscript += 'x<-c('+','.join(map(str,xInfo[:-1]))+')\n' # throw the last point which may be buggy
    for i in range(len(avgValues)):
        avgscores = avgValues[i]
        tmpname = 'y'+str(i)
        rscript += tmpname+'<-c('+','.join(map(str,avgscores[:-1]))+')\n' # throw the last point which may be buggy

    tmplist = []
    for i in range(len(avgValues)):
        tmplist.append( "y%d" % i )

    rscript += "ymax <- max("+ ",".join(tmplist) +")\n"
    rscript += "ymin <- min("+ ",".join(tmplist) +")\n"
    rscript += "yquart <- (ymax-ymin)/4\n"

    rscript += 'plot(x,y0,type="l",col=rainbow(%d)[1],main=\"%s\",xlab="Distance from the Center (bp)",ylab="Average Phastcons",ylim=c(ymin-yquart,ymax+yquart))\n' % (len(avgValues),title)
    for i in range(1,len(avgValues)):
        rscript += 'lines(x,y'+str(i)+',col=rainbow(%d)[%d])\n' % (len(avgValues),i+1)
    rscript += 'abline(v=0)\n'
    legend_list = map(lambda x:"'"+x+"'", bedlabel)
    rscript += 'legend("topright",c(%s),col=rainbow(%d),lty=c(%s))\n' % (','.join(legend_list),len(avgValues),','.join(['1']*len(avgValues)))

    rscript += 'dev.off()\n'
    rFile.write(rscript)
    rFile.close()
    #executing the R file and forming the pdf file
    data = subprocess.call(['Rscript',fileName+'.R'])

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me! ;-) See you!\n")
        sys.exit(0)

#!/usr/bin/env python2.7
"""
used to calc the score of each gene to be regulated by factor.
1. For each refseq gene in genome, input a distance (for example 100kb), then I get the peak center within 100kb from gene TSS.
2. filter the peaks by p-value < 1e-5 from MACS, and only get top 10,000 peaks if it's more than 10,000
3. Then calculate a sum of 'Score' for each gene use this formula:
  Sg = lambda ldx: sum([math.exp(-0.5-4*t) for t in ldx])
4. output is in bed format. the 5th column is score.

input bed dile should 'chr1' instead of 'chrI'
## Modified to be compatible with all species
"""
import sys, os, time
import sqlite3, math
from optparse import OptionParser

CHROM_CONVERT = {'chrI':'chr1','chrII':'chr2','chrIII':'chr3','chrIV':'chr4','chrV':'chr5','chrVI':'chr6',
                 'chrVII':'chr7','chrVIII':'chr8','chrIX':'chr9','chrX':'chr10','chrXI':'chr11','chrXII':'chr12',
                 'chrXIII':'chr13','chrXIV':'chr14','chrXV':'chr15','chrXVI':'chr16','chrXVII':'chr17','chrXVIII':'chr18',
                 'chrXIX':'chr19','chrXX':'chr20','chrXXI':'chr21','chrXXII':'chr22','chrX':'chrX','chrY':'chrY','chrM':'chrM'}

#Score calc function
Sg = lambda ldx: sum([math.exp(-0.5-4*t) for t in ldx])

# print current time and information on screen
def Info(infoStr):
    print("[%s] %s" %(time.strftime('%H:%M:%S'), infoStr))

def prepare_optparser():
    """
    Prepare optparser object and validation.
    New options will be added in this function first.
    """
    usage = "usage: %prog ......"
    description = "Input a peak file, and ...... "
    optparser = OptionParser(version="%prog v1.00", description=description, usage=usage, add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")
    optparser.add_option("-t","--treat",dest="peakfile",type="string",
                         help="Input the MACS's result peak file(.bed), it will recognize it by extension.")
    optparser.add_option("-n","--name",dest="name",type="string",
                         help="this argument is used to name the result file.")
    optparser.add_option("-d","--distance", dest="distance", type="int",
                         help="Set a number which unit is 'base'. It will get peaks within this distance from gene TSS. default:100000 (100kb)", default=100000)
    optparser.add_option("-g","--genome",dest="genome",type="string",
                         help="Select a genome file (sqlite3 file) to search refGenes.")
    optparser.add_option("--top",dest="top",type="float",
                         help="input a number between 0-1, so that the script will only output a percentage of top genes.\
                               input a number bigger than 1, for example, 2000. so that the script will only output top 2000 genes. Default: output All ~30000 genes", default = -1)

    (options,args) = optparser.parse_args()
    if not options.peakfile and not options.genome and not options.name:
        optparser.print_help()
        sys.exit(1)
    if not os.path.isfile(options.peakfile):
        Info('ERROR: Cannot find peak file, a tab-peak file must be given through -t (--treat).')
        sys.exit(1)

    if not os.path.isfile(options.genome):
        Info("ERROR: Genome file not found! A annottion file must be given through -g (--genome).")
        sys.exit(1)

    if not options.name:
        options.name = os.path.splitext(options.peakfile)[0] + "_result"

    if options.top > 1:
        options.output_flag = "number"
    elif 0 < options.top < 1:
        options.output_flag = "percentage"
    elif options.top == -1:
        options.output_flag = "all"
    else:
        Info("--top options error, please set a number 0-1 or bigger than 1")
        sys.exit(1)

    # print arguments
    Info("Argument List: ")
    Info("Name = " + options.name)
    Info("peak file = " + options.peakfile)
    Info("distance = %d bp" %options.distance)
    Info("genome = %s" %options.genome)
    Info("top = %f" %options.top)
    print
    return options

class PScore:
    # connect to sqlite, select genome.
    def __init__(self, options):
        self.peakfile = options.peakfile
        self.genome = options.genome
        self.db = open(self.genome)
        self.c = self.db.readlines()
        self.geneInfo = []
        self.output_flag = options.output_flag
        self.top = options.top
        self.opts_string = "# Argument List:\n" +\
                           "# Name = %s\n" %options.name +\
                           "# peak file = %s\n" %options.peakfile +\
                           "# distance = %d bp\n" %options.distance +\
                           "# top = %f\n" %options.top
        self.peaklist = {}

    def readfile(self): #reads the file and returns a vector: each element is a bed_row.
        peakf = open(self.peakfile)
        count = 0
        self.peaklist = {}
        for line in peakf:
            if line.startswith("#") or not line.strip(): #skip "#" lines and empty lines
                continue
            line = line.split() #.bed-> 0:chrom 1:pStart 2:pEnd 3:peakName 4:-10*log10(pvalue)
            line = [line[0], int(line[1]), int(line[2]), line[3], float(line[4])]
            try:
                line[0] = CHROM_CONVERT[line[0]]
            except KeyError:
                pass
            try:
                self.peaklist[line[0]].append(line)
            except KeyError:
                self.peaklist[line[0]] = [line]
            count += 1
        peakf.close()
        for i in self.peaklist.keys():
            self.peaklist[i].sort()
        Info("Read file <%s> OK! All <%d> peaks." %(self.peakfile, count))

    def ScoreCalc(self, distance):
        #get at most 10k sites, p-value < 1e-5. each gene's regulatory potential sg = ...
        # downloaded ucsc refseq from selected fields http://genome.ucsc.edu/cgi-bin/hgTables
        for line in self.c:
            if line.startswith("#"):
                continue
            line = line.strip().split()
            self.geneInfo.append([line[2], line[1], int(line[4]), int(line[5]), line[3], line[-4]]) # (0:chrom, 1:name, 2:txStart, 3:txEnd, 4:strand, 5:symbol)

        self.geneInfo.sort()
        self.geneInfo = [list(t) for t in self.geneInfo]

        count = 0
        pvaluecutoff = 5
        for igene in self.geneInfo:
            if igene[4] == '+':
                gTSS = igene[2]
            elif igene[4] == '-':
                gTSS = igene[3]
            try:
                peaks = self.peaklist[igene[0]]
            except KeyError:
                peaks = []
            #peaksInDistance = [t+[abs((t[1]+t[2])/2-gTSS)] for t in peaks if t[4]>5 and abs((t[1]+t[2])/2-gTSS) < distance ]
            #peaksInDistance.sort(key=lambda x:x[-1])
            peaksInDistance = [abs((t[1]+t[2])/2-gTSS)*1.0/distance for t in peaks if t[4]>pvaluecutoff and abs((t[1]+t[2])/2-gTSS) < distance ] #peak pvalue<1e-5, and distance < distance
            peaksInDistance.sort()
            if len(peaksInDistance) > 10000: # extract no more than 10k peaks
                peaksInDistance = peaksInDistance[:10000]
            #score = sum(math.exp(-0.5-4*t[-1]) for t in peaksInDistance)
            igene.append(Sg(peaksInDistance))
            count += 1
            if not count % 2000:
                Info('Process <%d> genes'%count)
        self.geneInfo.sort(key=lambda x:x[-1], reverse=True)

    def Output2File(self, name):
        #print self.output_flag
        if self.output_flag == "percentage":
            number = int(len(self.geneInfo)*self.top)
            sub_geneInfo = self.geneInfo[:number]
            #print number
        elif self.output_flag == "number":
            number = int(self.top)
            sub_geneInfo = self.geneInfo[:number]
        elif self.output_flag == "all":
            sub_geneInfo = self.geneInfo[:]

        outf = open("%s_gene_score.txt"%name, "w")
        outf.write(self.opts_string)
        outf.write('#chrom\ttxStart\ttxEnd\trefseq\tscore\tstrand\tsymbol\n')
        for line in sub_geneInfo:
            outf.write('%s\t%d\t%d\t%s\t%.3f\t%s\t%s\n'%(
                       line[0], line[2], line[3], line[1], line[6], line[4], line[5]))
        outf.close()
        Info("Finished! result output to <%s_gene_score.txt>"%name)

def main():
    opts=prepare_optparser()
    g = PScore(opts)
    g.readfile()
    g.ScoreCalc(opts.distance)
    g.Output2File(opts.name)

if __name__ == "__main__":
    main()

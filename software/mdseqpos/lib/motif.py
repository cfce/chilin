"""
SYNOPSIS: This module contains the Motif and MotifList classes.  Motif stores
and processes motif information.  MotifList is a simple list wrapper for
lists--it contains some helper functions that operate over a list of motifs.
"""
import os
import sys
import re
import math
import json
from xml.dom.minidom import parse, parseString

#NOTE: when seqscan is modified to take pssms that aren't numpy arrays, then i can remove this!
import numpy

try:
    import _seq
except:
    import mdseqpos._seq as _seq

import count
from database import MotifParser
import Prob

def calc_pssm(sequences):
    """Given a list of sequences, e.g.
    [['A','C','C','A','A','T','T'],
     ['A','A','C','C','T','A','T'],
     ...
     ]
     this function returns the pssm describing the sequences.
     """
    #Step 1: For each COLUMN, count the number of A's in that pos, C's etc
    pssm = []
    if sequences and (len(sequences) > 0):
        for col in range(len(sequences[0])):
            tmp = [0, 0, 0, 0]
            for row in range(len(sequences)):
                if sequences[row][col] == 'A':
                    tmp[0] += 1
                elif sequences[row][col] == 'C':
                    tmp[1] += 1
                elif sequences[row][col] == 'G':
                    tmp[2] += 1
                else: #T
                    tmp[3] += 1
            #Step 2: normalize the counts
            hit_sum = sum(tmp)
            pssm.append(map(lambda ct: float(ct)/float(hit_sum), tmp))
    return pssm

def revcomp(s):
    """Returns the reverse compliment of a sequence"""
    t = s[:]
    t.reverse()
    rc = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N'}
    t = [rc[x] for x in t]
    return t


class DuplicateMotif(Exception):
    pass

class BadPssm(Exception):
    pass

class Motif:
    """Class for sequence motifs.
    (NOT SO) Required Attributes:
    NOTE: in ying's old Motif class, these could be None
    self.id (string): A string containing a unique ID identifying the motif.
    self.pssm (float array, N x 4 cols):
    position specific scoring matrix (PSSM) for the motif. The array
    should contain N rows, each corresponding to a single position in
    the PSSM, and 4 columns, each corresponding to the 4 DNA basepairs
    A, C, G, T, in order.  Each element in the matrix should be a
    float, and each row of the matrix should sum to 1.0.

    self.seqpos_results (dict): {'numhits','cutoff','zscore','meanposition',
    'pvalue}
    NOTE: we have convenience fns, e.g. self.getpvalue() returns
    self.seqpos_results['pvalue']
    self.antisense (bool): whether the motif is on the antisense strand
    

    Optional Attributes (descriptive):
    self.source (string): the source of the motif, e.g. 'Transfac' or 'PBM'
    self.sourcefile (string): the path to the source file of the motif
    self.species (list of string): list of species names, eg ['mouse', 'human']
    self.fullname (string): fullname of the motif    
    self.pmid (integer): the Pubmed ID of the article in which the motif was
        published
    self.numseqs(integer): number of sequences used to derive the motif's
        position specific scoring matrix (PSSM).
    self.factors (list of string): list of the transcription factors that
        bind to the motif, e.g. ['ERa', 'Runx2'].
    self.entrezs (list of integers): list of integers representing the
        Entrez Gene IDs of the transcription factor that bind to the motif,
        such as [6908, 2099].
    self.refseqs (list of strings): list of the refseq strings, e.g.
        ['NM_001141919', 'NM_013464']
    self.dbd (string): the full name of the DNA binding domain of the
        transcription factor that binds to the motif, such as
        'estrogen receptor region C'.
    """
    _ATTRIBUTES = ['id', 'pssm', 'seqpos_results', 'antisense',
                   'source', 'sourcefile', 'species', 'fullname', 'pmid',
                   'numseqs', 'factors', 'entrezs', 'refseqs', 'dbd']
    _results_fields = ['numhits','cutoff','zscore','meanposition','pvalue']

    #maybe i'll make pssm the only required attribute
    def __init__(self):
        """initializes all attributes to None"""
        for attr in Motif._ATTRIBUTES:
            setattr(self, attr, None)
        #set the seqpos_result field to None
        #and convenience accessors/setters, e.g. self.getpvalue
        self.seqpos_results = {}
        for attr in Motif._results_fields:
            self.seqpos_results[attr] = None
            def curry(x):
                return lambda : self.seqpos_results[x]
            #setattr(self, "get"+attr, lambda : self.seqpos_results[attr])
            setattr(self, "get"+attr, curry(attr))
            #setattr(self, "set"+attr, lambda x: self.seqpos_results[attr] = x)

                
    @staticmethod
    def _validpssm(pssm):
        """Checks to see if its a valid pssm:
        1. each row must have 4 cols of floats
        2. each row must sum to 1.0
        """
        def feq(f1, f2):
            epsilon = 1e-6
            return abs(f1 - f2) < epsilon

        if pssm is None: raise BadPssm("None is not a valid pssm")
        if len(pssm) == 0: raise BadPssm("Empty pssm")
        
        for i, row in enumerate(pssm):
            if len(row) != 4:
                raise BadPssm("row %s is not length 4" % i)
            else:
                if not feq(sum(row), 1.0):
                    raise BadPssm("row %s sums to %s NOT 1.0" % (i, sum(row)))
        return True
                
    def setpssm(self, pssm):
        """The preferred method to set the motif's pssm. does error checking"""
        if Motif._validpssm(pssm):
            self.pssm = pssm

    def equals(self, other):
        """Returns true if self.id == other.id"""
        return (self.id == other.id)

    def seqpos_stat(self, start, end, score, seq_length):
        """Calculates the statistics (and other things) that go into
        seqpos_results.
        numhits - number of cumulative hits
        meanposition - mean of the positions where the motif was found
        cutoff - ???
        zscore - motif's zscore
        pvalue - motif's pvalue:= max(normal cum dist(zscore), 1e-30)
        """
        plotPosWin = 200 #sliding window used to calculate mean pos plot
        numBins = 50
        MINHITS = 20
        MUALPHA0 =  0.296
        MUALPHA1 =  0.641
        SDALPHA0 =  0.943
        SDALPHA1 = -0.111
        mu_bias = 0.5
        offset = len(self.pssm)

        #calculate relative positions of where the motif is found in sequence
        frac = [abs((float((s + e - offset))/(seq_length - offset)) - 1.0)
                for (s,e) in zip(start, end)]
        #generate list of indices sorted by largest score, e.g. idx[0] is
        #the index of the largest score and score[idx[0]] is the largest score
        #print frac
        idx = [i for (i,s) in
               sorted(enumerate(score), key=lambda x: x[1], reverse=True)]
        #print idx

        cumfrac = 0
        zscore_min = 99 #infinity
        kmin, numhits, cuthits, meanpos = -1, 0, 0, 0
        
        #calculate the sums of the relative positions of sites (cumfrac)
        #starting w/ highest score on downward
        
        #NOTE: this is NOT an error, numhits is the number of sites seen
        #so far--ie. simply the index, and i is the score/fracpos index
        for numhits, i in enumerate(idx):
            if frac[i] > 1.0 or frac[i] < 0.0:                
                raise(ValueError,"Fraction out of bounds %4.2f" % frac[i] )

            cumfrac += frac[i]
            if numhits > MINHITS: #a valid number of sites have been computed
                #NOTE: 12.0 is the magic number according to cliff, its
                #the expected std dev, and it "comes from the variance of a
                #convolution of uniform distributions"
                zscore = (cumfrac - numhits*mu_bias)/(math.sqrt(numhits/12.0))

                #take most significan site
                if zscore < zscore_min:
                    zscore_min = zscore
                    kmin = score[i]
                    cuthits = numhits
                    meanpos = cumfrac/numhits - 0.5
        frac_sort = [frac[i] for i in idx]
        #calculate a mean sliding average using plotPosWin as the window size
        #normalize using (seq_length-offset)*0.5
        if len(frac_sort) <= plotPosWin:
            tplotPosWin = len(frac_sort) - 1
        else:
            tplotPosWin = plotPosWin
        meanSldAvg=[sum(frac_sort[i:i+tplotPosWin])/tplotPosWin*((seq_length-offset)*0.5)\
                        for i in range(len(frac_sort) - tplotPosWin)]

        #NOTE: downscale the number of means (in meanSldAvg) b/c plotting
        #the full thing will just kill highcharts.js on the front end

        #sample numBins = 50 at regular intervals along meanSldAvg
        chunk = len(meanSldAvg)/numBins
        binAvg = [meanSldAvg[chunk*i] for i in range(numBins-1)]
        #last one:
        binAvg.append(meanSldAvg[-1])

        if cuthits > MINHITS:
            # NOTE: w/ the pvalue, we are assuming the normal distribution
            # get the raw p-value score
#BUG?        z_mean_adj = - (MUALPHA0 + MUALPHA1*math.log(math.log(numhits)))
#BUG?        z_sd_adj = SDALPHA0 + SDALPHA1*math.log(math.log(numhits))
            z_mean_adj = - (MUALPHA0 + MUALPHA1*math.log(math.log(cuthits)))
            z_sd_adj = SDALPHA0 + SDALPHA1*math.log(math.log(cuthits))
            zscore_min_adjusted = ( zscore_min - z_mean_adj )/z_sd_adj
            pvalue = max(Prob.normal_cdf(zscore_min_adjusted), 1e-30)
        else:
            zscore_min_adjusted = 0
            pvalue = 1.0

        self.seqpos_results = {'numhits': cuthits,
                               'meanposition': meanpos,
                               'cutoff': kmin,
                               'zscore': zscore_min_adjusted,
                               'pvalue': pvalue,
                               'plot': {'chunk':chunk,
                                        'width':seq_length,
                                        'bin_avg': binAvg}}
        #import func; func.debug(locals())
        

    def seqpos(self, chip_regions, width=600, margin=50):
        """Score motif on how centrally located they are within each ChIP
        region.

        ChIP regions should be given as a ChipRegions object.
        The results of SeqPos will be stored as properties of
        self.seqpos_results.
        
        Adapted from Cliff Meyer's ChIP_region.CentralMotifScan() method."""
        ANTISENSE = 1
        MOTIFMIN = 1e-3
        
        if not chip_regions.preprocessed_regions:
            chip_regions.preprocess(width, margin)
            #process the chip-regions
            chip_regions.read_sequence(True)
            
        bgseqprob_mat = count.count(chip_regions.sequence)
        markov_order = 2
        prob_option = _seq.MAX_OPTION

        #GRR: THIS IS THE ONLY numpy dependency, and it is b/c seqscan
        #expects self.pssm to be a numpy array!
        pssm = numpy.array(self.pssm, float)
        s_idx, start, end, orient, score = \
               _seq.seqscan(chip_regions.sequence, pssm, bgseqprob_mat,
                            markov_order, prob_option)
        #adjust score
        adj_score = map(lambda s: math.log(s + MOTIFMIN), score)

        #calculate the seqpos_results (stats)
        self.seqpos_stat(start, end, adj_score, width + margin)

        #generating the observed pssm
        #fracpos is the fractional position of each site/hit
        fracpos = [abs(0.5*(s + e) - (margin + width)/2.0) / (width/2.0)
                   for (s,e) in zip(start, end)]
        #retrieve sequences whose fracposition is in (0.0, 1.0]
        seq,dis = [],[]
        for j,elem in enumerate(fracpos):
            if elem <= 1.0:
                t = list(chip_regions.sequence[int(s_idx[j])])
                dis.append(int(start[j])-len(t)/2)
                t = t[int(start[j]):int(end[j])]
                if orient[j] == ANTISENSE:
                    seq.append(revcomp(t))
                else:
                    seq.append(t)
        self.seqpos_results['pssm'] = calc_pssm(seq)
        self.seqpos_results['seq'] = ["".join(t) for t in seq]
        self.seqpos_results['dis'] = dis
        #if self.id == "MA0104":
        #    import func
        #    func.debug(locals())
                                    
    @staticmethod
    def from_flat_file(path):
        """Read in a motif from a text file which describes the pssm, and
        returns a Motif object- the pssm is in the following format, example:
        [[0.1, 0.1, 0.7, 0.1],
        [0.1, 0.7, 0.1, 0.1],
        [0.7, 0.1, 0.1, 0.1],
        [0.1, 0.1, 0.1, 0.7]]
        """
        #NOTE: this inner fn is ugly, maybe i'll just lambda it
        def pre_process(input_txt):
            """Converts weird literals, e.g. '__ob__' --> '[' according to a
            dictionary of symbol translations.  For example, galaxy converts
            '[' --> '__ob__', this fn will convert them back"""
            sym_dict = [("__ob__", "["), ("__cb__", "]")]
            for s in sym_dict:
                input_txt = input_txt.replace(s[0], s[1])
            return input_txt
            
        tmp = Motif()
        f = open(path)
        try:
            #flatten file into a string
            pssm_txt = "\n".join([line for line in f])
            pssm_txt = pre_process(pssm_txt)
            # use regular expression to pick out the rows
            row_pattern = '\[\s*(0\.\d+)\s*\,\s*(0\.\d+)\s*\,\s*(0\.\d+)\s*\,\s*(0\\.\d+)\s*\]'
            #function that converts a tuple of float strings to float objs
            convert_tpl = lambda x: [float(i) for i in x]
            pssm = [convert_tpl(r) for r in re.findall(row_pattern, pssm_txt)]
            tmp.setpssm(pssm)
        finally:
            f.close()
        
        return tmp

    @staticmethod
    def pssm_from_xml(dom_node):
        """Reads in the xml data, and returns a Nx4 col matrix (hopefully).
        <pssm><pos num='1'><A/><C/><G/><T/></pos>...</pssm>
        """
        pos = dom_node.getElementsByTagName("pos")

        tags = ["A", "C", "T", "G"]
        pssm = []
        for p in pos:
            num = p.getAttribute("num")
            row = [float(Motif.getText(p.getElementsByTagName(tag)[0].childNodes))
                   for tag in tags]
            pssm.append((num, row))
        #can pos-s be out of order?--should we order by "num"??
        pssm = sorted(pssm, key=lambda p: p[0])
        return [p[1] for p in pssm]

    @staticmethod #???
    def getText(nodelist):
        """RIPPED from: http://docs.python.org/library/xml.dom.minidom.html
        """
        rc = []
        for node in nodelist:
            if node.nodeType == node.TEXT_NODE:
                rc.append(node.data)
        return (''.join(rc)).strip()

    #NOTE: the following method is obsolete b/c we should use MotifParser 
    #module to parse db xml!
    @staticmethod
    def from_xml(dom_node):
        """
        Tries to create a new motif object that is defined in the node xml
        tree.  EXPECTED XML struct:
           <motif id='M00347'>
              <source>...</source>
              <sourcefile>...</sourcefile>
              <fullname>...</fullname>
              <pmid>...</pmid>
              <pssm><pos num='1'><A/><C/><G/><T/></pos>...</pssm>
              <numseqs>...</numseqs>
              <specieslist><species>...</species></specieslist>
              <symbollist><symbol/>...</symbollist>
              <entrezlist><entrez/>...</entrezlist>
              <refseqlist><refseq/>...</refseqlist>
              <dbd>...</dbd>
           </motif>
        """
        #NOTE: these are DUPLICATED in to_xml--should move them to
        #Motif.literals EXCEPT for antisense
        literals = ['source', 'sourcefile', 'fullname', 'pmid', 'numseqs',
                    'dbd', 'antisense']
        lists = ['symbols', 'entrezs', 'refseqs', 'species']
        
        tmp = Motif()
        tmp.id = dom_node.getAttribute("id")

        for attr in literals:
            val = dom_node.getElementsByTagName(attr)
            if val:
                setattr(tmp, attr, Motif.getText(val[0].childNodes))

        for attr in lists:
            singular_name = attr[:-1] if attr != 'species' else attr
            ls = dom_node.getElementsByTagName(singular_name+"list")
            #these lists are all optional -- so ls might not be found
            elms_list = ls[0].getElementsByTagName(singular_name) if ls else None
            if elms_list:
                setattr(tmp, attr,
                        map(Motif.getText, [e.childNodes for e in elms_list]))

        pssm_node = dom_node.getElementsByTagName("pssm")
        if pssm_node:
            tmp.pssm = Motif.pssm_from_xml(pssm_node[0])
        
        return tmp

    @staticmethod
    def from_dict(d):
        """Returns a new motif with the values specified in the dictionary. 
        The dictionary is expected to follow the format specified by the
        MotifParser
        ref: http://cistrome.dfci.harvard.edu/trac/wiki/MotifParser
        """
        fields = ['id', 'status', 'source', 'sourcefile', 'dbd', 'pmid', ]
        lists = ['entrezs', 'refseqs', 'symbols', 'synonyms',]
        #special lists: species, numseqs
        #missing: fullname, curators 
        tmp = Motif()
        for f in fields:
            #EVERY item in d is a list, so we have to make them into values
            val = ','.join(d[f])
            setattr(tmp, f, val)

        for l in lists:
            #we want to set the value, chop off the s for the dictionary
            setattr(tmp, l, d[l[:-1]])


        tmp.pssm = numpy.array(d['pssm'][0])
        tmp.species = d['species']
        tmp.numseqs = d['numseqs']
        #NOTE: symbols is a synonym for factors...and until jian changes this
        #we will have to do the following
        tmp.factors = tmp.symbols

        return tmp


    def __str__(self):
        """How motifs should be printed out"""
        divider = "".join(["-" for x in range(79)])
        lines = ["%s: %s" % (attr,getattr(self, attr)) \
                 for attr in Motif._ATTRIBUTES]
        return "%s\n%s\n%s" % (divider, "\n".join(lines), divider)

    def pssm_to_xml(self):
        """returns a string representation of the pssm:
        <pssm>
           <pos num='0'><A>0.1</A><C>0.3</C><G>0.3</G><T>0.3</T></pos>
        </pssm>
        """
        try:
            if Motif._validpssm(self.pssm):
                tmp = ""
                for i, row in enumerate(self.pssm):
                    A, C, T, G = row
                    tmp += "<pos num=\"%s\"><A>%s</A><C>%s</C><G>%s</G><T>%s</T></pos>\n" % (i, A, C, G, T)
                return "<pssm>\n%s</pssm>\n" % tmp
        except BadPssm:
            return "<pssm/>"
        
        
    def to_xml(self, print_non_schema=False):
        """Returns the xml representation of the motif. see Motif.from_xml
        to see the xml structure we are assuming
        """
        literals = ['source', 'sourcefile', 'fullname', 'pmid', 'numseqs',
                    'dbd']
        lists = ['symbols', 'entrezs', 'refseqs', 'species']
        
        #attributes that aren't defined in the xml schema, but might be useful
        #anyway
        non_schema = ['antisense', ]#seqpos_results]
        
        def print_lit(attr):
            if getattr(self, attr):
                return "<%s>%s</%s>\n" % (attr, getattr(self, attr), attr)
            return ""
        
        def print_lists(attr):
            try:
                val = getattr(self, attr)
            except:
                val = None
            if val and len(val) > 0:
                name = attr[:-1] if attr != 'species' else attr #entrezs --> entrez
                tmp = ''.join(["<%s>%s</%s>\n" % (name,v,name) for v in val])
                return "<%s>\n%s</%s>\n" % (name+"list", tmp, name+"list")
            return ""
                
        lits = ''.join(map(print_lit, literals))
        lists = ''.join(map(print_lists, lists))
        non_s = ''
        if print_non_schema:
            non_s = ''.join(map(print_lit, non_schema))
        pssm = self.pssm_to_xml()
        id = "id=\"%s\"" % self.id if self.id else ""

        return "<motif %s>\n%s</motif>" % (id, pssm + lits + lists + non_s)

    def to_json(self):
        """Returns the json representation of the motif.
        """
        _aNumpyAr = numpy.array([])
        #GRAB all of the _ATTRIBUTES and dump them into a dictionary
        tmp = {}
        for attr in Motif._ATTRIBUTES:
            val = getattr(self, attr)
            if type(val) == type(_aNumpyAr):
                tmp[attr] = val.tolist()
            else:
                tmp[attr] = val
            
        #CHECK for logoImg
        if 'logoImg' in vars(self):
            tmp['logoImg'] = self.logoImg
            
        #Control the floating point precision, unless it is pvalue --we want to preserve
        #the precision
        def printFloats(f):
            return "%.4f" % f if abs(f) > 0.0001 else f.__str__()

        json.encoder.FLOAT_REPR = printFloats
        return json.dumps(tmp)

    #NOTE: moving Motif.logo --> motiflogo.make_logo (fn)
    
class MotifList(list):
    """Class to hold a list of unique motifs, and perform operations over them
    """
        
    def append(self, motif_obj):
        """Tries to append motif_obj to the end of the list;
        IF duplicate: raises DuplicateMotif exception
        """
        if motif_obj not in self:
            #note if you do self.append, that's infinite recursion
            unbnd_app = list.append
            unbnd_app(self, motif_obj)
        else:
            raise DuplicateMotif

    def cull(self, pvalcutoff=0.001, maxmotifs=100):
        """Filter the MotifList based on the pvalue and the max number
        of motifs allowed in the list. NOTE: this is non-destructive!
        """
        sig_motifs = filter(lambda m: m.getpvalue() <= pvalcutoff, self)
        if (maxmotifs != 0) and (len(sig_motifs) > maxmotifs):
            #sort and return top maxmotifs
            sorted_motifs = sorted(self, key = lambda elem: elem.getpvalue())
            del sorted_motifs[maxmotifs:]
            return MotifList(sorted_motifs)
        else:
            return MotifList(sig_motifs)
        
    def from_xml_file(self, path):
        """Tries to parse the xml file specified by 'path', and load in the
        motifs to the MotifList.
        Expected XML file struct:<motifs><motif>...</motif></motifs>
        """
        #NOTE: here is where we rely on jian's MotifParser class               
        p = MotifParser.MotifParser(path)
        motifs = [Motif.from_dict(p.motifs[m]) for m in p.motifs]
        self.extend(MotifList(motifs))
        
    def to_xml(self):
        """Prints out the motif list in xml format, i.e.:
        <motifs><motif>...</motif></motifs>
        """
        tmp = '\n'.join([x.to_xml() for x in self])
        #yea this might be overkill to convert it to a dom object first
        dom = parseString("<motifs>\n%s\n</motifs>\n" % tmp)
        return dom.toprettyxml()

    def to_json(self):
        """Prints out the motif list as json data
        """
        return "[%s]" % ',\n'.join([m.to_json() for m in self])

    def save_to_xml(self, path):
        """writes the motif list as an xml file"""
        #if _DEBUG: print "writing new motifs to file..."
        dir_name = os.path.dirname(path)
        if not os.path.exists(dir_name):
            os.makedirs(dir_name)

        motifs_file = open(path, 'w')
        try:
            motifs_file.write(self.to_xml())
        finally:
            motifs_file.close()

    def filterBySpecies(self, species_list):
        """Filter the list by the species in species_list"""
        mapSpecies = {'hs' : 'Homo sapiens', 'mm' : 'Mus musculus',
                      'dm' : 'Drosophila melanogaster',
                      'ce' : 'Caenorhabditis elegans',
                      'zv' : 'Danio rerio'}
        #MAP short names to long names
        sl = map(lambda s: mapSpecies[s].lower(), species_list)
        motifs = []
        for m in self:
            #take the motif IFF species = [], OR an element in species is
            #in species_list
            if m.species:
                #print m.species
                found = False
                i = 0
                while (i < len(m.species) and not found):
                    if (m.species[i].lower() in sl): found = True
                    i += 1
                if found: 
                    motifs.append(m)
            else:
                motifs.append(m)
        return MotifList(motifs)

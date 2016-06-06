import sys    
import math
#sys.path.append( '/cluster/homes/cliff/PYTHON/cistrome' )

#get complement of seq
def complement(seq):
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N':'N'}
        complseq = [complement[base] for base in seq]
        return complseq

#get reverse complement of seq
def reverse_complement(seq):
        seq = list(seq)
        seq.reverse()
        return ''.join(complement(seq))

#enumerate k-mers
def permutations(items, n):   
        if n == 0:   
            yield ''  
        else:   
            for i in range(len(items)):   
                for base in permutations(items, n - 1):   
                    yield str(items[i]) + str(base)

#sort a dictionary
def sortdict(adict):
        keys = adict.keys()
        keys.sort()
        return map(adict.get, keys)

nucleotides = ['A', 'C', 'G', 'T'] 


def count(seqs, length=3):
    #print seqs
    #for a missing key, the dict entry is initialized to zero   
    counts = {}

    #count the length-element subsequences in each sequence   
    for k in range(1, length+1):
        #this may be needed if there might be zero count for some k-mer
        for key in permutations(nucleotides, int(k)):
            counts[key] = 0
        for i in seqs:   
            i = str(i).upper()
            ri = reverse_complement(i)
            for n in range(200, len(i) - k - 199):   
                if i[n : n+k].find('N') == -1:
                    counts[i[n : n + k]] += 1
                if ri[n : n+k].find('N') == -1:    
                    counts[ri[n : n + k]] += 1
  
    #print out the sequences that occur more than once   
    #for count in counts:   
    #        print ''.join(count), counts[count]  

    #sort dictionary alphabetically and store into cond
    cond = sortdict(counts)
    #print cond

    #calculate conditional probabilities
    bgprob = [ ]
    cond2g1 = [ ]
    cond3g2 = [ ]
    tot1 = 0.0
    for i in range(4):
        ii = i * (pow(4 , int(length)) - 1) / 3
        bgprob.append(cond[ii])
        tot1 += cond[ii] + 1
        tot2 = 0.0
        for j in range(4):
            jj = ii + 1 + j * 5
            tot2 += cond[jj] + 1
        for j in range(4):
            jj = ii + 1 + j * 5
            cond2g1.append(int(100000 * math.log((cond[jj] + 1) / tot2)))
            tot3 = 0.0
            for k in range(4):
                tot3 += cond[jj + k + 1] + 1
            for k in range(4):    
                cond3g2.append(int(100000 * math.log((cond[jj + k + 1] + 1)/ tot3)))
        
    for i in range(4):
        bgprob[i] = int(100000 * math.log((bgprob[i] + 1)/ tot1))

    for i in range(4):
        for j in range(4):
            if j < i:
                tmp = cond2g1[4 * i + j]
                cond2g1[4 * i + j] = cond2g1[4 * j + i]
                cond2g1[4 * j + i] = tmp
            for k in range(4):
                if k < i:
                    tmp = cond3g2[16 * i + 4 * j + k] 
                    cond3g2[16 * i + 4 * j + k] = cond3g2[16 * k + 4 * j +i]
                    cond3g2[16 * k + 4 * j + i] = tmp
    return bgprob + cond2g1 + cond3g2


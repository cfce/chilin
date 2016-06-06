#!/usr/bin/env python
#coding: utf-8

"""Reference: http://nar.oxfordjournals.org/content/early/2013/12/23/nar.gkt1302.full

Created by: Jian Ma 2014-03-15
Modiied by: Jian Ma 2014-04-14
"""
import math
import numpy
import mdseqpos
#import MotifParser as mp

def sum_IC(m1):
    """sum of IC of each column in the matrix
    """
    ic = [IC(t) for t in m1]
    return sum(ic)

def IC(v1):
    """IC of a vector
    """
    total = sum(v1)
    v1 = [t * 1.0 / total for t in v1]
    return 2 + sum([t * math.log(t, 2) for t in v1])

def pcc_vector(v1, v2):
    """Pearson Correlation Coefficient for 2 vectors
    """
    len1 = len(v1)
    len2 = len(v2)
    if len1 != len2:
        return None
    else:
        length = len1
    avg1 = 1.0 * sum(v1) / len(v1)
    avg2 = 1.0 * sum(v2) / len(v2)
    dxy = [(v1[i] - avg1) * (v2[i] - avg2) for i in range(length)]
    dx2 = [(v1[i] - avg1) ** 2 for i in range(length)]
    dy2 = [(v2[i] - avg2) ** 2 for i in range(length)]
    return sum(dxy) / (sum(dx2) * sum(dy2)) ** 0.5

def pcc_matrix_IC(m1, m2):
    """IC weighted pcc score of each column, and get a sum
    
    The sharp of m1 and m2 should be same
    """
    if len(m1) != len(m2):
        return None
    weighted_pcc = 0
    for pos in range(len(m1)):
        ic_x = IC(m1[pos])
        ic_y = IC(m2[pos])
        #print ic_x, ic_y
        ic_xy = (ic_x * ic_y) ** 0.5
        pcc = pcc_vector(m1[pos], m2[pos])
        
        #if ic_x >= 0.2 and ic_y >= 0.2:
        #    ic = ic_xy
        #elif ic_x >= 0.2:
        #    ic = ic_x
        #elif ic_y >= 0.2:
        #    ic = ic_y
        #else:
        #    ic = 0
        #weighted_pcc += ic * pcc
        weighted_pcc += ic_xy * pcc
    return weighted_pcc

def similarity(m1, m2):
    len1 = len(m1)
    len2 = len(m2)
    if len1 < len2:
        m1, m2 = m2, m1 #ensure m1 > m2
        len1, len2 = len2, len1
        #print range(0, len2 - len1 + 1)
    max_score = 0
    shift_pos = None
    reverse = None
    for shift in range(1-len2, len1):
        if shift < 0:
            m1new = m1[:len2+shift]
            m2new = m2[-shift:]
            m2rev = [t[::-1] for t in m2[::-1]][-shift:] # important: reverse first, then slice
        elif shift <= len1 - len2:
            m1new = m1[shift:shift+len2]
            m2new = m2[:]
            m2rev = [t[::-1] for t in m2[::-1]][:]
        elif len1 - shift < len2:
            m1new = m1[shift:]
            m2new = m2[:len1-shift]
            m2rev = [t[::-1] for t in m2[::-1]][:len1-shift]
        
        ic1new = [IC(t) for t in m1new]
        ic2new = [IC(t) for t in m2new]
        ic2rev = [IC(t) for t in m2rev]
        
        # un-reverse
        weight_overlap = sum([(t1 * t2) ** 0.5 for t1, t2 in zip(ic1new, ic2new)])
        weight_nonoverlap = sum_IC(m1) + sum_IC(m2) - sum_IC(m1new) - sum_IC(m2new)
        weighted_pcc = pcc_matrix_IC(m1new, m2new)
        score = weighted_pcc / (weight_overlap + weight_nonoverlap)
        if score > max_score:
            max_score = score
            shift_pos = shift
            reverse = False
            
        # reverse
        weight_overlap = sum([(t1 * t2) ** 0.5 for t1, t2 in zip(ic1new, ic2rev)])
        weight_nonoverlap = sum_IC(m1) + sum_IC(m2) - sum_IC(m1new) - sum_IC(m2rev)
        weighted_pcc = pcc_matrix_IC(m1new, m2rev)
        score = weighted_pcc / (weight_overlap + weight_nonoverlap)
        if score > max_score:
            max_score = score
            shift_pos = shift
            reverse = True
    
    return max_score, shift_pos, reverse

'''
#test
inf = open('test.distance')
mm=[]
ii=0
for l in inf:
    ii+=1
    if ii%10000==0:
        print ii
    l = l.split('\t')
    id1 = l[0].split('_')[0]
    id2 = l[1].split('_')[0]
    score = float(l[2])
    cc = similarity(p.motifs[id1]['pssm'][0], p.motifs[id2]['pssm'][0])
    mm.append(cc[0]-score)
    if abs(cc[0]-score) > 0.2:
        print id1, id2, cc, score


mmm=sorted(mm, key=lambda x:abs(x))
mmm[0]
mmm[-1]
'''

class Cluster:
   def __init__(self):
       self.motif = {}
       self.nodes = []
       self.nodescount = 1
       self.score_for_node = 0.0

def flat(cluster):
    '''flat a Cluster class into a list.
    '''
    if not cluster.nodes:
        return [cluster.motif]
    else:
        return flat(cluster.nodes[0]) + flat(cluster.nodes[1])

def motif_hcluster(motif_list, cutoff):
    """Use complete distance, that mean get farest for each cluster.
    """
    #from mdseqpos import pwmclus_motif_comp as pmc
    #p = mp.MotifParser
    #keys = p.motifs.keys()
    #keys = keys[:20]
    #print keys
    clusters = [] # a list with clustered motifs.
    cluster_score_mat = numpy.array([[None for t in range(len(motif_list))] for m in range(len(motif_list))]) # 2d matrix, order as keys.
    
    # Adding motifs from mp to clusters
    #for k in keys:
    for m in motif_list:
        cl = Cluster()
        cl.motif.update(m)
        cl.motif['pssm'] = [[float(t1) for t1 in t] for t in m['pssm']]  #numpy.array(m['pssm'])
        cl.motif['id'] = m['id']
        clusters.append(cl)

    #Calc each 2 pair of motifs and fill (score, position, strand) in matrix like below
    # xxx
    #  xx
    #   x
    #    
    for i in range(len(clusters)-1):
        for j in range(i+1, len(clusters)):
            ts = similarity(clusters[i].motif['pssm'], clusters[j].motif['pssm']) # score, position, strand
            cluster_score_mat[i][j] = ts[0] # 1-distance score

    idcount = 0
    mround = 0
    
    #find index of the max score
    while 1:
        max_score = cluster_score_mat.max()
        mshape = cluster_score_mat.shape
        max_index = (999, 999)
        findit = False
        for i in range(mshape[0]):
            for j in range(i+1, mshape[1]):
                if abs(max_score - cluster_score_mat[i][j]) < 1e-5:
                    max_index = (i, j)
                    findit = True
                    break
            if findit:
                break
        
        # if the max is smaller than cutoff, end of clustering
        if max_score < cutoff:
            break
        
        # Build merged motif except the matrix, matrix will add later.
        max2 = clusters.pop(max_index[1])
        max1 = clusters.pop(max_index[0])
        merged = Cluster()
        merged.score_for_node = max_score
        merged.nodescount = max1.nodescount + max2.nodescount
        merged.nodes = [max1,max2]
        merged.motif['id'] = 'MT%d'%idcount
        
        idcount += 1
        
        # Append merged motif to cluster
        clusters.append(merged)
        print 'Finish round', mround
        mround += 1
        
        # Refine score matrix and list. Including cut 2 line and 2 col, and then add 1 line and 1 col.
        i, j = max_index
        score_y1 = numpy.append(cluster_score_mat[:,i][:i], cluster_score_mat[i][i:])
        score_y2 = numpy.append(cluster_score_mat[:,j][:j], cluster_score_mat[j][j:])
        score_ymin = []
        for x1, x2 in zip(score_y1, score_y2):
            if x1 is None or x2 is None:
                pass
            elif x1 > x2:
                score_ymin.append(x2)
            else:
                score_ymin.append(x1)
        cluster_score_mat = numpy.append(cluster_score_mat[:j], cluster_score_mat[j+1:],0)
        cluster_score_mat = numpy.append(cluster_score_mat[:i], cluster_score_mat[i+1:],0)
        cluster_score_mat = numpy.append(cluster_score_mat[:,:j], cluster_score_mat[:,j+1:], 1)
        cluster_score_mat = numpy.append(cluster_score_mat[:,:i], cluster_score_mat[:,i+1:], 1)
        x = cluster_score_mat.shape[0]
        
        if score_ymin:
            cluster_score_mat = numpy.append(cluster_score_mat, [[t] for t in score_ymin], 1)
            cluster_score_mat = numpy.append(cluster_score_mat, [[None for t in range(x+1)]], 0)
        else:
            break # all motif in it are 1 cluster.
    
    #print 'Cluster %d motifs into %d clusters' %(len(keys), len(clusters))
    flat_clusters = [flat(t) for t in clusters]
    return flat_clusters

def motif_hcluster2(motif_list, cutoff):
    """Use complete distance, that mean get farest for each cluster.
    
    This one is for the MDSeqPos.py to use, input motif_list is list of Motif() obj.
    """
    clusters = [] # a list with clustered motifs.
    cluster_score_mat = numpy.array([[None for t in range(len(motif_list))] for m in range(len(motif_list))]) # 2d matrix, order as keys.
    
    # Adding motifs from mp to clusters
    for m in motif_list:
        cl = Cluster()
        cl.motif = m
        clusters.append(cl)

    #Calc each 2 pair of motifs and fill (score, position, strand) in matrix like below
    # xxx
    #  xx
    #   x
    #    
    count = 0
    for i in range(len(clusters)-1):
        count+=1
        print count
        for j in range(i+1, len(clusters)):
            ts = similarity(clusters[i].motif.seqpos_results['pssm'], clusters[j].motif.seqpos_results['pssm']) # score, position, strand
            cluster_score_mat[i][j] = ts[0] # 1-distance score

    idcount = 0
    mround = 0
    
    #find index of the max score
    while 1:
        if not cluster_score_mat.any():
            break
        max_score = cluster_score_mat.max()
        mshape = cluster_score_mat.shape
        max_index = (999, 999)
        findit = False
        for i in range(mshape[0]):
            for j in range(i+1, mshape[1]):
                if abs(max_score - cluster_score_mat[i][j]) < 1e-5:
                    max_index = (i, j)
                    findit = True
                    break
            if findit:
                break
        
        # if the max is smaller than cutoff, end of clustering
        if max_score < cutoff:
            break
        
        # Build merged motif except the matrix, matrix will add later.
        max2 = clusters.pop(max_index[1])
        max1 = clusters.pop(max_index[0])
        merged = Cluster()
        merged.motif = mdseqpos.motif.Motif()
        merged.score_for_node = max_score
        merged.nodescount = max1.nodescount + max2.nodescount
        merged.nodes = [max1,max2]
        merged.motif.id = 'MT%d'%idcount
        
        idcount += 1
        
        # Append merged motif to cluster
        clusters.append(merged)
        print 'Clustering, finish round', mround
        mround += 1
        
        # Refine score matrix and list. Including cut 2 line and 2 col, and then add 1 line and 1 col.
        i, j = max_index
        score_y1 = numpy.append(cluster_score_mat[:,i][:i], cluster_score_mat[i][i:])
        score_y2 = numpy.append(cluster_score_mat[:,j][:j], cluster_score_mat[j][j:])
        score_ymin = []
        for x1, x2 in zip(score_y1, score_y2):
            if x1 is None or x2 is None:
                pass
            elif x1 > x2:
                score_ymin.append(x2)
            else:
                score_ymin.append(x1)
        cluster_score_mat = numpy.append(cluster_score_mat[:j], cluster_score_mat[j+1:],0)
        cluster_score_mat = numpy.append(cluster_score_mat[:i], cluster_score_mat[i+1:],0)
        cluster_score_mat = numpy.append(cluster_score_mat[:,:j], cluster_score_mat[:,j+1:], 1)
        cluster_score_mat = numpy.append(cluster_score_mat[:,:i], cluster_score_mat[:,i+1:], 1)
        x = cluster_score_mat.shape[0]
        
        if score_ymin:
            cluster_score_mat = numpy.append(cluster_score_mat, [[t] for t in score_ymin], 1)
            cluster_score_mat = numpy.append(cluster_score_mat, [[None for t in range(x+1)]], 0)
        else:
            break # all motif in it are 1 cluster.
    
    #print 'Cluster %d motifs into %d clusters' %(len(keys), len(clusters))
    flat_clusters = [flat(t) for t in clusters]
    return flat_clusters

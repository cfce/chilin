#!/usr/bin/env python
#coding: utf-8

"""Convert xml database to cluster

Require PWMclus
"""

import MotifParser as mp
import os

motifdb = 'cistrome.xml'
p = mp.MotifParser(motifdb)

print 'Create tmp file for PWMclus input'
outf = open('tmp_pwmclus.tmp', 'w')
for k in p.motifs.keys():
    outf.write('>%s\n' %(k, ))
    for line in p.motifs[k]['pssm'][0]:
        outf.write('\t'.join([str(mm) for mm in line]))
        outf.write('\n')
    outf.write('\n')
outf.close()

print 'Running PWMclus to generate clusters'
os.system('PWMclus -in tmp_pwmclus.tmp -out tmp_ -linkage c')

print 'Format output file'
count = 0
outf = open('cistrome.cluster', 'w')
for line in open('tmp_.clusters').readlines():
    cid, l = line.strip().split('\t', 1)
    l = l.split('\t')
    for il in l:
        outf.write("%s\t%s\n"%(il, cid))
    count += 1
outf.close()
print 'Clusters x%d' %count

print 'Cleaning tmp files'
os.system('rm tmp_*')



############################
# to do
# for each symbol in cistrome.xml, only keep one record if pssm are similar
# 1) do hcluster in each symbol group, 
# 2) for each group in hcluster, keep hs > mm, keep MC > MS > JASPAR > ..., 
# keep the smaller id to make sure motif do not change btween versions.
############################
'''
import mdseqpos.pwmclus_motif_comp as pmc
keys = p.motifs.keys()
symbol2key = {} #symbol2key['ATF3'] = mp
for k in keys:
    symbol = '|'.join(p.motifs[k]['symbol'])
    if symbol2key.has_key(symbol.upper()):
        symbol2key[symbol.upper()] += p.SearchMotif(id=k)
    else:
        symbol2key[symbol.upper()] = p.SearchMotif(id=k)

for each in symbol2key.values():
'''

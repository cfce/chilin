#!/usr/bin/env python

""" 
convert multiple dc json files into one
one directory has one json file
"""
import sys
reload(sys)
sys.setdefaultencoding( "utf-8" )

import glob
import re
import json
import sys, os

def sort_dc(l):
    """ process samples id to be uniform and ordered
    """
    reorder_treat = []
    reorder_control = []
    treat = []
    control = []
    for i in l:
        index=int(i.split("_rep")[1])
        if "treat" in i:
            treat.append(index-1)
            reorder_treat.append(i)
        elif "control" in i:
            control.append(index-1)
            reorder_control.append(i)
    tn, tc = 0, 0
    for t in treat:
        tn+=1
    for c in control:
        tc+=1
            
    reorder = [ reorder_treat[t] for t in treat ] + [ reorder_control[c] for c in control ]
    return (reorder, tn, tc)
            

def json_dump(json_dict, output): 
    """
    dump out uniform json files for collecting statistics
    :param json_dict: output python dict to json
    :return: json_file name
    """
    # if not (os.path.exists(output) and os.path.getsize(output) > 0):
    with open(output, "w") as f:
        json.dump(json_dict, f, indent=4)

def html_dump(content, output): 
    """ write html content
    """
    with open(output, "w") as f:
        f.write(content)

def judge(metric, cutoff):
    """ return binary pass or fail 
    """
    if metric >= cutoff:
        return 1
    return 0

def json_load(json_file):
    with open(json_file, "r") as f:
        json_dict = json.load(f)
    return json_dict

def summary_to_json(d):
    """ convert to one json
    path: dataset1195/attic/json
          dataset1195/json
    *map.json
    *fastqc.json
    *dhs.json,
    *frip.json
    *pbc.json
    *meta.json
    *macs2*.json
    """

    f = open(os.path.basename(d) + "_compiled_database.xls", 'w')
    print >>f, "DatsetId\tTreatOrControlReplicates\tReadLength\tFastQC\tRawReadsNumber\tUniquelyMappedReadsNumber\tUniquelyMappedRatio\tPBC\tFRiP\tPeaksNumber\tPeaksFoldChange>=10\tPeaksFoldChange>=20\tFragmentSize\tChIPReplicatesPairwiseLabels\tReplicatesOverlapRatio\tReplicatesWigCorrelation\tMergedPeaksNumber\tMergedPeaksFoldChange>=10\tMergedPeaksFoldChange>=20\tMergedPeaksPromoterPercentage\tMergedPeaksExonPercentage\tMergedPeaksIntronPercentage\tMergedPeaksIntergenicPercentage\tMergedPeaksUnionDHSRatio"
    with open(d) as inf:
        for line in inf:
            line = line.strip().split()
            def path_find(x):
                fin = glob.glob(os.path.join(os.path.abspath(line[1]), "attic", "json", x)) + glob.glob(os.path.join(os.path.abspath(line[1]), "json", x))
                if fin and os.path.getsize(fin[0]) > 0:
                    return fin
                return None
            if os.path.exists(line[1]):
                json_dict = {}
                explore = path_find("*map.json") ## try to search prefix by *map.json
                explore2 = path_find("*fastqc.json") ## try to search prefix by *map.json
                if explore:
                    prefix = os.path.basename('_'.join(explore[0].split('_')[:-1]))
                else:
                    if explore2:
                        prefix = os.path.basename('_'.join(explore2[0].split('_')[:-1]))
                print(prefix)
                maps = None
                if path_find(prefix+"_map.json"):
                    maps = json_load(path_find(prefix+"_map.json")[0])
                fastqs = None
                if path_find(prefix+"_fastqc.json"):
                    fastqs = json_load(path_find(prefix+"_fastqc.json")[0])

                if maps:
                    samples = maps['stat'].keys()
                if fastqs:
                    samples = fastqs['stat'].keys()

                # get samples labels
                samples, tn, tc = sort_dc(samples)
                if path_find(prefix+"_fastqc.json"):
                    fastqs = json_load(path_find(prefix+"_fastqc.json")[0])['stat']
                else:
                    fastqs = {}
                if path_find(prefix+"_pbc.json"):
                    pbcs = json_load(path_find(prefix+"_pbc.json")[0])['stat']
                else:
                    pbcs = {}
                if path_find(prefix+"_dhs.json"):
                    dhs = json_load(path_find(prefix+"_dhs.json")[0])['stat']
                else:
                    dhs = 0.0
                if path_find(prefix+"_frip.json"):
                    frip = json_load(path_find(prefix+"_frip.json")[0])['stat']
                else:
                    frip = {}
                if path_find(prefix+"_frag.json"):
                    frag = json_load(path_find(prefix+"_frag.json")[0])['stat']
                else:
                    frag = {}

                if path_find(prefix+"_macs2.json"):
                    mergedpeaks = json_load(path_find(prefix+"_macs2.json")[0])
                else:
                    mergedpeaks = {}

                reps = glob.glob(os.path.join(line[1], 'attic', 'json', '*_rep.json'))+glob.glob(os.path.join(line[1], 'json', '*_rep.json'))
                peaks = []
                f10 = []
                f20 = []
                FLAG = True
                pairlabels = ""
                repcors = ""
                repol = ""
                if reps:
                   for r in reps:
                       if 'macs' in r:
                           peaks = []
                           with open(r) as repf:
                               m = json.load(repf)
                           for i in samples:
                               if not 'control' in i:
                                   try:
                                       peaks.append( m['stat'][i]['totalpeak'] )
                                   except:
                                       FLAG = False
                                       break
                                   f10.append( m['stat'][i]['peaksge10'] )
                                   f20.append( m['stat'][i]['peaksge20'] )
                       else:
                           with open(r) as repf:
                               m = json.load(repf)
                               metric = m['stat']
                               meta = m['input']['overlap']
                               mets = []
                               for met in meta:
                                   mets.append(map(lambda x: str(int(x)+1),  met.replace('.overlap', '').split('_')[-2:]))

                           cor = metric['cor']
                           overlap = metric['overlap']


                   if not FLAG:
                       continue

                   for c, o, m in zip(cor, overlap, mets):
                       pairlabels += ','.join(m) + ';'
                       repcors += str(c) + ';'
                       largePeak = max([peaks[int(m[0])-1], peaks[int(m[1])-1]])
                       if largePeak > 0:
                           repol += str(float(o)/largePeak) + ';'
                       else:
                           repol += 'NA;'
                meta = path_find(prefix+"_meta.json")
                if meta:
                    for i in meta:
                        if "enrich" in i:
                            continue
                        else:
                            meta = json_load(i)
                n = 0
                for sam in samples:
                    content = ''
                    content += line[0] + '\t'
                    content += sam + '\t'
                    if fastqs:
                        content += str(fastqs.get(sam, 0)['sequence_length']) + '\t' + str(fastqs.get(sam, 0)['median']) + '\t'
                    else:
                        content += 'NA\t'
                        continue
                    if maps:
                        content += str(maps['stat'][sam]['total']) + '\t'

                        content += str(maps['stat'][sam]['mapped']) + '\t'
                        content += str(round(float(maps['stat'][sam]["mapped"])/maps['stat'][sam]["total"],4)) + '\t'
                    else:
                        content += 'NA\t'
                        content += 'NA\t'
                        content += 'NA\t'
                        continue

                    if pbcs:
                        content += str(pbcs.get(sam, 0)['PBC'])+'\t'
                    else:
                        content += 'NA\t'
                    if frip:
                        content += str(frip.get(sam, 0)['frip']) + '\t'
                    else:
                        content += 'NA\t'


                    if 'control' in sam:
                        content += 'NA\t'
                        content += 'NA\t'
                        content += 'NA\t'
                        content += 'NA\t'
                    else:
                        if tn == 1:
                            if mergedpeaks:
                                content += str(mergedpeaks['stat'].get('totalpeak', 0)) + '\t'
                                content += str(mergedpeaks['stat'].get('peaksge10', 0)) + '\t'
                                content += str(mergedpeaks['stat'].get('peaksge20', 0)) + '\t'
                            else:
                                content += 'NA\t'
                                content += 'NA\t'
                                content += 'NA\t'
                        else:
                            if peaks:
                                content += "%s\t"%(peaks[n])
                                content += "%s\t"%f10[n]
                                content += "%s\t"%f20[n]
                            else:
                                content += 'NA\t'
                                content += 'NA\t'
                                content += 'NA\t'
                            n += 1
                        if frag:
                            content += "%s\t"%(frag[sam])
                        else:
                            content += 'NA\t'

                    if reps and peaks:
                        if pairlabels:
                            content += pairlabels.rstrip(';') + '\t'
                            content += repol.rstrip(';') + '\t'
                            content += repcors.rstrip(';') + '\t'
                        else:
                            content += 'NA\t'
                            content += 'NA\t'
                            content += 'NA\t'
                    else:
                        content += 'NA\t'
                        content += 'NA\t'
                        content += 'NA\t'

                    if mergedpeaks:
                        content += str(mergedpeaks['stat'].get('totalpeak', 0)) + '\t'
                        content += str(mergedpeaks['stat'].get('peaksge10', 0)) + '\t'
                        content += str(mergedpeaks['stat'].get('peaksge20', 0)) + '\t'
                    else:
                        content += 'NA\t'
                        content += 'NA\t'
                        content += 'NA\t'
                    if meta:
                        if meta.has_key('stat'):
                            meta = meta['stat']
                        else:
                            meta = meta
                        content += str(meta['promoter']) + '\t'
                        content += str(meta['exon']) + '\t'
                        content += str(meta['intron']) + '\t'
                        content += str(meta['inter']) + '\t'
                    else:
                        content += 'NA\t'
                        content += 'NA\t'
                        content += 'NA\t'
                        content += 'NA\t'

                    if dhs:
                        if hasattr(dhs, 'get'):
                            if dhs.get('number',0)>0:
                                content += str(float(dhs.get('overlap',0))/dhs.get('number',0))
                            else:
                                content += 'NA'
                        else:
                            content += 'NA'
                    else:
                        content += 'NA'
                    print >>f, content
                f.flush()
    f.close()

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print >>sys.stderr, "input dc directory annotation files like:"
        print >>sys.stderr, "id\tdirectory"
        print >>sys.stderr, "1\tabsolute_path1"
        print >>sys.stderr, "2\tabsolute_path2"
        sys.exit(1)

    summary_to_json(sys.argv[1])


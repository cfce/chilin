#!/usr/bin/env python

""" 
convert multiple dc json files into one
one directory has one json file
"""
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
    dataset1195/attic/json
    dataset1195/json
    /data1/DC_results/hg38
    /data1/DC_results/mm10
    /data2/DC_results/hg38
    /data2/DC_results/mm10
    /data3/DC_results/hg38
    /data3/DC_results/mm10
    *map.json, *fastqc.json, *dhs.json,
    *frip.json, *pbc.json, *meta.json
    *macs2*.json
    """
    n = 0
    def path_find(x):
        fin = glob.glob(os.path.join(os.path.abspath(d), "attic", "json", x)) + glob.glob(os.path.join(os.path.abspath(d), "json", x))
        if fin and os.path.getsize(fin[0]) > 0:
            return fin
        return None
    json_dict = {}
    explore = path_find("*map.json") ## try to search prefix by *map.json
    explore2 = path_find("*fastqc.json") ## try to search prefix by *map.json

    if explore:
        prefix = os.path.basename('_'.join(explore[0].split('_')[:-1]))
    else:
        if explore2:
            prefix = os.path.basename('_'.join(explore2[0].split('_')[:-1]))
    
    print prefix
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

    samples, tn, tc = sort_dc(samples)
    if path_find(prefix+"_fastqc.json"):
        fastqs = json_load(path_find(prefix+"_fastqc.json")[0])
        fastqc_judge = any(map(judge,[fastqs['stat'][z]["median"] for z in fastqs['stat'].keys()], [25]*len(samples)))
        fastqc = [fastqs['stat'][z]["median"] for z in fastqs['stat'].keys()]
    else:
        fastqc_judge = False
        fastqc = [ 0 ]
    if path_find(prefix+"_pbc.json"):
        pbcs = json_load(path_find(prefix+"_pbc.json")[0])

    if path_find(prefix+"_dhs.json"):
        dhs = json_load(path_find(prefix+"_dhs.json")[0])
    else:
        dhs = 0.0

    if path_find(prefix+"_frip.json"):
        frip = json_load(path_find(prefix+"_frip.json")[0])
    meta = path_find(prefix+"_meta.json")
    if path_find(prefix+"_macs2.json"):
        peaks = json_load(path_find(prefix+"_macs2.json")[0])
    if meta:
        for i in meta:
            if "enrich" in i:
                e_meta = json_load(i)
            else:
                meta = json_load(i)
    if path_find(prefix+"_seqpos.json"):
        motif = json_load(path_find(prefix+"_seqpos.json")[0])
    macs_reps = []
    reps = []
    try:
        for i in glob.glob(os.path.join(os.path.abspath(d), "attic", "json", "*rep.json")):
            if "macs" in i:
                macs_reps = json_load(i)
            else:
                reps = json_load(i)
    except:
        try:
            for i in glob.glob(os.path.join(os.path.abspath(d), "json", "*rep.json")):
                if "macs" in i:
                    macs_reps = json_load(i)
                else:
                    reps = json_load(i)
        except:
            pass
   # for each metric,
   # fastqc: 25
   # uniquely mapped reads: 60%
   # uniquely locations of 4M reads: 70%
   # PBC: 80%
   # FRiP: 1%
   # Fold change peaks 10, 500
   # DHS: 70%
   # overlap ratio: 50%(in both replicates) and correlation: 0.7
   # motif: occur in top 10
    json_dict["table"] = {}
    json_dict["judge"] = {}
    if glob.glob(os.path.join(os.path.abspath(d), "json", "*pbc.json")) + glob.glob(os.path.join(os.path.abspath(d), "attic", "json", "*pbc.json")):
        pbc_judge = any(map(judge,[float(pbcs['stat'][z]["PBC"]) for z in pbcs['stat']], [0.8]*len(samples)))
        pbc = [round(float(pbcs['stat'][z]["PBC"]),3) for z in pbcs['stat']]
    else:
        pbc_judge = False
        pbc = [ 0 for z in samples]

    if glob.glob(os.path.join(os.path.abspath(d), "attic", "json", "*macs2.json"))+glob.glob(os.path.join(os.path.abspath(d), "json", "*macs2.json")):
        peaks_judge = any([judge(int(peaks['stat']["peaksge10"]), 500)])
        peaks = [ peaks['stat']["totalpeak"], peaks['stat']["peaksge10"], peaks['stat']["peaksge20"] ]
    else:
        peaks_judge = False
        peaks = [ 0, 0, 0 ]
    if maps:
        ratio_judge = any(map(judge,[float(maps['stat'][z]["mapped"])/maps['stat'][z]["total"] for z in maps['stat']], [0.6]*len(samples)))
        mapstr = maps['stat'][z]["mapped"]
        maptstr = maps['stat'][z]["total"]
    else:
        mapstr = " "
        maptstr = " "
    if glob.glob(os.path.join(os.path.abspath(d), "attic", "json", "*frip.json")) + glob.glob(os.path.join(os.path.abspath(d), "json", "*frip.json")):
        frip_judge = any(map(judge,[frip['stat'][z]["frip"] for z in frip['stat']], [0.01]*len(samples)))
        frip = [round(frip['stat'][z]["frip"],3) for z in frip['stat'].keys()]
    else:
        frip_judge = False
        frip = [ 0 for z in samples]

    try:
        dhs_judge = any([judge(float(dhs['stat']['overlap'])/dhs['stat']['number'], 0.7)])
    except:
        dhs_judge = False
    
    if len(prefix.split('_')) < 2:
       factor = prefix.split('_')[0]
    else:
       factor = prefix.split('_')[1]
    exclude_motif = [re.match(r"^H\d+\w*\d*",factor.upper()), re.match(r"^POL",factor.upper()), factor.upper()=="DNASE"]
    if any(exclude_motif):
        motif = "NA"
    else:
        try: 
            motif = factor.upper() in [motif['stat']['satisfied_motifs'][i]['factors'][0].upper() for i in range(len(motif['stat']['satisfied_motifs']))]
        except:
            motif = False

    if maps:
        ratio = [round(float(maps['stat'][z]["mapped"])/maps['stat'][z]["total"],3) for z in samples]
        map_num = [maps['stat'][z]["mapped"] for z in samples]
        raw_num = [maps['stat'][z]["total"] for z in samples]

    trans = lambda x: str(round(x, 3)*100)+"%"
    try:
        meta = meta['stat']
    except:
        meta = ''
    # table content with number
    json_dict["table"]["sample"] = samples
    json_dict["table"]["treat_number"] = tn
    json_dict["table"]["control_number"] = tc
    json_dict["table"]['fastqc'] = fastqc
    if maps:
        json_dict["table"]['map'] = map(trans, ratio)
        json_dict["table"]['map_number'] = map_num
        json_dict["table"]['raw_number'] = raw_num
    json_dict["table"]['pbc'] = map(trans, pbc)
    json_dict["table"]['peaks'] = peaks

    # qc judgement with true or false
    json_dict["judge"]  = {}
    json_dict["judge"]["fastqc"]  = fastqc_judge
    if maps:
        json_dict["judge"]['map'] = ratio_judge
    json_dict["judge"]['pbc'] = pbc_judge
    json_dict["judge"]['peaks'] = peaks_judge
    # json_dict[d.strip("dataset")]['conservation'] = 'NA'
    if reps:
        reps_judge = any(map(judge, reps["stat"]['cor'], [0.6] * len(samples)))
        json_dict["table"]['reps'] = reps["stat"]['cor']
        json_dict["judge"]['reps'] = reps_judge
        json_dict["overlap"] = ' / '.join(map(str,reps["stat"]['overlap']))
    if macs_reps:
        peaks_reps = [ ' / '.join(map(str, [macs_reps['stat'][z]["totalpeak"], macs_reps['stat'][z]["peaksge10"], macs_reps['stat'][z]["peaksge20"]])) for z in samples if macs_reps['stat'].has_key(z) ]
        json_dict["table"]['peaks_reps'] = peaks_reps
    json_dict["table"]['frip'] = map(trans, frip)
    json_dict["judge"]['frip'] = frip_judge
    
    json_dict["table"]['motif'] = motif
    json_dict["judge"]['motif_judge'] = motif

    if type(dhs) != float and dhs['stat']['number'] !=0 :
	    dhs = round(float(dhs['stat']['overlap'])/dhs['stat']['number'], 3)
    else:
        dhs = 0

    json_dict["table"]['dhs'] = str(dhs*100)+"%"
    json_dict["judge"]['dhs'] = dhs_judge
    if meta:
        json_dict["table"]['meta'] = [ ' / '.join(map(trans, [meta['promoter'],meta['exon'],meta['intron'],meta['inter']]))]
        json_dict["table"]['meta_orig'] = meta
    # if e_meta:
    #     json_dict["table"]['e_meta'] = [ str(round(e_meta[s]['promoter'],3))+"%"+" / "+str(round(e_meta[s]['exon'],3)*100)+"%"+" / "+str(round(e_meta[s]['dhs'],3)*100)+"%" for s in samples ]
    #     json_dict["table"]['emeta_orig'] = e_meta
    print(json_dict)
    json_dump(json_dict, os.path.basename(d).replace("dataset", '') + ".json")
    f = open(os.path.basename(d).replace("dataset", '') + ".txt", 'w')
    print >>f, "SampleIDs\tFastQC\tRawReadsNumber\tUniquelyMappedRatio\tPBC\tFRiP\tUnionDHSRatio\tAllPeaksNumber\tPeaksFoldChange>=10\tPeaksFoldChange>=20\tPeaksPromoterPercentage\tPeaksExonPercentage\tPeaksIntronPercentage\tPeaksIntergenicPercentage"
    print >>f, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" %(os.path.basename(d).replace("dataset", ''), ','.join(map(str, json_dict['table']['fastqc'])), ','.join(map(str, json_dict['table']['map_number'])), ','.join(json_dict['table']['map']), ','.join(json_dict['table']['pbc']), ','.join(json_dict['table']['frip']), json_dict['table']['dhs'], '\t'.join(map(str, json_dict['table']['peaks'])), json_dict['table']['meta'][0].replace('/','\t'))
    f.close()

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print >>sys.stderr, "input dc directory"
        sys.exit(1)
    # env = Environment(loader=FileSystemLoader("."), trim_blocks=True, lstrip_blocks=True)
    # template = env.get_template("child.html")
    summary_to_json(sys.argv[1])


import json
import re
import os
import math
from samflow.command import ShellCommand, PythonCommand
from samflow.workflow import attach_back

from pkg_resources import resource_filename

from chilin2.modules.config.helpers import JinjaTemplateCommand, template_dump, json_load, json_dump

def stat_conservation(workflow, conf):
    collect = attach_back(workflow,
                PythonCommand(
                    json_conservation,
                    input={"score": conf.prefix + "_conserv.txt"},
                    output={"json": conf.json_prefix + "_conserv.json"},
                    param={"atype": conf.get("basics", "factor", "TF"), "id": conf.id},
                    name = "conservation score"))
    collect.allow_dangling = True
    collect.allow_fail = True

    if conf.long:  ## cluster figures, obsolete, keep for compatible
        fig = attach_back(workflow,
                    PythonCommand(conservation_figures,
                                  input ={"conservationR": conf.prefix + "_conserv.R",
                                          "historical_conservation_cluster_text": resource_filename("chilin2.modules.dbaccessor", "Histone_centers.txt")},


                                  output = {"R": conf.prefix+"_conserv_cluster.R",
                                            "compare_pdf": conf.prefix + "_conserv_compare.pdf"},
                                  param = {"id": conf.id}))
        fig.allow_fail = True
        fig.allow_dangling = True

def json_conservation(input={"score": ""}, output={"json": ""}, param={}):
    """
    collect conservation_plot output Phastcon score
    """
    json_dict = {"stat": [], "input": input, "output": output, "param": ""}
    rd = lambda x: str(round(float(x), 3))
    json_dict['stat'] = map(rd, open(input['score']).read().strip().split())
    json_dump(json_dict)

## remove obsolete clustering later
def _euclidian_distance(x, y):
   assert len(x) == len(y)
   sum_of_square = 0
   for x_i, y_i in zip(x, y):
       sum_of_square += pow((x_i - y_i), 2)
   distance = round(math.sqrt(sum_of_square), 4)
   return distance

def conservation_figures(input={"conservationR": "", "historical_conservation_cluster_text": ""},
                     output={"R": "", "compare_pdf": "", "pdf": ""},
                     param={"id": ""}):
   """
   For TFcenters data 1,2,3 pass, 4,5,6 fail
   For Histone center data 1,2,3,4 pass, 5,6,7,8 fail.
   """
   with open(input["conservationR"]) as conservation_r_file:
       for line in conservation_r_file:
           if re.findall(r'y0<-\S*\)', line):
               # TODO: Use more friendly regex
               value = re.findall(r'y0<-\S*\)', line)[0][6:-1]
               value = value.split(',')
           elif re.findall(r'x<-c\S*\)', line):
               xlab = line
       value = [float(i) for i in value]
       sumvalue = sum(value)
       value = [i / sumvalue for i in value]

   with open(input["historical_conservation_cluster_text"]) as historic_data:
       historyData = historic_data.readlines()

   less_than_this_id_means_pass_qc = len(historyData) / 2

   scoreList = []
   for i in range(len(historyData)):
       temp = historyData[i].strip()
       line = temp.split(' ')
       line = [float(j) for j in line]
       score = _euclidian_distance(value, line)
       scoreList.append(score)
   mindist = scoreList.index(min(scoreList)) # return the index of minimum distance group

   judgevalue = historyData[mindist].strip().split(' ')
   judgevalue = [str(i) for i in judgevalue]
   value = [str(i) for i in value]
   ymax = max(value + judgevalue)
   ymin = min(value + judgevalue)
   with open(output["R"], 'w') as f:
       f.write("pdf('%s',height=8.5,width=8.5)\n" % output["compare_pdf"])
       f.write("%s\n" % xlab)
       f.write("y1<-c(%s)\n" % ','.join(judgevalue))
       f.write("y2<-c(%s)\n" % ','.join(value))
       f.write("ymax<-(%s)\n" % ymax)
       f.write("ymin<-(%s)\n" % ymin)
       f.write("yquart <- (ymax-ymin)/4\n")
       f.write(
           "plot(x,y2,type='l',col='red',main='Normalized Conservation_at_summits',xlab='Distance from the Center (bp)',ylab='Normalized Average Phastcons',ylim=c(ymin-yquart,ymax+yquart))\n")
       f.write("points(x,y1,type='l',col='blue',lty=2)\n")
       f.write("legend('topleft',c('original group','compared group'),lty=c(1,2),col=c('red', 'blue'))\n")
       f.write("dev.off()\n")
   os.system('Rscript %s' % output["R"])

   # if mindist <= less_than_this_id_means_pass_qc:
   #     judge = 'Pass'
   # else:
   #     judge = 'Fail'

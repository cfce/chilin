import os
from chilin2.modules.config.helpers import JinjaTemplateCommand, template_dump, r_exec, json_load, json_dump, underline_to_space
## parse macs2 model R script
def get_size(rscript):
    values = {}
    with open(rscript) as model:
        for line in model:
            if line.startswith("p <- "):
                values["positive"] = line
            if line.startswith("m <- "):
                values["minus"] = line
            if line.startswith("x <- "):
                values["x"] = line
            if line.startswith("xcorr"):
                values['xcorr'] = line
            if line.startswith("ycorr"):
                values['ycorr'] = line
    return values

## calculate fragment mean length and standard variation
def stat_frag_std(input = {"r": "", "insert": ""}, output = {"json": "", "r": ""}, param = {"samples": "", "frag_tool": ""}):
    """ parse macs2 predictd r file into json file
    """
    json_dict = {"input": input, "output": output, "param": param, "stat": {}}
    for rin, rout, s in zip(input["r"], output["r"], param["samples"]):
        values = get_size(rin)
        with open(rout, 'w') as f:
            f.write(values['positive'])
            f.write(values['minus'])
            f.write(values['xcorr'])
            f.write(values['ycorr'])
            f.write("xcorr.max = xcorr[which(ycorr==max(ycorr))]\n")
            f.write(values['x'])
            f.write("p.expect = sum(x * p/100) \n")
            f.write("m.expect = sum(x * m/100) \n")
            f.write("p.sd = sqrt(sum(((x-p.expect)^2)*p/100)) \n")
            f.write("m.sd = sqrt(sum(((x-m.expect)^2)*m/100)) \n")
            f.write("cat(paste((p.sd + m.sd)/2, '\t', xcorr.max)) \n")
        f.close()
        std_frag = os.popen("Rscript %s" % rout).read().strip().split()
        json_dict["stat"][s] = "%s" % (int(float(std_frag[1])))
    json_dump(json_dict)

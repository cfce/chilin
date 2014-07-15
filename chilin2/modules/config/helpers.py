"""
separate latex template to several object
"""

import json
import random
from jinja2 import Environment, FileSystemLoader
from samflow.command import AbstractCommand, ShellCommand, PythonCommand
from pkg_resources import resource_filename

env = Environment(loader = FileSystemLoader("/"),
    block_start_string = '\BLOCK{',
    block_end_string = '}',
    variable_start_string = '\VAR{',
    variable_end_string = '}',
    comment_start_string = '\#{',
    comment_end_string = '}',
    line_statement_prefix = '%-',
    line_comment_prefix = '%#',
    trim_blocks = True,
    autoescape = False,)

def surround_by_quote(a_list):
    return ['"%s"' % an_element for an_element in a_list]

def count_in_million(x):
    if type(x) == int:
        return str(round(x/1000000., 1)) + "M"

def decimal_to_latex_percent(dec):
    if type(dec) == float:
        return str(round(dec*100, 2))    + "\%"

def underline_to_space(x):
    if type(x) == str:
        return x.replace("_", " ")
    return x

env.filters["surround_by_quote"] = surround_by_quote


class JinjaTemplateCommand(AbstractCommand):
    def __init__(self, template, tool=None, param = {}, input=[], output=[], name = ""):
        AbstractCommand.__init__(self, template=template, tool=None, param = param, input=[], output=[], name = name)
        self.env = env

        self._t = self.env.get_template(self.template)

    def _execute(self):
        """ Real-run current command"""
        self.result = self._t.render(input = self.input, output = self.output, **self.param)
        return True

    def _simulate(self):
        """ Dry-run current command: Pretend to run but not invoke anything """
        print("Rendering Latex part %s" % self.name, self.template)
        return True

def template_dump(jinja_template):
    jinja_template.invoke()
    with open(jinja_template.param["render_dump"], "w") as f:
        f.write(jinja_template.result)

def r_exec(jinja_template_r):
    ShellCommand(template="Rscript {input}",
        input=jinja_template_r.param["render_dump"],
        output=jinja_template_r.param["pdf"]).invoke()

def json_dump(json_dict):   # json
    """
    dump out uniform json files for collecting statistics
    :param json_dict: output python dict to json
    :return: json_file name
    """
    json_file = json_dict["output"]["json"]
    with open(json_file, "w") as f:
        json.dump(json_dict, f, indent=4)
    return json_file


def json_load(json_file):
    with open(json_file, "r") as f:
        json_dict = json.load(f)
    return json_dict


def latex_start(input = {"template": ""}, output = {"latex": ""}, param = {"id": ""}):
    end_latex = JinjaTemplateCommand(
        name = "end of latex document",
        template = input["template"],
        param = {"section_name": "begin",
                 "prefix_dataset_id": underline_to_space(param["id"]),
                 "logo": resource_filename("chilin2.modules", "summary/CFCE_Logo_Final.jpg"),
                 "render_dump": output["latex"],
                 "bmcard": param["bmcard"],
                 "version": param["version"],
                 "user": param["user"],
        })

    template_dump(end_latex)


def latex_end(input={"template": ""}, output = {"latex": ""}, param = {}):
    end_latex = JinjaTemplateCommand(
        name="end of latex document",
        template=input["template"],
        param={"section_name": "ending",
                 "render_dump": output["latex"]})

    template_dump(end_latex)


def make_link_command(orig, dest):
    """
    link original input to destination files
    :param orig: input
    :param dest: link symbol
    :return: ShellCommand Class
    ln has machine type problem
    """ ## not use symbol link
    return ShellCommand("cp -fr {input} {output}",
                        input=orig,
                        output=dest,
                        name="copy")


def sampling(orig, dest, rand, format, conf): # call fastq_sampling
    """
    prepare sampling fastq files for library contamination and fastqc
    rand: the number of random selected fastq reads
    use lh3's https://github.com/lh3/seqtk/ to sample fastq and fastq.gz
    """
    if format == "fastq":
        #return PythonCommand(fastq_sampling,
        #                     input=orig,
        #                     output=dest,
        #                     param={"random_number": rand})
        ## faster and support fastq.gz
        ## if paired end, we must use same -s
        return ShellCommand("{tool} sample -s 11 {input[fastq]} {param[rand]} > {output[fastq_sample]}",
                            tool = "seqtk",
                            input = orig,
                            output = dest,
                            param = {"rand": 100000})

    elif format == "sam":
        ## samtools sampling
        ## add judge condition
        return ShellCommand("""
                            count=$({tool} view -Sc {input[sam]})
                            ## judge mapped reads number less than sampling number
                            if [ $count -le {param[random_number]} ]
                            then
                                ln -f {input[sam]} {input[sam]}.{param[random_number]}
                                {tool} view -bS {input[sam]}.{param[random_number]} > {output[samp]}
                            else
                                sampling_pe_sam.py {input[sam]} {param[random_number]}
                                {tool} view -bS {input[sam]}.{param[random_number]} > {output[samp]}
                            fi
                            """,
                            tool = "samtools",
                            input={"sam": orig},
                            output={"samp": dest},
                            param={"random_number": rand},
                            name = "sampling bam")


def fastq_sampling(input, output, param):   ## a sampling helper
    """
    get N random headers from a fastq
    file without reading the whole thing
    into memory
    """
    num_lines = int(sum(1 for _ in open(input["fastq"])) / 4)

    rand_nums = sorted([random.randint(0, num_lines - 1) for _ in range(param["random_number"])])

    fastq = open(input["fastq"],"rU")
    fastq_sample = open(output["fastq_sample"], "w")

    cur_num = - 1
    for rand_num in rand_nums:
        while cur_num < rand_num:
            for i in range(4):
                fastq.readline()
            cur_num += 1

        for i in range(4):
            fastq_sample.write(fastq.readline())
        cur_num += 1

    fastq_sample.close()
    fastq.close()

    print("wrote to %s" % (output["fastq_sample"]))

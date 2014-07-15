
def _browser_hub(workflow, conf):   # link to washU epigenomics washU
    """
    write out json files for linking to epigenomics brower,
    this should include encodePeak, bed and bigWig files annotation
    :param workflow: samflow defined class
    :param conf: parsed config files
    :return: void
    """
    # hub = attach_back(workflow,
    #                   PythonCommand(
    #                       ,
    #                   ))

#
# def _seqpos_latex(workflow, conf):
#     attach_back(workflow,
#         PythonCommand(
#             latex_seqpos,
#             input={"json": conf.json_prefix + "_seqpos.json",
#                    "template": latex_template},
#             output={"latex": conf.latex_prefix + "_seqpos.latex"}))

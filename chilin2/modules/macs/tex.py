
def latex_macs2(input, output, param):
    # TODO: qian work out the peaks_summary_result part
    json_dict = json_load(input["json"])

    summary = [underline_to_space(json_dict["param"]["id"]),
               json_dict["stat"]["qvalue"],
               json_dict["stat"]["totalpeak"],
               json_dict["stat"]["peaksge10"],
               json_dict["stat"]["shiftsize"]]

    high_confident_latex = JinjaTemplateCommand(
        name = "high confident latex",
        template = input["template"],
        param = {"section_name": "high_confident_peaks",
                 "peak_summary_table": summary,
                 "high_confident_peak_graph": json_dict["output"]["pdf"],
                 "render_dump": output["latex"]})

    template_dump(high_confident_latex)


def latex_macs2_on_sample(input, output, param):
    json_dict = json_load(input["json"])
    latex = JinjaTemplateCommand(
        name="redunRateQC",
        template=input["template"],
        param = {"section_name": "redundant",
                 "redundant_ratio_graph": json_dict["output"]["pdf"],
                 "render_dump": output["latex"]})

    template_dump(latex)


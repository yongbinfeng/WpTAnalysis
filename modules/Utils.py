"""
collecting some util scripts here
"""

import pandas as pd

def FormatOutputForWZ(istring: str):
    labelmaps = {}
    labelmaps['muplus'] = '$\\mathrm{W}^{+}\\rightarrow\\mu^{+}\\nu$'
    labelmaps['muminus'] = '$\\mathrm{W}^{-}\\rightarrow\\mu^{-}\\bar{\\nu}$'
    labelmaps['mumu'] = '$\\mathrm{Z}\\rightarrow\\mu^{+}\\mu^{-}$'
    labelmaps['ee'] = '$\\mathrm{Z}\\rightarrow e^{+}e^{-}$'
    labelmaps['eplus'] = '$\\mathrm{W}^{+}\\rightarrow e^{+}\\nu$'
    labelmaps['eminus'] = '$\\mathrm{W}^{-}\\rightarrow e^{-}\\bar{\\nu}$'

    procmaps = {}
    procmaps['data'] = 'Data'
    procmaps['ewk'] = "EWK"
    procmaps['qcd'] = "QCD"
    procmaps['ttbar'] = "$t\\bar{t}$"

    for key in labelmaps.keys():
        istring = istring.replace(key, labelmaps[key])
    
    for key in procmaps.keys():
        istring = istring.replace(key, procmaps[key])

    return istring


def FormatTable(pdict: str, columns: list = None, caption: str = None, label: str = None):
    """
    given a dictionary, print the latex version of the table
    """
    df = pd.DataFrame(pdict, columns=columns)
    output = df.to_latex(float_format="{:.1f}".format, caption = caption, label = label)
    output = output.replace('\\toprule', '\\hline').replace('\\midrule', '\\hline').replace('\\bottomrule','\\hline')

    output = FormatOutputForWZ(output)

    return output
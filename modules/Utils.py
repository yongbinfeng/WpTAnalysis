"""
collecting some util scripts here
"""

import pandas as pd
from collections import OrderedDict
import re


def GetValues(pdict: OrderedDict, key: str):
    """
    given a dictionary, return all the values with keys matching the given key
    """
    results = []
    for x, val in pdict.items():
        if re.match(key, x):
            if type(val) is not list:
                val = [val]
            results += val
    return results

def FormatROOTInput(istring: str):
    labelmaps = {}
    
    labelmaps = {}
    labelmaps['lepplus']  = 'W^{+}#rightarrow l^{+}#nu'
    labelmaps['lepminus'] = 'W^{-}#rightarrow l^{-}#bar{#nu}'
    labelmaps['leplep']   = 'Z#rightarrow l^{+} l^{-}'
    labelmaps['lepinc']   = 'W^{#pm}#rightarrow l^{#pm}#nu'
    
    labelmaps["WOverZ"] = "W^{#pm}/Z"
    labelmaps["WchgRatio"] = "W^{+}/W^{-}"
    labelmaps['WpOverWm'] = labelmaps['WchgRatio']
    labelmaps['WZRatio'] = labelmaps['WOverZ']
    
    labelmaps['Wplus'] = labelmaps['lepplus']
    labelmaps['Wminus'] = labelmaps['lepminus']
    labelmaps['Zinc'] = labelmaps['leplep']
    
    for key in labelmaps.keys():
        istring = istring.replace(key, labelmaps[key])
        
    return istring
    

def FormatOutputForWZ(istring: str):
    labelmaps = {}
    labelmaps['muplus'] = '$\\mathrm{W}^{+}\\rightarrow\\mu^{+}\\nu$'
    labelmaps['muminus'] = '$\\mathrm{W}^{-}\\rightarrow\\mu^{-}\\bar{\\nu}$'
    labelmaps['mumu'] = '$\\mathrm{Z}\\rightarrow\\mu^{+}\\mu^{-}$'
    labelmaps['ee'] = '$\\mathrm{Z}\\rightarrow e^{+}e^{-}$'
    labelmaps['eplus'] = '$\\mathrm{W}^{+}\\rightarrow e^{+}\\nu$'
    labelmaps['eminus'] = '$\\mathrm{W}^{-}\\rightarrow e^{-}\\bar{\\nu}$'
    labelmaps['lepplus'] = '$\\mathrm{W}^{+}\\rightarrow \\ell^{+}\\nu$'
    labelmaps['lepminus'] = '$\\mathrm{W}^{-}\\rightarrow \\ell^{-}\\bar{\\nu}$'
    labelmaps['leplep'] = '$\\mathrm{Z}\\rightarrow \\ell^{+}\\ell^{-}$'
    labelmaps['Winc'] = '$\\mathrm{W}^{\\pm}\\rightarrow \\ell^{\\pm}\\nu$'
    labelmaps['WOverZ'] = '$\\mathrm{W}^{\pm}/\\mathrm{Z}$'
    labelmaps['WpOverWm'] = '$\\mathrm{W}^{+}/\\mathrm{W}^{-}$'
    labelmaps['WZRatio'] = '$\\mathrm{W}^{\pm}/\\mathrm{Z}$'
    labelmaps['WchgRatio'] = '$\\mathrm{W}^{+}/\\mathrm{W}^{-}$'
    
    labelmaps['Wplus'] = labelmaps['lepplus']
    labelmaps['Wminus'] = labelmaps['lepminus']
    labelmaps['Zinc'] = labelmaps['leplep']

    procmaps = {}
    procmaps['data'] = 'Data'
    procmaps['sig'] = "Signal"
    procmaps['ewk'] = "EWK"
    procmaps['qcd'] = "QCD"
    procmaps['ttbar'] = "$t\\bar{t}$"

    sysmaps = {}
    sysmaps['lumi'] = 'Lumi'
    sysmaps['recoil'] = 'Recoil'
    sysmaps['QCDbkg'] = 'Bkg QCD'
    sysmaps['effstat'] = 'Efficiency Stat.'
    sysmaps['prefire'] = 'Prefire'
    sysmaps['QCDscale'] = 'QCD Scale'

    sysmaps['effsys'] = 'Efficiency Syst.'
    sysmaps['pdfalphaS'] = 'PDF + $\\alpha_\\mathrm{S}$'
    sysmaps['mcsec'] = 'MC Norm'
    
    sysmaps['resummFSR'] = "Resum. + FSR"
    
    xsecNames = {}
    xsecNames['xsec'] = '$\\sigma$'
    
    hacks = {}
    hacks['lepplus\_diff'] = 'Diff $[\%]$'
    hacks['lepminus\_diff'] = 'Diff $[\%]$'
    hacks['leplep\_diff'] = 'Diff $[\%]$'
    hacks['lepinc\_diff'] = 'Diff $[\%]$'
    hacks['WchgRatio\_diff'] = 'Diff $[\%]$'
    hacks['WZRatio\_diff'] = 'Diff $[\%]$'
    
    for key in hacks.keys():
        istring = istring.replace(key, hacks[key])
        
    for key in labelmaps.keys():
        istring = istring.replace(key, labelmaps[key])
    
    for key in procmaps.keys():
        istring = istring.replace(key, procmaps[key])

    for key in sysmaps.keys():
        istring = istring.replace(key, sysmaps[key])
        
    for key in xsecNames.keys():
        istring = istring.replace(key, xsecNames[key])

    return istring


def FormatTable(pdict: str, columns: list = None, caption: str = None, label: str = None, precision: int=1, escape: bool = True):
    """
    given a dictionary, print the latex version of the table
    """
    df = pd.DataFrame(pdict, columns=columns)
    df = df.round(precision)
    #output = df.to_latex(float_format="{:.1f}".format, caption = caption, label = label)
    #output = df.to_latex(caption = caption, label = label)
    output = df.to_latex(escape = escape)
    output = output.replace('\\toprule', '\\hline').replace('\\midrule', '\\hline').replace('\\bottomrule','\\hline')

    output = FormatOutputForWZ(output)

    return output

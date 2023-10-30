"""
collecting some util scripts here
"""

import pandas as pd
from collections import OrderedDict
import re
import copy
import math

def findPrecision(values):
    err = max(values[1:])
    if err!= 0:
        precision = int(math.log10(float(f"{err:.2g}")))-1
    else:
        precision = 0
    isInt = False
    if precision < 0 and err < 1.0:
        # err < 1.0
        precision = precision -1
    elif precision >= 0:
        # err > 10.0
        isInt = True 
    return precision, isInt

def roundToError(values, precision = None, isInt = None):
    if precision == None:
        precision, isInt = findPrecision(values)
    values_new = []
    for val in values:
        if not isInt:
            values_new.append(f"{round(val, -precision):.{-precision}f}")
        else:
            #values_new.append(f"{round(val, -precision):0f}")
            values_new.append(str(int(round(val, -precision))))
    return tuple(values_new)


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
    
    labelmaps["WOverZ"] = "W^{#pm}#rightarrow l^{#pm}#nu / Z#rightarrow l^{+} l^{-}"
    labelmaps["WchgRatio"] = "W^{+}#rightarrow l^{+}#nu / W^{-}#rightarrow l^{-}#bar{#nu}"
    labelmaps['WpOverWm'] = labelmaps['WchgRatio']
    labelmaps['WZRatio'] = labelmaps['WOverZ']
    
    labelmaps['winc'] = labelmaps['lepinc']
    labelmaps['Winc'] = labelmaps['lepinc']
    labelmaps['Wplus'] = labelmaps['lepplus']
    labelmaps['Wminus'] = labelmaps['lepminus']
    labelmaps['Zinc'] = labelmaps['leplep']
    
    for key in labelmaps.keys():
        istring = istring.replace(key, labelmaps[key])
        
    return istring
    

def FormatOutputForWZ(istring: str, isYield: bool = False):
    labelmaps = {}
    labelmaps['muplus']    = '$\\PW^{+}\\rightarrow\\PGm^{+}\\PGn$'
    labelmaps['muminus']   = '$\\PW^{-}\\rightarrow\\PGm^{-}\\PAGn$'
    labelmaps['mumu']      = '$\\PZ\\rightarrow\\PGm^{+}\\PGm^{-}$'
    labelmaps['ee']        = '$\\PZ\\rightarrow\\Pe^{+}\\Pe^{-}$'
    labelmaps['eplus']     = '$\\PW^{+}\\rightarrow\\Pe^{+}\\PGn$'
    labelmaps['eminus']    = '$\\PW^{-}\\rightarrow\\Pe^{-}\\PAGn$'
    labelmaps['lepplus']   = '$\\PW^{+}\\rightarrow \\Pell^{+}\\PGn$'
    labelmaps['lepminus']  = '$\\PW^{-}\\rightarrow \\Pell^{-}\\PAGn$'
    labelmaps['leplep']    = '$\\PZ\\rightarrow \\Pell^{+}\\Pell^{-}$'
    labelmaps['Winc']      = '$\\PW^{\\pm}\\rightarrow \\Pell^{\\pm}\\PGn$'
    labelmaps['WOverZ']    = '$\\PW^{\pm}/\\PZ$'
    labelmaps['WpOverWm']  = '$\\PW^{+}/\\PW^{-}$'
    labelmaps['WZRatio']   = '$\\PW^{\pm}/\\PZ$'
    labelmaps['WchgRatio'] = '$\\PW^{+}/\\PW^{-}$'
    
    labelmaps['Wplus'] = labelmaps['lepplus']
    labelmaps['Wminus'] = labelmaps['lepminus']
    labelmaps['Zinc'] = labelmaps['leplep']
    labelmaps['winc'] = labelmaps['Winc']

    procmaps = {}
    procmaps['Measured'] = 'Data'
    procmaps['data'] = 'Data'
    procmaps['sig'] = "Signal"
    procmaps['ewk'] = "Electroweak"
    procmaps['qcd'] = "QCD multijet"
    procmaps['ttbar '] = "$\\ttbar$"

    sysmaps = {}
    sysmaps['lumi'] = 'Lumi'
    sysmaps['recoil'] = 'Hadronic recoil calibration'
    sysmaps['QCDbkg'] = 'Bkg QCD'
    sysmaps['effstat'] = 'Efficiency (stat)'
    sysmaps['prefire'] = 'Trigger prefire correction'
    sysmaps['qcdscale'] = '$\\mu_{\\mathrm{R}}$ and $\\mu_{\\mathrm{F}}$ scales'
    sysmaps['effsys'] = 'Efficiency (syst)'
    sysmaps['pdfalphaS'] = 'PDF + $\\alpha_\\mathrm{S}$'
    sysmaps['mcsec'] = 'EWK+$\\ttbar$ cross section'
    sysmaps['qcdsys'] = 'QCD multijet (syst)'
    sysmaps['qcdstat'] = 'QCD multijet (stat)'
    sysmaps['binByBinStat'] = 'MC sim. stat'
    sysmaps['resummFSR'] = "Resum. + FSR"
    sysmaps['QCDsys'] = sysmaps['qcdsys']
    sysmaps['QCDstat'] = sysmaps['qcdstat']
    sysmaps['QCDscale'] = sysmaps['qcdscale']
    
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
        # multiple columns for yield tables
        val = labelmaps[key]
        if isYield:
            val = "\\multicolumn{2}{c}{" + val + "}"
        istring = istring.replace(key, val)
    
    for key in sysmaps.keys():
        istring = istring.replace(key, sysmaps[key])
        
    for key in procmaps.keys():
        istring = istring.replace(key, procmaps[key])
        
    for key in xsecNames.keys():
        istring = istring.replace(key, xsecNames[key])

    return istring


def FormatTable(pdict: str, columns: list = None, precision: int=1, escape: bool = False, doTranspose = False, isYield : bool = False):
    """
    given a dictionary, print the latex version of the table
    """
    results = copy.deepcopy(pdict)
    for bkey, bval in results.items():
        print("bkey ", bkey)
        if isYield:
            # for event yield, no need to sync precision for different rows
            for key, val in bval.items():
                if key == 'data':
                    rval = roundToError(val, 0, True)
                    val = f"{rval[0]}"
                    val = "\\multicolumn{2}{c}{" + val + "}"
                else:
                    vprecision, visInt = findPrecision(val) 
                    visInt = True
                    rval = roundToError(val, vprecision, visInt)
                    val = f"{rval[0]} & {rval[1]}"
                bval[key] = val
        else:
            # loop over results first; find the lowest precision
            # such that all the values can be rounded to the same precision for one column
            vprecision = -100
            visInt = False
            for key, val in bval.items():
                if type(val) is tuple:
                    vtmp_precision, vtmp_isInt = findPrecision(val)
                    vprecision = max(vprecision, vtmp_precision)
                    visInt = visInt or vtmp_isInt
            for key, val in bval.items():
                print("key ", key)
                if type(val) is tuple:
                    rval = roundToError(val, vprecision, visInt)
                    if key == 'Measured':
                        if not "Ratio" in bkey and "Over" not in bkey:
                            val = f"{rval[0]}\pm{rval[1]}_\mathrm{{stat}}\pm{rval[2]}_\mathrm{{syst}}\pm{rval[3]}_\mathrm{{lum}}"
                        else:
                            val = f"{rval[0]}\pm{rval[1]}_\mathrm{{stat}}\pm{rval[2]}_\mathrm{{syst}}"
                    else:
                        # theory predictions
                        val = f"{rval[0]}^{{+{rval[2]}}}_{{-{rval[1]}}}"
                    val = "$" + val +"$"
                    bval[key] = val
        print("bval ", bval)
        results[bkey] = bval
        
    pd.set_option('display.max_colwidth', None)
    df = pd.DataFrame(results, columns=columns)
    if doTranspose:
        df = df.transpose()
    df = df.round(precision)
    #output = df.to_latex(float_format="{:.1f}".format, caption = caption, label = label)
    #output = df.to_latex(caption = caption, label = label)
    output = df.to_latex(escape = escape)
    output = output.replace('\\toprule', '\\hline').replace('\\midrule', '\\hline').replace('\\bottomrule','\\hline').replace('\\textbackslash pm', '\\pm').replace("\$", "$")
    
    print("OUtput before Format: ", output)

    output = FormatOutputForWZ(output, isYield)

    return output

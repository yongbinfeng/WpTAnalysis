import json
import ROOT
from collections import OrderedDict
import math

f_xsecs_MC = "data/xsecs.json"

with open(f_xsecs_MC) as f:
    xsecs_MC = json.load(f)
    
print(xsecs_MC)

maps = {
    "lepplus": "wmp",
    "lepminus": "wmm",
    "leplep": "zmm",
}

lumiUncs = {
    "13TeV": 0.017,
    "5TeV": 0.019,
    "sqrtS": math.sqrt(0.015**2 + 0.017**2)
}

f_fit = "forCombine_WithVpT/Fiducial/test0/commands/scripts/card_combined.root"

f_root = ROOT.TFile(f_fit)
tree = f_root.Get("fitresults")
tree.GetEntry(0)
print(tree)

xsecs = OrderedDict()
xsecs_ratios = OrderedDict()
xsecs_ratios_sqrtS = OrderedDict()
xsecs_ratios_sqrtS_ratios = OrderedDict()

def getRelError(bname):
    return getattr(tree, f"{bname}_err") / getattr(tree, bname)

for sqrtS in ["13TeV", "5TeV"]:
    xsecs[sqrtS] = OrderedDict()
    xsecs_ratios[sqrtS] = OrderedDict()
    
    for ch in ["lepplus", "lepminus", "leplep"]:
        "lepplus_13TeV_sig_mu"
        mu = getattr(tree, f"{ch}_{sqrtS}_sig_mu")
        mu_err = getattr(tree, f"{ch}_{sqrtS}_sig_mu_err")
        print("mu value for", ch, sqrtS, " is ", mu, " +/- ", mu_err)
        
        xsec_MC = xsecs_MC[sqrtS]['born'][maps[ch]]
        print("xsec value for", ch, sqrtS, "is", xsec_MC)
        
        xsecs[sqrtS][ch] = (mu * xsec_MC, mu_err * xsec_MC, lumiUncs[sqrtS] * xsec_MC)
        print("xsec fit value for", ch, sqrtS, "is", xsecs[sqrtS][ch])
        
    mu_err = getRelError(f"Winc_{sqrtS}_sig_sumxsec")
    xsec_tot_W = xsecs[sqrtS]["lepplus"][0] + xsecs[sqrtS]["lepminus"][0]
    xsecs[sqrtS]["Winc"] = (xsec_tot_W, mu_err * xsec_tot_W)
        
    # cross section ratios
    # W+ / W-
    mu_err = getRelError(f"WchgRatio_{sqrtS}_ratiometaratio")
    xsec_ratio_Wchg = xsecs[sqrtS]["lepplus"][0] / xsecs[sqrtS]["lepminus"][0]
    xsecs_ratios[sqrtS]["WchgRatio"] = (xsec_ratio_Wchg, mu_err * xsec_ratio_Wchg, 0.)
    # W / Z
    mu_err = getRelError(f"WZRatio_{sqrtS}_ratiometaratio")
    xsec_ratio_WZ = xsecs[sqrtS]["Winc"][0] / xsecs[sqrtS]["leplep"][0]
    xsecs_ratios[sqrtS]["WZRatio"] = (xsec_ratio_WZ, mu_err * xsec_ratio_WZ, 0.)
    
# measure the ratios between sqrtS
mu_err = getRelError("sqrtS_Wplus_ratio_ratiometaratio")
xsec_ratio_sqrtS_Wplus = xsecs["13TeV"]["lepplus"][0] / xsecs["5TeV"]["lepplus"][0]
xsecs_ratios_sqrtS["lepplus"] = (xsec_ratio_sqrtS_Wplus, mu_err * xsec_ratio_sqrtS_Wplus, lumiUncs["sqrtS"] * xsec_ratio_sqrtS_Wplus)

mu_err = getRelError("sqrtS_Wminus_ratio_ratiometaratio")
xsec_ratio_sqrtS_Wminus = xsecs["13TeV"]["lepminus"][0] / xsecs["5TeV"]["lepminus"][0]
xsecs_ratios_sqrtS["lepminus"] = (xsec_ratio_sqrtS_Wminus, mu_err * xsec_ratio_sqrtS_Wminus, lumiUncs["sqrtS"] * xsec_ratio_sqrtS_Wminus)

mu_err = getRelError("sqrtS_Zinc_ratio_ratiometaratio")
xsec_ratio_sqrtS_Z = xsecs["13TeV"]["leplep"][0] / xsecs["5TeV"]["leplep"][0]
xsecs_ratios_sqrtS["leplep"] = (xsec_ratio_sqrtS_Z, mu_err * xsec_ratio_sqrtS_Z, lumiUncs["sqrtS"] * xsec_ratio_sqrtS_Z)

mu_err = getRelError("sqrtS_Winc_ratio_ratiometaratio")
xsec_ratio_sqrtS_Winc = xsecs["13TeV"]["Winc"][0] / xsecs["5TeV"]["Winc"][0]
xsecs_ratios_sqrtS["lepinc"] = (xsec_ratio_sqrtS_Winc, mu_err * xsec_ratio_sqrtS_Winc, lumiUncs["sqrtS"] * xsec_ratio_sqrtS_Winc)

# double ratios between sqrtS
#mu_err = getRelError("sqrtS_WchgRatio_ratio_doubleratiometatotalxsec")
mu_err = getRelError("sqrtS_WchgRatio_ratio_doubleratiometaratio")
xsec_ratio_sqrtS_Wchg = xsecs_ratios["13TeV"]["WchgRatio"][0] / xsecs_ratios["5TeV"]["WchgRatio"][0]
xsecs_ratios_sqrtS_ratios["WchgRatio"] = (xsec_ratio_sqrtS_Wchg, mu_err * xsec_ratio_sqrtS_Wchg, 0.)

#mu_err = getRelError("sqrtS_WZRatio_ratio_doubleratiometatotalxsec")
mu_err = getRelError("sqrtS_WZRatio_ratio_doubleratiometaratio")
xsec_ratio_sqrtS_WZ = xsecs_ratios["13TeV"]["WZRatio"][0] / xsecs_ratios["5TeV"]["WZRatio"][0]
xsecs_ratios_sqrtS_ratios["WZRatio"] = (xsec_ratio_sqrtS_WZ, mu_err * xsec_ratio_sqrtS_WZ, 0.)

    
print("results xsecs: ", xsecs)
print("results xsecs ratios: ", xsecs_ratios)
print("results xs ratios sqrtS: ", xsecs_ratios_sqrtS)
print("results xs ratios sqrtS ratios: ", xsecs_ratios_sqrtS_ratios)

pdfsets = ["nnpdf4.0", "nnpdf3.1", "ct14", "msht20"]
    
# compare the theory predictions with the measured values
from theoryResults import *        
results = OrderedDict()
for sqrtS in ["13TeV", "5TeV"]:
    results[sqrtS] = OrderedDict()
    
    for ch in ["lepplus", "lepminus", "leplep"]:
        results[sqrtS][ch] = OrderedDict()
        
        results[sqrtS][ch]['Measured'] = xsecs[sqrtS][ch][0]
        for pdf in pdfsets:
            results[sqrtS][ch][pdf.upper()] = xsec_fid[pdf].GetXsec(ch, sqrtS)
            
        results[sqrtS][f"{ch}_diff"] = OrderedDict()
        results[sqrtS][f"{ch}_diff"]['Measured'] = math.sqrt((xsecs[sqrtS][ch][1]/xsecs[sqrtS][ch][0])**2 + (xsecs[sqrtS][ch][2]/xsecs[sqrtS][ch][0])**2) * 100.0
        for pdf in pdfsets:
            results[sqrtS][f"{ch}_diff"][pdf.upper()] = (xsec_fid[pdf].GetXsec(ch, sqrtS) / xsecs[sqrtS][ch][0] - 1.) * 100.0
    
    # measure the V/V ratios at the same sqrtS
    results[f"{sqrtS}_ratios"] = OrderedDict()
    for ch in ["WchgRatio", "WZRatio"]:
        results[f"{sqrtS}_ratios"][ch] = OrderedDict()
        
        results[f"{sqrtS}_ratios"][ch]['Measured'] = xsecs_ratios[sqrtS][ch][0]
        for pdf in pdfsets:
            results[f"{sqrtS}_ratios"][ch][pdf.upper()] = getattr(xsec_fid[pdf], f"Get{ch}")(sqrtS == "13TeV")
            
        results[f"{sqrtS}_ratios"][f"{ch}_diff"] = OrderedDict()
        results[f"{sqrtS}_ratios"][f"{ch}_diff"]['Measured'] = math.sqrt((xsecs_ratios[sqrtS][ch][1]/xsecs_ratios[sqrtS][ch][0])**2 + (xsecs_ratios[sqrtS][ch][2]/xsecs_ratios[sqrtS][ch][0])**2) * 100.0
        
        for pdf in pdfsets:
            results[f"{sqrtS}_ratios"][f"{ch}_diff"][pdf.upper()] = (getattr(xsec_fid[pdf], f"Get{ch}")(sqrtS == "13TeV") / xsecs_ratios[sqrtS][ch][0] - 1.) * 100.0
            
# measure the V xsec ratios between sqrtS
results["sqrtS_ratios"] = OrderedDict()
for ch in ["lepplus", "lepminus", "leplep"]:
    results["sqrtS_ratios"][ch] = OrderedDict()
    
    results['sqrtS_ratios'][ch]['Measured'] = xsecs_ratios_sqrtS[ch][0]
    
    for pdf in pdfsets:
        results['sqrtS_ratios'][ch][pdf.upper()] = xsec_fid[pdf].GetSqrtSRatio(ch)
    
    results['sqrtS_ratios'][f"{ch}_diff"] = OrderedDict()
    results['sqrtS_ratios'][f"{ch}_diff"]['Measured'] = math.sqrt((xsecs_ratios_sqrtS[ch][1]/xsecs_ratios_sqrtS[ch][0])**2 + (xsecs_ratios_sqrtS[ch][2]/xsecs_ratios_sqrtS[ch][0])**2) * 100.0
    
    for pdf in pdfsets:
        results['sqrtS_ratios'][f"{ch}_diff"][pdf.upper()] = (xsec_fid[pdf].GetSqrtSRatio(ch) / xsecs_ratios_sqrtS[ch][0] - 1.) * 100.0
        
# measure the V/V ratios between sqrtS
results["sqrtS_double_ratios"] = OrderedDict()
for ch in ["WchgRatio", "WZRatio"]:
    results["sqrtS_double_ratios"][ch] = OrderedDict()
    
    results['sqrtS_double_ratios'][ch]['Measured'] = xsecs_ratios_sqrtS_ratios[ch][0]
    
    for pdf in pdfsets:
        results['sqrtS_double_ratios'][ch][pdf.upper()] = xsec_fid[pdf].GetSqrtSDoubleRatio(ch)
    
    results['sqrtS_double_ratios'][f"{ch}_diff"] = OrderedDict()
    results['sqrtS_double_ratios'][f"{ch}_diff"]['Measured'] = math.sqrt((xsecs_ratios_sqrtS_ratios[ch][1]/xsecs_ratios_sqrtS_ratios[ch][0])**2 + (xsecs_ratios_sqrtS_ratios[ch][2]/xsecs_ratios_sqrtS_ratios[ch][0])**2) * 100.0
    
    for pdf in pdfsets:
        results['sqrtS_double_ratios'][f"{ch}_diff"][pdf.upper()] = (xsec_fid[pdf].GetSqrtSDoubleRatio(ch) / xsecs_ratios_sqrtS_ratios[ch][0] - 1.) * 100.0
        
        
outdir = "plots/xsecResults/"
    
from modules.Utils import FormatTable
from modules.postFitScripts import WriteOutputToText
outputs = FormatTable(results['13TeV'])
print(outputs)
WriteOutputToText(outputs, f"{outdir}/tables/results_13TeV_xsecs.tex") 

outputs = FormatTable(results['5TeV'])
print(outputs)
WriteOutputToText(outputs, f"{outdir}/tables/results_5TeV_xsecs.tex")

outputs = FormatTable(results['13TeV_ratios'], precision = 4)
print(outputs)
WriteOutputToText(outputs, f"{outdir}/tables/results_13TeV_ratios.tex")

outputs = FormatTable(results['5TeV_ratios'], precision = 4)
print(outputs)
WriteOutputToText(outputs, f"{outdir}/tables/results_5TeV_ratios.tex")

outputs = FormatTable(results['sqrtS_ratios'], precision = 4)
print(outputs)
WriteOutputToText(outputs, f"{outdir}/tables/results_sqrtS_ratios.tex")

outputs = FormatTable(results['sqrtS_double_ratios'], precision = 4)
print(outputs)
WriteOutputToText(outputs, f"{outdir}/tables/results_sqrtS_double_ratios.tex")

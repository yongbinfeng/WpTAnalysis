import json
import ROOT
from collections import OrderedDict
import math
from CMSPLOTS.myFunction import DrawHistos
from modules.postFitScripts import GetPOIValue
from modules.Utils import roundToError
from copy import deepcopy
import sys

doFiducial = True

f_xsecs_MC = "data/xsecs_fid.json" if doFiducial else "data/xsecs_inc.json"

with open(f_xsecs_MC) as f:
    xsecs_MC = json.load(f)
    
#print(xsecs_MC)

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

if doFiducial:
    f_fit = "forCombine_WithVpT/Fiducial/test0/commands/scripts/card_combined.root"
else:
    f_fit = "forCombine_WithVpT/Inclusive/test0/commands/scripts/card_combined.root"
f_fit_fid = "forCombine_WithVpT/Fiducial/test0/commands/scripts/card_combined.root"

f_root = ROOT.TFile(f_fit)
tree = f_root.Get("fitresults")
tree.GetEntry(0)

f_root_fid = ROOT.TFile(f_fit_fid)
tree_fid = f_root.Get("fitresults")
tree_fid.GetEntry(0)

xsecs = OrderedDict()
xsecs_ratios = OrderedDict()
xsecs_ratios_sqrtS = OrderedDict()
xsecs_ratios_sqrtS_ratios = OrderedDict()

def getRelError(ifilename, poiname):
    val_poi, err_poi = GetPOIValue(ifilename, poiname)
    return err_poi / val_poi

def getRelStatSystError(ifilename, poiname, hname):
    ifile = ROOT.TFile(ifilename)
    himpact_grouped = ifile.Get(hname)

    # find the POI bin for poiname
    ibinX = -1
    for binX in range(1, himpact_grouped.GetNbinsX()+1):
        poi = himpact_grouped.GetXaxis().GetBinLabel(binX)
        if poi == poiname:
            ibinX = binX
            break
    assert ibinX >=0, "Can not find the POI {} in the postfit file {}. Please check.".format(poiname, ifilename)

    ibinY = -1
    for binY in range(1, himpact_grouped.GetNbinsY()+1):
        nuis = himpact_grouped.GetYaxis().GetBinLabel(binY)
        if nuis == "stat":
            ibinY = binY
            break
    assert ibinY >=0, "Can not find the stat nuisance in the postfit file {}. Please check.".format(ifilename)
    
    err_stat = himpact_grouped.GetBinContent(ibinX, ibinY)
    
    val_poi, err_poi = GetPOIValue(ifilename, poiname)
    
    err_syst = math.sqrt(err_poi**2 - err_stat**2)
    
    return (err_stat / val_poi, err_syst / val_poi)

def GetTotError(xsecs):
    return math.sqrt(xsecs[1]**2 + xsecs[2]**2 + xsecs[3]**2)

def GetTotRelError(xsecs):
    return GetTotError(xsecs) / xsecs[0]


for sqrtS in ["13TeV", "5TeV"]:
    xsecs[sqrtS] = OrderedDict()
    xsecs_ratios[sqrtS] = OrderedDict()
    
    for ch in ["lepplus", "lepminus", "leplep"]:
        mu, mu_err = GetPOIValue(f_fit, f"{ch}_{sqrtS}_sig_mu")
        #print("mu value for", ch, sqrtS, " is ", mu, " +/- ", mu_err)

        if not doFiducial:
            # use the fiducial central val as the inclusive central value
            mu, _ = GetPOIValue(f_fit_fid, f"{ch}_{sqrtS}_sig_mu")
        
        xsec_MC = xsecs_MC[sqrtS]['born'][maps[ch]]
        #print("xsec value for", ch, sqrtS, "is", xsec_MC)
        
        #xsecs[sqrtS][ch] = (mu * xsec_MC, math.sqrt(mu_err**2 + lumiUncs[sqrtS]**2) * xsec_MC)
        
        mu_stat_err, mu_syst_err = getRelStatSystError(f_fit, f"{ch}_{sqrtS}_sig_mu", "nuisance_group_impact_mu")
        xsecs[sqrtS][ch] = (mu * xsec_MC, mu_stat_err * xsec_MC, mu_syst_err * xsec_MC, lumiUncs[sqrtS] * xsec_MC)
        #print("xsec fit value for", ch, sqrtS, "is", xsecs[sqrtS][ch])
        #print("stat, syst unc ", mu_stat_err, mu_syst_err)
        
    #sys.exit(0)
        
    mu_err = getRelError(f_fit, f"Winc_{sqrtS}_sig_sumxsec")
    mu_stat_err, mu_syst_err = getRelStatSystError(f_fit, f"Winc_{sqrtS}_sig_sumxsec", "nuisance_group_impact_sumpois")
    xsec_tot_W = xsecs[sqrtS]["lepplus"][0] + xsecs[sqrtS]["lepminus"][0]
    xsecs[sqrtS]["winc"] = (xsec_tot_W, mu_stat_err * xsec_tot_W, mu_syst_err * xsec_tot_W, lumiUncs[sqrtS] * xsec_tot_W)
    #xsecs[sqrtS]["Winc"] = (xsec_tot_W, math.sqrt(mu_err**2 + lumiUncs[sqrtS]**2) * xsec_tot_W)
        
    # cross section ratios
    # W+ / W-
    mu_err = getRelError(f_fit, f"WchgRatio_{sqrtS}_ratiometaratio")
    mu_stat_err, mu_syst_err = getRelStatSystError(f_fit, f"WchgRatio_{sqrtS}_ratiometaratio", "nuisance_group_impact_ratiometapois")
    xsec_ratio_Wchg = xsecs[sqrtS]["lepplus"][0] / xsecs[sqrtS]["lepminus"][0]
    #xsecs_ratios[sqrtS]["WchgRatio"] = (xsec_ratio_Wchg, mu_err * xsec_ratio_Wchg, 0.)
    xsecs_ratios[sqrtS]["WchgRatio"] = (xsec_ratio_Wchg, mu_stat_err * xsec_ratio_Wchg, mu_syst_err * xsec_ratio_Wchg,  0.0 * xsec_ratio_Wchg)
    
    # W / Z
    mu_err = getRelError(f_fit, f"WZRatio_{sqrtS}_ratiometaratio")
    mu_stat_err, mu_syst_err = getRelStatSystError(f_fit, f"WZRatio_{sqrtS}_ratiometaratio", "nuisance_group_impact_ratiometapois")
    xsec_ratio_WZ = xsecs[sqrtS]["winc"][0] / xsecs[sqrtS]["leplep"][0]
    #xsecs_ratios[sqrtS]["WZRatio"] = (xsec_ratio_WZ, mu_err * xsec_ratio_WZ, 0.)
    xsecs_ratios[sqrtS]["WZRatio"] = (xsec_ratio_WZ, mu_stat_err * xsec_ratio_WZ, mu_syst_err * xsec_ratio_WZ, 0.0 * xsec_ratio_WZ)
    
# measure the ratios between sqrtS
mu_err = getRelError(f_fit, "sqrtS_Wplus_ratio_ratiometaratio")
mu_stat_err, mu_syst_err = getRelStatSystError(f_fit, "sqrtS_Wplus_ratio_ratiometaratio", "nuisance_group_impact_ratiometapois")
xsec_ratio_sqrtS_Wplus = xsecs["13TeV"]["lepplus"][0] / xsecs["5TeV"]["lepplus"][0]
xsecs_ratios_sqrtS["lepplus"] = (xsec_ratio_sqrtS_Wplus, mu_stat_err * xsec_ratio_sqrtS_Wplus, mu_syst_err * xsec_ratio_sqrtS_Wplus,lumiUncs["sqrtS"] * xsec_ratio_sqrtS_Wplus)
#xsecs_ratios_sqrtS["lepplus"] = (xsec_ratio_sqrtS_Wplus, math.sqrt(mu_err**2 + lumiUncs["sqrtS"]**2) * xsec_ratio_sqrtS_Wplus)

mu_err = getRelError(f_fit, "sqrtS_Wminus_ratio_ratiometaratio")
mu_stat_err, mu_syst_err = getRelStatSystError(f_fit, "sqrtS_Wminus_ratio_ratiometaratio", "nuisance_group_impact_ratiometapois")
xsec_ratio_sqrtS_Wminus = xsecs["13TeV"]["lepminus"][0] / xsecs["5TeV"]["lepminus"][0]
xsecs_ratios_sqrtS["lepminus"] = (xsec_ratio_sqrtS_Wminus, mu_stat_err * xsec_ratio_sqrtS_Wminus, mu_syst_err * xsec_ratio_sqrtS_Wminus, lumiUncs["sqrtS"] * xsec_ratio_sqrtS_Wminus)
#xsecs_ratios_sqrtS["lepminus"] = (xsec_ratio_sqrtS_Wminus, math.sqrt(mu_err**2 + lumiUncs["sqrtS"]**2) * xsec_ratio_sqrtS_Wminus)

mu_err = getRelError(f_fit, "sqrtS_Zinc_ratio_ratiometaratio")
mu_stat_err, mu_syst_err = getRelStatSystError(f_fit, "sqrtS_Zinc_ratio_ratiometaratio", "nuisance_group_impact_ratiometapois")
xsec_ratio_sqrtS_Z = xsecs["13TeV"]["leplep"][0] / xsecs["5TeV"]["leplep"][0]
xsecs_ratios_sqrtS["leplep"] = (xsec_ratio_sqrtS_Z, mu_stat_err * xsec_ratio_sqrtS_Z, mu_syst_err * xsec_ratio_sqrtS_Z, lumiUncs["sqrtS"] * xsec_ratio_sqrtS_Z)
#xsecs_ratios_sqrtS["leplep"] = (xsec_ratio_sqrtS_Z, math.sqrt(mu_err**2 + lumiUncs["sqrtS"]**2) * xsec_ratio_sqrtS_Z)

mu_err = getRelError(f_fit, "sqrtS_Winc_ratio_ratiometaratio")
mu_stat_err, mu_syst_err = getRelStatSystError(f_fit, "sqrtS_Winc_ratio_ratiometaratio", "nuisance_group_impact_ratiometapois")
xsec_ratio_sqrtS_Winc = xsecs["13TeV"]["winc"][0] / xsecs["5TeV"]["winc"][0]
xsecs_ratios_sqrtS["winc"] = (xsec_ratio_sqrtS_Winc, mu_stat_err * xsec_ratio_sqrtS_Winc, mu_syst_err * xsec_ratio_sqrtS_Winc, lumiUncs["sqrtS"] * xsec_ratio_sqrtS_Winc)
#xsecs_ratios_sqrtS["winc"] = (xsec_ratio_sqrtS_Winc, math.sqrt(mu_err**2 + lumiUncs["sqrtS"]**2) * xsec_ratio_sqrtS_Winc)


# double ratios between sqrtS
mu_err = getRelError(f_fit, "sqrtS_WchgRatio_ratio_doubleratiometaratio")
mu_stat_err, mu_syst_err = getRelStatSystError(f_fit, "sqrtS_WchgRatio_ratio_doubleratiometaratio", "nuisance_group_impact_doubleratiometapois")
xsec_ratio_sqrtS_Wchg = xsecs_ratios["13TeV"]["WchgRatio"][0] / xsecs_ratios["5TeV"]["WchgRatio"][0]
xsecs_ratios_sqrtS_ratios["WchgRatio"] = (xsec_ratio_sqrtS_Wchg, mu_stat_err * xsec_ratio_sqrtS_Wchg, mu_syst_err * xsec_ratio_sqrtS_Wchg, 0.)
#xsecs_ratios_sqrtS_ratios["WchgRatio"] = (xsec_ratio_sqrtS_Wchg, mu_err * xsec_ratio_sqrtS_Wchg)


mu_err = getRelError(f_fit, "sqrtS_WZRatio_ratio_doubleratiometaratio")
mu_stat_err, mu_syst_err = getRelStatSystError(f_fit, "sqrtS_WZRatio_ratio_doubleratiometaratio", "nuisance_group_impact_doubleratiometapois")
xsec_ratio_sqrtS_WZ = xsecs_ratios["13TeV"]["WZRatio"][0] / xsecs_ratios["5TeV"]["WZRatio"][0]
xsecs_ratios_sqrtS_ratios["WZRatio"] = (xsec_ratio_sqrtS_WZ, mu_stat_err * xsec_ratio_sqrtS_WZ, mu_syst_err * xsec_ratio_sqrtS_WZ, 0.)
#xsecs_ratios_sqrtS_ratios["WZRatio"] = (xsec_ratio_sqrtS_WZ, mu_err * xsec_ratio_sqrtS_WZ)

    
print("results xsecs: ", xsecs)
print("results xsecs ratios: ", xsecs_ratios)
print("results xs ratios sqrtS: ", xsecs_ratios_sqrtS)
print("results xs ratios sqrtS ratios: ", xsecs_ratios_sqrtS_ratios)

pdfsets = ["nnpdf4.0", "nnpdf3.1", "ct18", "msht20"]
    
# compare the theory predictions with the measured values
from theoryResults import *        
if doFiducial:
    xsec_th = xsec_fid 
else:
    xsec_th = xsec_inc
results = OrderedDict()
for sqrtS in ["13TeV", "5TeV"]:
    results[sqrtS] = OrderedDict()
    
    for ch in ["lepplus", "lepminus", "winc", "leplep"]:
        results[sqrtS][ch] = OrderedDict()
        
        results[sqrtS][ch]['Measured'] = xsecs[sqrtS][ch]
        for pdf in pdfsets:
            results[sqrtS][ch][pdf.upper()] = xsec_th[pdf].GetXsec(ch, sqrtS)
            
        results[sqrtS][f"{ch}_diff"] = OrderedDict()
        #results[sqrtS][f"{ch}_diff"]['Measured'] = (1., xsecs[sqrtS][ch][1] / xsecs[sqrtS][ch][0])
        results[sqrtS][f"{ch}_diff"]["Measured"] = (1., GetTotRelError(xsecs[sqrtS][ch]))
        for pdf in pdfsets:
            results[sqrtS][f"{ch}_diff"][pdf.upper()] = ((xsec_th[pdf].GetXsec(ch, sqrtS)[0] / xsecs[sqrtS][ch][0]),(xsec_th[pdf].GetXsec(ch, sqrtS)[1] / xsecs[sqrtS][ch][0]), (xsec_th[pdf].GetXsec(ch, sqrtS)[2] / xsecs[sqrtS][ch][0]))
    
    # measure the V/V ratios at the same sqrtS
    results[f"{sqrtS}_ratios"] = OrderedDict()
    for ch in ["WchgRatio", "WZRatio"]:
        results[f"{sqrtS}_ratios"][ch] = OrderedDict()
        
        results[f"{sqrtS}_ratios"][ch]['Measured'] = xsecs_ratios[sqrtS][ch]
        for pdf in pdfsets:
            results[f"{sqrtS}_ratios"][ch][pdf.upper()] = getattr(xsec_th[pdf], f"Get{ch}")(sqrtS == "13TeV")
            
        results[f"{sqrtS}_ratios"][f"{ch}_diff"] = OrderedDict()
        #results[f"{sqrtS}_ratios"][f"{ch}_diff"]['Measured'] = (1.0, xsecs_ratios[sqrtS][ch][1] / xsecs_ratios[sqrtS][ch][0])
        results[f"{sqrtS}_ratios"][f"{ch}_diff"]["Measured"] = (1.0, GetTotRelError(xsecs_ratios[sqrtS][ch]))
        
        for pdf in pdfsets:
            results[f"{sqrtS}_ratios"][f"{ch}_diff"][pdf.upper()] = ((getattr(xsec_th[pdf], f"Get{ch}")(sqrtS == "13TeV")[0] / xsecs_ratios[sqrtS][ch][0] ), (getattr(xsec_th[pdf], f"Get{ch}")(sqrtS == "13TeV")[1] / xsecs_ratios[sqrtS][ch][0]), (getattr(xsec_th[pdf], f"Get{ch}")(sqrtS == "13TeV")[2] / xsecs_ratios[sqrtS][ch][0]))
            
# measure the V xsec ratios between sqrtS
results["sqrtS_ratios"] = OrderedDict()
for ch in ["lepplus", "lepminus", "winc", "leplep"]:
    results["sqrtS_ratios"][ch] = OrderedDict()
    
    results['sqrtS_ratios'][ch]['Measured'] = xsecs_ratios_sqrtS[ch]
    
    for pdf in pdfsets:
        results['sqrtS_ratios'][ch][pdf.upper()] = xsec_th[pdf].GetSqrtSRatio(ch)
        #print("sqrtS ratios", results['sqrtS_ratios'][ch][pdf.upper()])
    
    results['sqrtS_ratios'][f"{ch}_diff"] = OrderedDict()
    #results['sqrtS_ratios'][f'{ch}_diff']['Measured'] = (1., xsecs_ratios_sqrtS[ch][1] / xsecs_ratios_sqrtS[ch][0])
    results['sqrtS_ratios'][f'{ch}_diff']['Measured'] = (1., GetTotRelError(xsecs_ratios_sqrtS[ch]))
    
    for pdf in pdfsets:
        results['sqrtS_ratios'][f"{ch}_diff"][pdf.upper()] = ((xsec_th[pdf].GetSqrtSRatio(ch)[0] / xsecs_ratios_sqrtS[ch][0] ), (xsec_th[pdf].GetSqrtSRatio(ch)[1] / xsecs_ratios_sqrtS[ch][0] ), (xsec_th[pdf].GetSqrtSRatio(ch)[2] / xsecs_ratios_sqrtS[ch][0] ))
        
# measure the V/V ratios between sqrtS
results["sqrtS_double_ratios"] = OrderedDict()
for ch in ["WchgRatio", "WZRatio"]:
    results["sqrtS_double_ratios"][ch] = OrderedDict()
    
    results['sqrtS_double_ratios'][ch]['Measured'] = xsecs_ratios_sqrtS_ratios[ch]
    
    for pdf in pdfsets:
        results['sqrtS_double_ratios'][ch][pdf.upper()] = xsec_th[pdf].GetSqrtSDoubleRatio(ch)
    
    results['sqrtS_double_ratios'][f"{ch}_diff"] = OrderedDict()
    #results['sqrtS_double_ratios'][f"{ch}_diff"]['Measured'] = (1., xsecs_ratios_sqrtS_ratios[ch][1] / xsecs_ratios_sqrtS_ratios[ch][0])
    results['sqrtS_double_ratios'][f"{ch}_diff"]['Measured'] = (1., GetTotRelError(xsecs_ratios_sqrtS_ratios[ch]))
    
    for pdf in pdfsets:
        results['sqrtS_double_ratios'][f"{ch}_diff"][pdf.upper()] = ((xsec_th[pdf].GetSqrtSDoubleRatio(ch)[0] / xsecs_ratios_sqrtS_ratios[ch][0]),  (xsec_th[pdf].GetSqrtSDoubleRatio(ch)[1] / xsecs_ratios_sqrtS_ratios[ch][0]), (xsec_th[pdf].GetSqrtSDoubleRatio(ch)[2] / xsecs_ratios_sqrtS_ratios[ch][0]))
        
        
outdir = "plots/xsecResults_fid/" if doFiducial else "plots/xsecResults_inc/"
    
from modules.Utils import FormatTable, FormatROOTInput
from modules.postFitScripts import WriteOutputToText

def dict2Table(list_results):
    results_new = OrderedDict()
    for results in list_results:
        for key, val in results.items():
            if "diff" in key:
                continue
            results_new[key] = val
    return results_new

outputs = FormatTable(dict2Table([results['13TeV'],results['13TeV_ratios']]), doTranspose=True)
print(outputs)
WriteOutputToText(outputs, f"{outdir}/tables/results_13TeV_all.tex")

outputs = FormatTable(dict2Table([results['5TeV'],results['5TeV_ratios']]), doTranspose=True)
print(outputs)
WriteOutputToText(outputs, f"{outdir}/tables/results_5TeV_all.tex")

outputs = FormatTable(dict2Table([results['sqrtS_ratios'], results['sqrtS_double_ratios']]), doTranspose=True)
print(outputs)
WriteOutputToText(outputs, f"{outdir}/tables/results_sqrtS_all.tex")

def DrawCompGraph(xsecs_diffs, outputname, ymin = -0.95, ymax = 1.05, is5TeV = False, doCombineYear = False, canH = 600, canW=1200):
    markers = [2, 3, 4, 5, 25, 26, 27]
    colors = [3, 7, 4, 6, 9, 46]
    tgraphs = OrderedDict()
    
    tgraphs["Measured"] = ROOT.TGraphErrors()
    tgraphs["Measured"].SetFillColorAlpha(15, 0.5)
    
    for j,pdf in enumerate(pdfsets):
        pdf = pdf.upper()
        tgraphs[pdf] = ROOT.TGraphAsymmErrors()
        tgraphs[pdf].SetMarkerStyle(markers[j])
        tgraphs[pdf].SetMarkerSize(2.0)
        tgraphs[pdf].SetMarkerColor(colors[j])
        tgraphs[pdf].SetLineColor(colors[j])
        tgraphs[pdf].SetLineWidth(2)
    
    i = 0    
    binlabels = []
    for ch in xsecs_diffs:
        if "diff" not in ch: 
            continue
        tgraphs["Measured"].SetPoint(i, i + 1, xsecs_diffs[ch]["Measured"][0])
        tgraphs["Measured"].SetPointError(i, 0.5, xsecs_diffs[ch]["Measured"][1])
        binlabels.append(FormatROOTInput(ch.replace("_diff", "")))
        for j, pdf in enumerate(pdfsets):
            pdf = pdf.upper()
            tgraphs[pdf].SetPoint(i, i + 1 + j * 0.1, xsecs_diffs[ch][pdf][0])
            #tgraphs[pdf].SetPointError(i, 0., xsecs_diffs[ch][pdf][1])
            tgraphs[pdf].SetPointEXlow(i, 0.)
            tgraphs[pdf].SetPointEXhigh(i, 0.)
            tgraphs[pdf].SetPointEYlow(i, xsecs_diffs[ch][pdf][1])
            tgraphs[pdf].SetPointEYhigh(i, xsecs_diffs[ch][pdf][2])
            #print("error ", ch, pdf, xsecs_diffs[ch][pdf][1])
        i += 1
            
    DrawHistos(list(tgraphs.values()), list(tgraphs.keys()), 0.5, len(binlabels) + 0.5, "POIs", ymin, ymax, "ratio (Theory / Measured)", outputname, dology=False, binlabels = binlabels, drawoptions = ["E2"] + ["PE1"] * len(pdfsets), legendoptions = ["F"] + ["P"] * len(pdfsets), legendPos = [0.67, 0.7, 0.87, 0.9], is5TeV = is5TeV, doCombineYear = doCombineYear, canH = canH, canW = canW)
    

def DrawHorizontalCompGraph(xsecs_diffs, outputname, xmin = 0.95, xmax = 1.05, is5TeV = False, doCombineYear = False, canH = 600, canW=800, showXsecs = True, pdfsets_plot = ["NNPDF4.0"], poiname = "#sigma^{tot}"):
    markers = [21, 22, 23, 24, 25, 26, 27]
    colors = [4, 3, 7, 6, 9, 46]
    tgraphs = OrderedDict()
    
    counts = 0
    for ch in xsecs_diffs:
        if "diff" not in ch:
            continue
        counts += 1
    
    tgraphs["Measured"] = ROOT.TGraphErrors()
    tgraphs["Measured"].SetFillColorAlpha(15, 0.6)
    
    for j,pdf in enumerate(pdfsets_plot):
        pdf = pdf.upper()
        tgraphs[pdf] = ROOT.TGraphAsymmErrors()
        tgraphs[pdf].SetMarkerStyle(markers[j])
        tgraphs[pdf].SetMarkerSize(1.0)
        tgraphs[pdf].SetMarkerColor(colors[j])
        tgraphs[pdf].SetLineColor(colors[j])
        tgraphs[pdf].SetLineWidth(2)
    
    i = 0    
    procNames = []
    splitLines = []
    valMeasureds = []
    valTheorys = []
    
    splitLine = ROOT.TLine(xmin, counts - i + 0.5, xmax, counts - i + 0.5)
    splitLine.SetLineColor(ROOT.kBlack)
    splitLine.SetLineStyle(1)
    splitLine.SetLineWidth(1)
    splitLines.append(splitLine)
    
    for ch in xsecs_diffs:
        if "diff" not in ch: 
            continue
        tgraphs["Measured"].SetPoint(i, xsecs_diffs[ch]["Measured"][0],counts - i)
        tgraphs["Measured"].SetPointError(i, xsecs_diffs[ch]["Measured"][1], 0.5)
        for j, pdf in enumerate(pdfsets_plot):
            pdf = pdf.upper()
            tgraphs[pdf].SetPoint(i, xsecs_diffs[ch][pdf][0], counts - i - j * 0.15 + 0.2 )
            tgraphs[pdf].SetPointEXlow(i, xsecs_diffs[ch][pdf][1])
            tgraphs[pdf].SetPointEXhigh(i, xsecs_diffs[ch][pdf][2])
            tgraphs[pdf].SetPointEYlow(i, 0.0)
            tgraphs[pdf].SetPointEYhigh(i, 0.0)
            
        splitLine = ROOT.TLine(xmin, counts - i - 0.5, xmax, counts - i - 0.5)
        splitLine.SetLineColor(ROOT.kBlack)
        splitLine.SetLineStyle(2)
        splitLine.SetLineWidth(1)
        splitLines.append(splitLine)
        
        procName = ROOT.TPaveText(0.09, (counts - i - 0.18) / (counts + 1.5) * 0.81 + .08, 0.29, (counts - i + 0.18) / (counts+1.5)*0.81+0.08, "NDC")
        #print("y1, y2", (counts - i - 0.2) / (counts + 1) - .05, (counts - i + 0.2) / (counts + 1) - 0.05)
        procName.SetFillStyle(0)
        procName.SetBorderSize(0)
        procName.SetTextAlign(12)
        procName.AddText(FormatROOTInput(ch.replace("_diff", "")))
        procName.SetTextColor(ROOT.kBlack)
        procName.SetTextSize(0.03)
        procNames.append(procName)
        
        if showXsecs:
            unit = " pb" if (not ("Over" in ch or "Ratio" in ch) and not doCombineYear) else ""
        
            valMeasured = ROOT.TPaveText(0.52, (counts - i ) / (counts + 1.5) * 0.81 + .07, 0.89, (counts - i +0.38) / (counts+1.5)*0.81+0.07, "NDC")
            valMeasured.SetFillStyle(0)
            valMeasured.SetBorderSize(0)
            valMeasured.SetTextAlign(12)
            val = xsecs_diffs[ch.replace("_diff","")]["Measured"]
            rval = roundToError(val)
            #print("measured : ", val)
            if not "Ratio" in ch and "Over" not in ch:
                valMeasured.AddText(f"{rval[0]}#pm{rval[1]}_{{stat}}#pm{rval[2]}_{{syst}}#pm{rval[3]}_{{lum}}{unit}")
            else:
                valMeasured.AddText(f"{rval[0]}#pm{rval[1]}_{{stat}}#pm{rval[2]}_{{syst}}{unit}")
            valMeasured.SetTextColor(ROOT.kBlack)
            valMeasured.SetTextSize(0.03)
            valMeasureds.append(valMeasured)
        
            valTheory = ROOT.TPaveText(0.52, (counts - i - 0.38) / (counts + 1.5) * 0.81 + .07, 0.89, (counts - i ) / (counts+1.5)*0.81+0.07, "NDC")
            valTheory.SetFillStyle(0)
            valTheory.SetBorderSize(0)
            valTheory.SetTextAlign(12)
            val = xsecs_diffs[ch.replace("_diff","")][pdfsets_plot[0]]
            rval = roundToError(val)
            #print("theory ", val)
            valTheory.AddText(f"{rval[0]}^{{+{rval[2]}}}_{{-{rval[1]}}}{unit}")
            valTheory.SetTextColor(colors[0])
            valTheory.SetTextSize(0.03)
            valTheorys.append(valTheory)
        
        i += 1
    
    theoLine = ROOT.TLine(1., 0.5, 1., counts + 0.5)
    theoLine.SetLineColor(ROOT.kRed)
    theoLine.SetLineStyle(1)
    theoLine.SetLineWidth(3)
    
    #headerText = ROOT.TPaveText(0.15, (counts +0.9 ) / (counts + 1.5) * 0.81 + .08, 0.30, (counts + 1.4) / (counts+1.5)*0.81+0.08, "NDC")
    #headerText.SetFillStyle(0)
    #headerText.SetBorderSize(0)
    #headerText.SetTextAlign(12)
    #headerText.AddText(header)
    #headerText.SetTextColor(ROOT.kBlack)
    #headerText.SetTextSize(0.035)
            
    DrawHistos(list(tgraphs.values()), list(tgraphs.keys()), xmin, xmax, "Ratio (Theory / Measured) of " + poiname, 0.5, len(procNames) + 0.5 + 1.5, "", outputname, dology=False, drawoptions = ["E2"] + ["1PE"] * len(pdfsets_plot), legendoptions = ["F"] + ["PLF"] * len(pdfsets_plot), legendPos = [0.10, 0.80, 0.90, 0.9], is5TeV = is5TeV, doCombineYear = doCombineYear, canH = canH, canW = canW, yndivisions=0, additionalToDraw=[theoLine]+splitLines+procNames+valMeasureds+valTheorys, legendNCols = 5, leftmargin=0.05)
    
if doFiducial:
    poiname = "#sigma^{fid}"
else:
    poiname = "#sigma^{tot}"
toplot = deepcopy(results['13TeV'])
toplot.update(results['13TeV_ratios'])
#DrawHorizontalCompGraph(toplot, f"{outdir}/plots/results_13TeV", 0.85, 1.15)
DrawHorizontalCompGraph(toplot, f"{outdir}/plots/results_13TeV", 0.85, 1.20, pdfsets_plot=["NNPDF3.1", "NNPDF4.0", "CT18", "MSHT20"], showXsecs = True, poiname = poiname)


toplot = deepcopy(results['5TeV'])
toplot.update(results['5TeV_ratios'])
DrawHorizontalCompGraph(toplot, f"{outdir}/plots/results_5TeV", 0.85, 1.20, is5TeV = True, pdfsets_plot=["NNPDF3.1", "NNPDF4.0", "CT18", "MSHT20"], showXsecs = True, poiname = poiname)

toplot = deepcopy(results['sqrtS_ratios'])
toplot.update(results['sqrtS_double_ratios'])
DrawHorizontalCompGraph(toplot, f"{outdir}/plots/results_sqrtS", 0.85, 1.20, doCombineYear = True, pdfsets_plot=["NNPDF3.1", "NNPDF4.0", "CT18", "MSHT20"], showXsecs = True, poiname = poiname)
            
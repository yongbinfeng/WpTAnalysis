"""
script to measure FR and carry out closure test
"""
import ROOT
import numpy as np
import sys
from collections import OrderedDict
from modules.Binnings import mass_bins_test
from CMSPLOTS.myFunction import DrawHistos, MultiplyH2, PositiveProtection2D, IntegralAndError2D, CombineOneBin2D, GetRatioPanel, LHistos2Hist, IncludeOverflow2D

ROOT.gROOT.SetBatch(True)

def DrawDataMCStack(hdata, hmc, hpred, xmin = 0, xmax = 140, xlabel = "m_{T} GeV", outputname = "test", is5TeV = False):
    hdata.SetLineColor(1)
    hdata.SetMarkerColor(1)
    hmc.SetLineColor(92)
    hmc.SetFillColor(92)
    hmc.SetMarkerColor(92)
    hpred.SetLineColor(226)
    hpred.SetMarkerColor(226)
    hpred.SetFillColor(226)
    hsmc = ROOT.THStack(f"hsmc_{outputname}", f"hsmc_{outputname}")
    hsmc.Add(hpred)
    hsmc.Add(hmc) 
    hratiopanel = GetRatioPanel(hsmc) 
    DrawHistos([hdata, hsmc], ["Obs", "MC", "Pred QCD"], xmin, xmax, xlabel, 0., None, "Events / 0.1", f"{outputname}", dology=False, showratio=True, yrmin = 0.89, yrmax = 1.11, ratiobase=1, redrawihist=0, hratiopanel=hratiopanel, is5TeV=is5TeV)
    
def GetFR(hnum, hden, projX = False, projY = False):
    assert projX + projY < 2, "only one projection can be done at a time"
    hFR = hnum.Clone(f"{hnum.GetName()}_FR")
    
    if projX:
        # 1D histo along x axis
        hnumN = hnum.ProjectionX(f"{hnum.GetName()}_ProjX")
        hdenN = hden.ProjectionX(f"{hden.GetName()}_ProjX")
        hnumN.Divide(hdenN)
        for iy in range(1, hFR.GetNbinsY()+1):
            for ix in range(1, hFR.GetNbinsX()+1):
                hFR.SetBinContent(ix, iy, hnumN.GetBinContent(ix))
                hFR.SetBinError(ix, iy, hnumN.GetBinError(ix))
    elif projY:
        # 1D histo along y axis
        hnumN = hnum.ProjectionY(f"{hnum.GetName()}_ProjY")
        hdenN = hden.ProjectionY(f"{hden.GetName()}_ProjY")
        hnumN.Divide(hdenN)
        for ix in range(1, hFR.GetNbinsX()+1):
            for iy in range(1, hFR.GetNbinsY()+1):
                hFR.SetBinContent(ix, iy, hnumN.GetBinContent(iy))
                hFR.SetBinError(ix, iy, hnumN.GetBinError(iy))
    else:
        hFR.Divide(hden)
        
    return hFR
    
doPtVsEta = False
doPtVsMet = True
is5TeV = False
doElectron = False

doPtVsDeltaPhi = not doPtVsEta and not doPtVsMet
lepname = "mu" if not doElectron else "e"
sqrtS = "13TeV" if not is5TeV else "5TeV"

if doPtVsEta:
    iname = f"root/output_qcdLepPtVsEtaMean_{lepname}nu_{sqrtS}.root"
    var1 = f"{lepname}_pt"
    var2 = "eta"
    x1label = "p_{T}^{l} [GeV]"
    x1min = 25
    x1max = 70
    x2label = "#eta^{l}"
    x2min = -2.4
    x2max = 2.4  
elif doPtVsMet:
    iname = f"root/output_qcdLepPtVsMetMean_{lepname}nu_{sqrtS}.root"
    var1 = f"{lepname}_pt"
    var2 = "met"
    x1label = "p_{T}^{l} [GeV]"
    x1min = 25
    x1max = 70
    x2label = "p_{T}^{miss} [GeV]"
    x2min = 0
    x2max = 70
else:
    iname = f"root/output_qcdLepPtVsDeltaPhiMean_{lepname}nu_{sqrtS}.root"
    var1 = f"{lepname}_pt"
    var2 = "deltaPhi"
    x1label = "p_{T}^{l} [GeV]"
    x1min = 25
    x1max = 70
    x2label = "#Delta#phi"
    x2min = 0.
    x2max = ROOT.TMath.Pi()
    
var2D = f"{var1}_vs_{var2}"
outdir = f"plots/{sqrtS}/FRAndClosure/{var2D}/"

f = ROOT.TFile.Open(iname) 

isogroups = OrderedDict()
if not doElectron:
    isogroups["SR"] = ["iso0", "iso1", "iso2", "iso3"]
    isogroups["CR0"] = ["iso6", "iso7", "iso8", "iso9", "iso10", "iso11", "iso12", "iso13"]
    #isogroups["CR1"] = ["iso9", "iso10", "iso11", "iso12"]
    #isogroups["CR2"] = ["iso13", "iso14", "iso15", "iso16"]
else:
    isogroups["SR"] = ["iso3", "iso4"]
    isogroups["CR1"] = ["iso6", "iso7", "iso8", "iso9", "iso10"]
    
etabins = ["lepEta_bin0"]

wptbins = ["WpT_bin0"]

chgbins = [lepname+"plus", lepname + "minus"]

mass_bins = mass_bins_test[0]
mtbins = []
for imt in range(len(mass_bins)-1):
    mtbins.append(f"mt{imt}")
mtgroups = OrderedDict()
mtgroups["CR1MT"] = ["mt0"] 
for mtbin in mtbins:
    if mtbin in mtgroups["CR1MT"]:
        continue
    mtgroups[mtbin] = [mtbin]
mass_bins_new = np.array([0.])
igroup = 0
for imt in range(1, len(mass_bins)):
    if imt < len(mtgroups["CR1MT"]):
        continue
    mass_bins_new = np.append(mass_bins_new, mass_bins[imt])
print(mass_bins_new)
    
    
mc_suffix = "MC"
data_suffix = "Data"

# get the histograms
histos_count = OrderedDict()
histos_count_Data = OrderedDict()
histos_count_MC = OrderedDict()
histos_count_QCD = OrderedDict()
histos_FR = OrderedDict()
histos_CT = OrderedDict()
histos_ToSave = []

for wpt in wptbins:
    histos_count[wpt] = OrderedDict()
    histos_count_Data[wpt] = OrderedDict()
    histos_count_MC[wpt] = OrderedDict()
    histos_count_QCD[wpt] = OrderedDict()
    histos_FR[wpt] = OrderedDict()
    histos_CT[wpt] = OrderedDict()
    for lepeta in etabins:
        histos_count[wpt][lepeta] = OrderedDict()
        histos_count_Data[wpt][lepeta] = OrderedDict()
        histos_count_MC[wpt][lepeta] = OrderedDict()
        histos_count_QCD[wpt][lepeta] = OrderedDict()
        histos_FR[wpt][lepeta] = OrderedDict()
        histos_CT[wpt][lepeta] = OrderedDict()
        for chg in chgbins:
            histos_count[wpt][lepeta][chg] = OrderedDict()
            histos_count_Data[wpt][lepeta][chg] = OrderedDict()
            histos_count_MC[wpt][lepeta][chg] = OrderedDict()
            histos_count_QCD[wpt][lepeta][chg] = OrderedDict()
            histos_FR[wpt][lepeta][chg] = OrderedDict()
            histos_CT[wpt][lepeta][chg] = OrderedDict()
            for mtg in mtgroups:
                histos_count[wpt][lepeta][chg][mtg] = OrderedDict()
                histos_count_Data[wpt][lepeta][chg][mtg] = OrderedDict()
                histos_count_MC[wpt][lepeta][chg][mtg] = OrderedDict()
                histos_count_QCD[wpt][lepeta][chg][mtg] = OrderedDict()
                histos_FR[wpt][lepeta][chg][mtg] = OrderedDict()
                histos_CT[wpt][lepeta][chg][mtg] = OrderedDict()
               
                # read and combine all histograms in different mt and iso groups
                for mt in mtgroups[mtg]: 
                    for isog in isogroups:
                        hdata = None
                        hmc = None
                        for iso in isogroups[isog]:
                            strname = "weight_{}_{}_{}_{}_{}".format(chg, iso, wpt, lepeta, mt)
                            hname = f"histo_wjets_{var2D}_{strname}"

                            if not hdata:
                                print("Getting " + hname + f"_{data_suffix}")
                                hdata = f.Get(hname + f"_{data_suffix}").Clone(hname+"_PureData")
                            else:
                                hdata.Add(f.Get(hname + f"_{data_suffix}"))
                            if not hmc:
                                hmc = f.Get(hname + f"_{mc_suffix}").Clone(hname+"_PureMC")
                            else:
                                hmc.Add(f.Get(hname + f"_{mc_suffix}"))

                        # inlcude overflow and underflow bins 
                        IncludeOverflow2D(hdata, True)
                        IncludeOverflow2D(hmc, True)

                        # count = count_Data - count_MC
                        h = hdata.Clone(hname + f"_Subtracted")
                        h.Add(hmc, -1)
                        #PositiveProtection2D(h)

                        if isog not in histos_count[wpt][lepeta][chg][mtg]:
                            histos_count[wpt][lepeta][chg][mtg][isog] = h
                            histos_count_Data[wpt][lepeta][chg][mtg][isog] = hdata
                            histos_count_MC[wpt][lepeta][chg][mtg][isog] = hmc       
                        else:
                            histos_count[wpt][lepeta][chg][mtg][isog].Add(h) 
                            histos_count_Data[wpt][lepeta][chg][mtg][isog].Add(hdata)
                            histos_count_MC[wpt][lepeta][chg][mtg][isog].Add(hmc)
                    
            # calculate FR in different mt and iso groups, 
            # as a function of var1 and var2        
            for mtg in mtgroups:
                for isog in isogroups: 
                    PositiveProtection2D(histos_count[wpt][lepeta][chg][mtg][isog])
                    # calculate FR / transfer factor
                    hnum = histos_count[wpt][lepeta][chg][mtg]["SR"]
                    #hnum = h.Clone(h.GetName() + "_FR")
                    hden = histos_count[wpt][lepeta][chg][mtg][isog] 
                    #hnum.Divide(hden)

                    hFR = GetFR(hnum, hden, projY=0, projX=0)
                    histos_FR[wpt][lepeta][chg][mtg][isog] = hFR

                    DrawHistos([hFR], [], x1min, x1max, x1label, x2min, x2max, x2label, f"{outdir}/{isog}/FR/histo_wjets_{var2D}_FR_{chg}_{wpt}_{lepeta}_{mtg}_{isog}", False, False, False, dologz=False, doth2=True, drawoptions="COLZ,text", is5TeV=is5TeV) 
                    
            hdatas_var1 = []
            hmcs_var1 = []
            hpreds_var1 = []
            hdatas_var2 = []
            hmcs_var2 = []
            hpreds_var2 = []
            for isog in isogroups: 
                if isog == "SR":
                    continue
                for mtg in mtgroups:
                    hfr = histos_FR[wpt][lepeta][chg]["CR1MT"][isog]
                    hcounts_cr = histos_count[wpt][lepeta][chg][mtg][isog]
                    hcounts_sr = histos_count[wpt][lepeta][chg][mtg]["SR"]

                    h = hcounts_cr.Clone(hcounts_cr.GetName() + "_CT")

                    hpred = MultiplyH2(h, hfr)
                    histos_count_QCD[wpt][lepeta][chg][mtg]['SR'] = hpred

                    hpred_var1 = hpred.ProjectionX(f"{hpred.GetName()}_Proj_{var1}")
                    hcounts_sr_var1 = hcounts_sr.ProjectionX(f"{hcounts_sr.GetName()}_Proj_{var1}")

                    hpred_var2 = hpred.ProjectionY(f"{hpred.GetName()}_Proj_{var2}")
                    hcounts_sr_var2 = hcounts_sr.ProjectionY(f"{hcounts_sr.GetName()}_Proj_{var2}")

                    hratio = hpred.Clone(hpred.GetName() + "_ratio")
                    hratio.Divide(hcounts_sr)

                    # compare predicted QCD with data-MC
                    DrawHistos([hcounts_sr_var1, hpred_var1], ["Obs", "Pred"], x1min, x1max, x1label, 0, None, "Events", f"{outdir}/{isog}/Closure_mT/histos_wjets_{var1}_CT_{chg}_{wpt}_{lepeta}_{mtg}", dology=False, showratio=True, yrmin = 0.8, yrmax = 1.2, is5TeV = is5TeV, mycolors=[ROOT.kBlack, ROOT.kRed])

                    DrawHistos([hcounts_sr_var2, hpred_var2], ["Obs", "Pred"], x2min, x2max, x2label, 0., None, "Events", f"{outdir}/{isog}/Closure_mT/histos_wjets_{var2}_CT_{chg}_{wpt}_{lepeta}_{mtg}", dology=False, showratio=True, yrmin = 0.8, yrmax = 1.2, is5TeV = is5TeV, mycolors=[ROOT.kBlack, ROOT.kRed])

                    # make the stacked histograms 
                    # so that the signal contribution can be visualized
                    hdata = histos_count_Data[wpt][lepeta][chg][mtg]["SR"]
                    hmc = histos_count_MC[wpt][lepeta][chg][mtg]["SR"]  

                    hdata_var2 = hdata.ProjectionY(f"{hdata.GetName()}_Proj_{var2}_for{isog}")
                    hmc_var2 = hmc.ProjectionY(f"{hmc.GetName()}_Proj_{var2}_for{isog}")
                    DrawDataMCStack(hdata_var2, hmc_var2, hpred_var2, x2min, x2max, x2label, f"{outdir}/{isog}/Closure_mT_Stacked/histos_wjets_{var2}_CT_{chg}_{wpt}_{lepeta}_{mtg}_stacked", is5TeV = is5TeV)
                    hdatas_var2.append(hdata_var2)
                    hmcs_var2.append(hmc_var2)
                    hpreds_var2.append(hpred_var2)

                    hdata_var1 = hdata.ProjectionX(f"{hdata.GetName()}_Proj_{var1}_for{isog}")
                    hmc_var1 = hmc.ProjectionX(f"{hmc.GetName()}_Proj_{var1}_for{isog}")
                    DrawDataMCStack(hdata_var1, hmc_var1, hpred_var1, x1min, x1max, x1label, f"{outdir}/{isog}/Closure_mT_Stacked/histos_wjets_{var1}_CT_{chg}_{wpt}_{lepeta}_{mtg}_stacked", is5TeV = is5TeV)
                    hdatas_var1.append(hdata_var1)
                    hmcs_var1.append(hmc_var1)
                    hpreds_var1.append(hpred_var1)

                    # draw 2D closure ratio
                    DrawHistos([hratio], [], x1min, x1max, x1label, x2min, x2max, x2label, f"{outdir}/{isog}/Closure_mT_2D/histo_wjets_{var2D}_CT_{chg}_{wpt}_{lepeta}_{mtg}", False, False, False, dologz=False, doth2=True, drawoptions="COLZ,text", zmin=0.5, zmax=1.2, is5TeV=is5TeV)

                # check lepton pt distribution
                strname = f"{lepname}_{var1}_{chg}_{wpt}_{lepeta}_{isog}"
                hdata_var1 = LHistos2Hist(hdatas_var1, f"hdata_{strname}_Data")
                hmc_var1 = LHistos2Hist(hmcs_var1, f"hmc_{strname}_MC")
                hpred_var1 = LHistos2Hist(hpreds_var1, f"hqcd_{strname}_QCD")
                DrawDataMCStack(hdata_var1, hmc_var1, hpred_var1, x1min, x1max, x1label, f"{outdir}/{isog}/histos_wjets_{var1}_CT_{chg}_{wpt}_{lepeta}_stacked", is5TeV = is5TeV)
            
                # check lepton eta distribution
                strname = f"{lepname}_{var2}_{chg}_{wpt}_{lepeta}_{isog}"
                hdata_var2 = LHistos2Hist(hdatas_var2, f"hdata_{strname}_Data")
                hmc_var2 = LHistos2Hist(hmcs_var2, f"hmc_{strname}_MC")
                hpred_var2 = LHistos2Hist(hpreds_var2, f"hqcd_{strname}_QCD")
                DrawDataMCStack(hdata_var2, hmc_var2, hpred_var2, x2min, x2max, x2label, f"{outdir}/{isog}/histos_wjets_{var2}_CT_{chg}_{wpt}_{lepeta}_stacked", is5TeV = is5TeV)

                # check the mT distribution
                hdata_mt = ROOT.TH1D(f"hdata_{lepname}_mT_{chg}_{wpt}_{lepeta}_{isog}_Data", f"hdata_{lepname}_mT_{chg}_{wpt}_{lepeta}_Data", len(mass_bins_new)-1, mass_bins_new) 
                hdata_mt.Sumw2()
                hmc_mt = ROOT.TH1D(f"hmc_{lepname}_mT_{chg}_{wpt}_{lepeta}_{isog}_MC", f"hmc_{lepname}_mT_{chg}_{wpt}_{lepeta}_MC", len(mass_bins_new)-1, mass_bins_new)
                hmc_mt.Sumw2()
                hpred_mt = ROOT.TH1D(f"hqcd_{lepname}_mT_{chg}_{wpt}_{lepeta}_{isog}_QCD", f"hqcd_{lepname}_mT_{chg}_{wpt}_{lepeta}_QCD", len(mass_bins_new)-1, mass_bins_new)
                hpred_mt.Sumw2()
                for imtg, mtg in enumerate(mtgroups.keys()):
                    val_data, err_data = IntegralAndError2D(histos_count_Data[wpt][lepeta][chg][mtg]["SR"])
                    val_mc, err_mc = IntegralAndError2D(histos_count_MC[wpt][lepeta][chg][mtg]["SR"])
                    val_pred, err_pred = IntegralAndError2D(histos_count_QCD[wpt][lepeta][chg][mtg]["SR"])
                    hdata_mt.SetBinContent(imtg + 1, val_data)
                    hdata_mt.SetBinError(imtg + 1, err_data)
                    hmc_mt.SetBinContent(imtg + 1, val_mc)
                    hmc_mt.SetBinError(imtg + 1, err_mc)
                    hpred_mt.SetBinContent(imtg + 1, val_pred)
                    hpred_mt.SetBinError(imtg + 1, err_pred)

                DrawDataMCStack(hdata_mt, hmc_mt, hpred_mt, 0, 140, "m_{T}", f"{outdir}/{isog}/histos_wjets_{lepname}_mt_CT_{chg}_{wpt}_{lepeta}_stacked", is5TeV = is5TeV)
            
                histos_ToSave.append(hdata_mt)
                histos_ToSave.append(hmc_mt)
                histos_ToSave.append(hpred_mt)
            
ofile = ROOT.TFile.Open(f"{outdir}/histos_qcdFR_{lepname}_{sqrtS}.root", "RECREATE")
for h in histos_ToSave:
    h.SetDirectory(ofile)
    h.Write()
ofile.Close()
"""
script to measure FR and carry out closure test
"""
import ROOT
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
    ymax = hdata.GetMaximum() * 1.3
    hratiopanel = GetRatioPanel(hsmc) 
    DrawHistos([hdata, hsmc], ["Obs", "MC", "Pred QCD"], xmin, xmax, xlabel, 0., ymax, "Events / 0.1", f"{outdir}/{outputname}", dology=False, showratio=True, yrmin = 0.89, yrmax = 1.11, ratiobase=1, redrawihist=0, hratiopanel=hratiopanel, is5TeV=is5TeV)
        
doPtVsEta = True
is5TeV = False
doElectron = True

doPtVsMet = not doPtVsEta
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
else:
    iname = f"root/output_qcdLepPtVsMetMean_{lepname}nu_{sqrtS}.root"
    var1 = f"{lepname}_pt"
    var2 = "met"
    x1label = "p_{T}^{l} [GeV]"
    x1min = 25
    x1max = 70
    x2label = "p_{T}^{miss} [GeV]"
    x2min = 0
    x2max = 70
    
var2D = f"{var1}_vs_{var2}"
outdir = f"plots/{sqrtS}/FRAndClosure/{var2D}/"

# base mt for calculate FR
basemt = "mt0"

f = ROOT.TFile.Open(iname) 

isoCuts = [0.0,0.001, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50,
           0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95]
isobins = []
for isobin in range(0, 20):
    isobins.append(f"iso{isobin}")
    
isogroups = OrderedDict()
if not doElectron:
    isogroups["SR"] = ["iso0", "iso1", "iso2", "iso3"]
    isogroups["antiIso1"] = ["iso9", "iso10", "iso11", "iso12", "iso13", "iso14", "iso15"]
else:
    isogroups["SR"] = ["iso3", "iso4"]
    isogroups["antiIso1"] = ["iso7", "iso8", "iso9", "iso10"]
    
etabins = ["lepEta_bin0"]

wptbins = ["WpT_bin0"]

chgbins = [lepname+"plus", lepname + "minus"]

mass_bins = mass_bins_test[0]
mtbins = []
for imt in range(len(mass_bins)-1):
    mtbins.append(f"mt{imt}")
    
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
            for mt in mtbins:
                histos_count[wpt][lepeta][chg][mt] = OrderedDict()
                histos_count_Data[wpt][lepeta][chg][mt] = OrderedDict()
                histos_count_MC[wpt][lepeta][chg][mt] = OrderedDict()
                histos_count_QCD[wpt][lepeta][chg][mt] = OrderedDict()
                histos_FR[wpt][lepeta][chg][mt] = OrderedDict()
                histos_CT[wpt][lepeta][chg][mt] = OrderedDict()
                
                for isogroup in isogroups:
                    hdata = None
                    hmc = None
                    for iso in isogroups[isogroup]:
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
                    
                    h = hdata.Clone(hname + f"_Subtracted")
                    h.Add(hmc, -1)
                    PositiveProtection2D(h)
                    
                    histos_count[wpt][lepeta][chg][mt][isogroup] = h
                    histos_count_Data[wpt][lepeta][chg][mt][isogroup] = hdata
                    histos_count_MC[wpt][lepeta][chg][mt][isogroup] = hmc
               
                for isogroup in isogroups: 
                    # calculate FR
                    h = histos_count[wpt][lepeta][chg][mt]["SR"]
                    hnum = h.Clone(h.GetName() + "_FR")
                    
                    hden = histos_count[wpt][lepeta][chg][mt][isogroup] 
                    hnum.Divide(hden)
                    
                    histos_FR[wpt][lepeta][chg][mt][isogroup] = hnum 
                    
                    DrawHistos([hnum], [], x1min, x1max, x1label, x2min, x2max, x2label, f"{outdir}/histo_wjets_{var2D}_FR_{chg}_{wpt}_{lepeta}_{mt}_{isogroup}", False, False, False, dologz=False, doth2=True, drawoptions="COLZ,text", is5TeV=is5TeV) 
                    
            hdatas_lep_pt = []
            hmcs_lep_pt = []
            hpreds_lep_pt = []
            hdatas_lep_eta = []
            hmcs_lep_eta = []
            hpreds_lep_eta = []
            
            for mt in mtbins:
                hfr = histos_FR[wpt][lepeta][chg][basemt]["antiIso1"]
                hcounts_cr = histos_count[wpt][lepeta][chg][mt]["antiIso1"]
                
                h = hcounts_cr.Clone(hcounts_cr.GetName() + "_CT")
                
                hpred = MultiplyH2(h, hfr)
                histos_count_QCD[wpt][lepeta][chg][mt]['SR'] = hpred
                
                hcounts_sr = histos_count[wpt][lepeta][chg][mt]["SR"]
                hpred_lep_pt = hpred.ProjectionX(f"{hpred.GetName()}_Proj_lep_pt")
                hcounts_sr_pt = hcounts_sr.ProjectionX(f"{hcounts_sr.GetName()}_Proj_lep_pt")
                
                hpred_lep_eta = hpred.ProjectionY(f"{hpred.GetName()}_Proj_lep_eta")
                hcounts_sr_eta = hcounts_sr.ProjectionY(f"{hcounts_sr.GetName()}_Proj_lep_eta")
                
                hratio = hpred.Clone(hpred.GetName() + "_ratio")
                hratio.Divide(hcounts_sr)
                
                hcounts_sr_eta.SetLineColor(ROOT.kBlack)
                hcounts_sr_eta.SetMarkerColor(ROOT.kBlack)
                hpred_lep_eta.SetLineColor(ROOT.kRed)
                hpred_lep_eta.SetMarkerColor(ROOT.kRed)
                
                # compare predicted QCD with data-MC
                ymax = max(hcounts_sr_eta.GetMaximum(), hpred_lep_eta.GetMaximum()) * 1.2
                DrawHistos([hcounts_sr_eta, hpred_lep_eta], ["Obs", "Pred"], x2min, x2max, x2label, 0., ymax, "Events / 0.1", f"{outdir}/histos_wjets_{lepname}_lep_eta_CT_{chg}_{wpt}_{lepeta}_{mt}", dology=False, showratio=True, yrmin = 0.8, yrmax = 1.2, is5TeV = is5TeV)
                
                hcounts_sr_pt.SetLineColor(ROOT.kBlack)
                hcounts_sr_pt.SetMarkerColor(ROOT.kBlack)
                hpred_lep_pt.SetLineColor(ROOT.kRed)
                hpred_lep_pt.SetMarkerColor(ROOT.kRed)
                ymax = max(hcounts_sr_pt.GetMaximum(), hpred_lep_pt.GetMaximum()) * 1.2
                DrawHistos([hcounts_sr_pt, hpred_lep_pt], ["Obs", "Pred"], x1min, x1max, x1label, 0, ymax, "Events / GeV", f"{outdir}/histos_wjets_{lepname}_lep_pt_CT_{chg}_{wpt}_{lepeta}_{mt}", dology=False, showratio=True, yrmin = 0.8, yrmax = 1.2, is5TeV = is5TeV)
                
                # make the stacked histograms 
                # so that the signal contribution can be visualized
                hdata = histos_count_Data[wpt][lepeta][chg][mt]["SR"]
                hmc = histos_count_MC[wpt][lepeta][chg][mt]["SR"]  
                
                hdata_lep_eta = hdata.ProjectionY(f"{hdata.GetName()}_Proj_lep_eta")
                hmc_lep_eta = hmc.ProjectionY(f"{hmc.GetName()}_Proj_lep_eta")
                DrawDataMCStack(hdata_lep_eta, hmc_lep_eta, hpred_lep_eta, x2min, x2max, x2label, f"histos_wjets_{lepname}_lep_eta_CT_{chg}_{wpt}_{lepeta}_{mt}_stacked", is5TeV = is5TeV)
                hdatas_lep_eta.append(hdata_lep_eta)
                hmcs_lep_eta.append(hmc_lep_eta)
                hpreds_lep_eta.append(hpred_lep_eta)
                
                hdata_lep_pt = hdata.ProjectionX(f"{hdata.GetName()}_Proj_lep_pt")
                hmc_lep_pt = hmc.ProjectionX(f"{hmc.GetName()}_Proj_lep_pt")
                DrawDataMCStack(hdata_lep_pt, hmc_lep_pt, hpred_lep_pt, x1min, x1max, x1label, f"histos_wjets_{lepname}_lep_pt_CT_{chg}_{wpt}_{lepeta}_{mt}_stacked", is5TeV = is5TeV)
                hdatas_lep_pt.append(hdata_lep_pt)
                hmcs_lep_pt.append(hmc_lep_pt)
                hpreds_lep_pt.append(hpred_lep_pt)
                
                # draw 2D closure ratio
                DrawHistos([hratio], [], x1min, x1max, x1label, x2min, x2max, x2label, f"{outdir}/histo_wjets_{var2D}_CT_{chg}_{wpt}_{lepeta}_{mt}", False, False, False, dologz=False, doth2=True, drawoptions="COLZ,text", zmin=0.5, zmax=1.2, is5TeV=is5TeV)
                
            # check lepton pt distribution
            strname = f"{lepname}_lep_pt_{chg}_{wpt}_{lepeta}"
            hdata_lep_pt = LHistos2Hist(hdatas_lep_pt, f"hdata_{strname}_Data")
            hmc_lep_pt = LHistos2Hist(hmcs_lep_pt, f"hmc_{strname}_MC")
            hpred_lep_pt = LHistos2Hist(hpreds_lep_pt, f"hqcd_{strname}_QCD")
            DrawDataMCStack(hdata_lep_pt, hmc_lep_pt, hpred_lep_pt, x1min, x1max, x1label, f"histos_wjets_{lepname}_lep_pt_CT_{chg}_{wpt}_{lepeta}_stacked", is5TeV = is5TeV)
            
            # check lepton eta distribution
            strname = f"{lepname}_lep_eta_{chg}_{wpt}_{lepeta}"
            hdata_lep_eta = LHistos2Hist(hdatas_lep_eta, f"hdata_{strname}_Data")
            hmc_lep_eta = LHistos2Hist(hmcs_lep_eta, f"hmc_{strname}_MC")
            hpred_lep_eta = LHistos2Hist(hpreds_lep_eta, f"hqcd_{strname}_QCD")
            DrawDataMCStack(hdata_lep_eta, hmc_lep_eta, hpred_lep_eta, x2min, x2max, x2label, f"histos_wjets_{lepname}_lep_eta_CT_{chg}_{wpt}_{lepeta}_stacked", is5TeV = is5TeV)
                
            # check the mT distribution
            hdata_mt = ROOT.TH1D(f"hdata_{lepname}_mT_{chg}_{wpt}_{lepeta}_Data", f"hdata_{lepname}_mT_{chg}_{wpt}_{lepeta}_Data", len(mass_bins)-1, mass_bins) 
            hdata_mt.Sumw2()
            hmc_mt = ROOT.TH1D(f"hmc_{lepname}_mT_{chg}_{wpt}_{lepeta}_MC", f"hmc_{lepname}_mT_{chg}_{wpt}_{lepeta}_MC", len(mass_bins)-1, mass_bins)
            hmc_mt.Sumw2()
            hpred_mt = ROOT.TH1D(f"hqcd_{lepname}_mT_{chg}_{wpt}_{lepeta}_QCD", f"hqcd_{lepname}_mT_{chg}_{wpt}_{lepeta}_QCD", len(mass_bins)-1, mass_bins)
            hpred_mt.Sumw2()
            for imt, mt in enumerate(mtbins):
                val_data, err_data = IntegralAndError2D(histos_count_Data[wpt][lepeta][chg][mt]["SR"])
                val_mc, err_mc = IntegralAndError2D(histos_count_MC[wpt][lepeta][chg][mt]["SR"])
                val_pred, err_pred = IntegralAndError2D(histos_count_QCD[wpt][lepeta][chg][mt]["SR"])
                hdata_mt.SetBinContent(imt + 1, val_data)
                hdata_mt.SetBinError(imt + 1, err_data)
                hmc_mt.SetBinContent(imt + 1, val_mc)
                hmc_mt.SetBinError(imt + 1, err_mc)
                hpred_mt.SetBinContent(imt + 1, val_pred)
                hpred_mt.SetBinError(imt + 1, err_pred)

            DrawDataMCStack(hdata_mt, hmc_mt, hpred_mt, 0, 140, "m_{T}", f"histos_wjets_{lepname}_mt_CT_{chg}_{wpt}_{lepeta}_stacked", is5TeV = is5TeV)
            
            histos_ToSave.append(hdata_mt)
            histos_ToSave.append(hmc_mt)
            histos_ToSave.append(hpred_mt)
            
ofile = ROOT.TFile.Open(f"{outdir}/histos_qcdFR_{lepname}_{sqrtS}.root", "RECREATE")
for h in histos_ToSave:
    h.SetDirectory(ofile)
    h.Write()
ofile.Close()
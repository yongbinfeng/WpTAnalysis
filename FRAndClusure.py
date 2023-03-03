"""
script to measure FR and carry out closure test
"""
import ROOT
from collections import OrderedDict
from modules.Binnings import mass_bins_test
from CMSPLOTS.myFunction import DrawHistos

ROOT.gROOT.SetBatch(True)

def MultiplyH2(h1, h2):
    assert h1.GetNbinsX() == h2.GetNbinsX(), "h1 and h2 have different number of bins in x"
    assert h1.GetNbinsY() == h2.GetNbinsY(), "h1 and h2 have different number of bins in y"
    
    for ix in range(1, h1.GetNbinsX()+1):
        for iy in range(1, h1.GetNbinsY()+1):
            val1 = h1.GetBinContent(ix, iy)
            val2 = h2.GetBinContent(ix, iy)
            err1 = h1.GetBinError(ix, iy)
            err2 = h2.GetBinError(ix, iy)
            h1.SetBinContent(ix, iy, val1*val2)
            h1.SetBinError(ix, iy, ROOT.TMath.Sqrt((val1*err2)**2 + (val2*err1)**2))
    return h1

def IntegralAndError(h2):
    val = 0
    err2 = 0
    for ibinx in range(1, h2.GetNbinsX()+1):
        for ibiny in range(1, h2.GetNbinsY()+1):
            val += h2.GetBinContent(ibinx, ibiny)
            err2 += h2.GetBinError(ibinx, ibiny)**2
    return val, ROOT.TMath.Sqrt(err2)

def CombineOneBin2D(h, ix, iy, jx, jy):
    """
    combien the j-th bin to i-th bin, and clean j-th bin
    """
    h.SetBinContent(ix, iy, h.GetBinContent(ix, iy) + h.GetBinContent(jx, jy))
    h.SetBinError(ix, iy, ROOT.TMath.Sqrt(h.GetBinError(ix, iy)**2 + h.GetBinError(jx, jy)**2))
    # clean the original bin so that there would not be double counting
    h.SetBinContent(jx, jy, 0)
    h.SetBinError(jx, jy, 0)
    
def GetRatioPanel(hs):
    """
    build the ratio error bars used for THStack
    Assuming different histograms are uncorrelated
    """
    hlist = hs.GetHists()
    hratio = None
    for h in hlist:
        if not hratio:
            hratio = h.Clone(h.GetName() + "_ratio")
        else:
            hratio.Add(h)
    hratio.Divide(hratio)
    hratio.SetFillColor(15)
    return hratio

def IncludeOverflow(h2, doUnderflow=False):
    """
    this might not work for one th2 with only one bin in one direction
    """
    nbinsX = h2.GetNbinsX()
    nbinsY = h2.GetNbinsY()
    
    for ix in range(1, nbinsX+1):
        CombineOneBin2D(h2, ix, nbinsY, ix, nbinsY+1)
    for iy in range(1, nbinsY+1):
        CombineOneBin2D(h2, nbinsX, iy, nbinsX+1, iy)
    # correct the corner
    CombineOneBin2D(h2, nbinsX, nbinsY, nbinsX+1, nbinsY+1)
        
    if doUnderflow:
        for ix in range(1, nbinsX+1):
            CombineOneBin2D(h2, ix, 1, ix, 0)
        for iy in range(1, nbinsY+1):
            CombineOneBin2D(h2, 1, iy, 0, iy)
        # correct the corner
        CombineOneBin2D(h2, 1, 1, 0, 0)
        CombineOneBin2D(h2, 1, nbinsY, 0, nbinsY+1)
        CombineOneBin2D(h2, nbinsX, 1, nbinsX+1, 0)

iname = "root/output_qcdLepPtVsEtaMean_munu_13TeV.root"
outdir = "plots/13TeV/FRAndClosure/"

basemt = "mt0"

f = ROOT.TFile.Open(iname) 

isoCuts = [0.0,0.001, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50,
           0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95]
isobins = []
for isobin in range(0, 20):
    isobins.append(f"iso{isobin}")
    
isogroups = OrderedDict()
isogroups["SR"] = ["iso0", "iso1", "iso2", "iso3"]
isogroups["antiIso1"] = ["iso9", "iso10", "iso11", "iso12", "iso13", "iso14", "iso15"]
    
etabins = ["lepEta_bin0"]

wptbins = ["WpT_bin0"]

lepname = "mu" 
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
                        hname = f"histo_wjets_{lepname}_pt_vs_eta_{strname}"
                    
                        if not hdata:
                            hdata = f.Get(hname + f"_{data_suffix}").Clone(hname+"_PureData")
                        else:
                            hdata.Add(f.Get(hname + f"_{data_suffix}"))
                        if not hmc:
                            hmc = f.Get(hname + f"_{mc_suffix}").Clone(hname+"_PureMC")
                        else:
                            hmc.Add(f.Get(hname + f"_{mc_suffix}"))
                            
                    # inlcude overflow and underflow bins 
                    IncludeOverflow(hdata, True)
                    IncludeOverflow(hmc, True)
                    
                    h = hdata.Clone(hname + f"_Subtracted")
                    h.Add(hmc, -1)
                    
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
                    
                    DrawHistos([hnum], [], 25, 70, "p^{l}_{T} [GeV]", -2.4, 2.4, "#eta ", f"{outdir}/histo_wjets_{lepname}_pt_vs_eta_FR_{chg}_{wpt}_{lepeta}_{mt}_{isogroup}", False, False, False, dologz=False, doth2=True, drawoptions="COLZ,text") 
            
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
                
                ymax = max(hcounts_sr_eta.GetMaximum(), hpred_lep_eta.GetMaximum()) * 1.2
                DrawHistos([hcounts_sr_eta, hpred_lep_eta], ["Obs", "Pred"], -2.4, 2.4, "#eta", 0., ymax, "Events / 0.1", f"{outdir}/histos_wjets_{lepname}_lep_eta_CT_{chg}_{wpt}_{lepeta}_{mt}", dology=False, showratio=True, yrmin = 0.8, yrmax = 1.2)
                
                hcounts_sr_pt.SetLineColor(ROOT.kBlack)
                hcounts_sr_pt.SetMarkerColor(ROOT.kBlack)
                hpred_lep_pt.SetLineColor(ROOT.kRed)
                hpred_lep_pt.SetMarkerColor(ROOT.kRed)
                ymax = max(hcounts_sr_pt.GetMaximum(), hpred_lep_pt.GetMaximum()) * 1.2
                DrawHistos([hcounts_sr_pt, hpred_lep_pt], ["Obs", "Pred"], 25, 70, "p^{l}_{T} [GeV]", 0, ymax, "Events / GeV", f"{outdir}/histos_wjets_{lepname}_lep_pt_CT_{chg}_{wpt}_{lepeta}_{mt}", dology=False, showratio=True, yrmin = 0.8, yrmax = 1.2)
                
                # make the stacked histograms so that the signal contribution can be visualized
                hdata = histos_count_Data[wpt][lepeta][chg][mt]["SR"]
                hmc = histos_count_MC[wpt][lepeta][chg][mt]["SR"]  
                
                hdata_lep_eta = hdata.ProjectionY(f"{hdata.GetName()}_Proj_lep_eta")
                hdata_lep_eta.SetLineColor(1)
                hdata_lep_eta.SetMarkerColor(1)
                hmc_lep_eta = hmc.ProjectionY(f"{hmc.GetName()}_Proj_lep_eta")
                hmc_lep_eta.SetLineColor(92)
                hmc_lep_eta.SetFillColor(92)
                hmc_lep_eta.SetMarkerColor(92)
                hpred_lep_eta.SetLineColor(226)
                hpred_lep_eta.SetMarkerColor(226)
                hpred_lep_eta.SetFillColor(226)
                hsmc_lep_eta = ROOT.THStack("hsmc_lep_eta", "hsmc_lep_eta")
                hsmc_lep_eta.Add(hpred_lep_eta)
                hsmc_lep_eta.Add(hmc_lep_eta)
                ymax = hdata_lep_eta.GetMaximum() * 1.3
                hratiopanel = GetRatioPanel(hsmc_lep_eta) 
                DrawHistos([hdata_lep_eta, hsmc_lep_eta], ["Obs", "MC", "Pred QCD"], -2.4, 2.4, "#eta", 0., ymax, "Events / 0.1", f"{outdir}/histos_wjets_{lepname}_lep_eta_CT_{chg}_{wpt}_{lepeta}_{mt}_stacked", dology=False, showratio=True, yrmin = 0.8, yrmax = 1.2, ratiobase=1, redrawihist=0, hratiopanel=hratiopanel)
                
                hdata_lep_pt = hdata.ProjectionX(f"{hdata.GetName()}_Proj_lep_pt")
                hdata_lep_pt.SetLineColor(1)
                hdata_lep_pt.SetMarkerColor(1)
                hmc_lep_pt = hmc.ProjectionX(f"{hmc.GetName()}_Proj_lep_pt")
                hmc_lep_pt.SetLineColor(92)
                hmc_lep_pt.SetFillColor(92)
                hmc_lep_pt.SetMarkerColor(92)
                hpred_lep_pt.SetLineColor(226)
                hpred_lep_pt.SetMarkerColor(226)
                hpred_lep_pt.SetFillColor(226)
                hsmc_lep_pt = ROOT.THStack("hsmc_lep_pt", "hsmc_lep_pt")
                hsmc_lep_pt.Add(hpred_lep_pt)
                hsmc_lep_pt.Add(hmc_lep_pt)
                ymax = hdata_lep_pt.GetMaximum() * 1.3
                hratiopanel = GetRatioPanel(hsmc_lep_pt)
                DrawHistos([hdata_lep_pt, hsmc_lep_pt], ["Obs", "MC", "Pred QCD"], 25, 70, "p^{l}_{T} [GeV]", 0, ymax, "Events / GeV", f"{outdir}/histos_wjets_{lepname}_lep_pt_CT_{chg}_{wpt}_{lepeta}_{mt}_stacked", dology=False, showratio=True, yrmin = 0.8, yrmax = 1.2, ratiobase=1, redrawihist=0, hratiopanel=hratiopanel)
                
                
                # draw 2D closure ratio
                DrawHistos([hratio], [], 25, 70, "p^{l}_{T} [GeV]", -2.4, 2.4, "#eta ", f"{outdir}/histo_wjets_{lepname}_pt_vs_eta_CT_{chg}_{wpt}_{lepeta}_{mt}", False, False, False, dologz=False, doth2=True, drawoptions="COLZ,text", zmin=0.5, zmax=1.2)
                
            # check the mT distribution
            hdata_mt = ROOT.TH1D(f"hdata_{lepname}_mT_{chg}_{wpt}_{lepeta}_Data", f"hdata_{lepname}_mT_{chg}_{wpt}_{lepeta}_Data", len(mass_bins)-1, mass_bins) 
            hdata_mt.Sumw2()
            hmc_mt = ROOT.TH1D(f"hmc_{lepname}_mT_{chg}_{wpt}_{lepeta}_MC", f"hmc_{lepname}_mT_{chg}_{wpt}_{lepeta}_MC", len(mass_bins)-1, mass_bins)
            hmc_mt.Sumw2()
            hpred_mt = ROOT.TH1D(f"hqcd_{lepname}_mT_{chg}_{wpt}_{lepeta}_QCD", f"hqcd_{lepname}_mT_{chg}_{wpt}_{lepeta}_QCD", len(mass_bins)-1, mass_bins)
            hpred_mt.Sumw2()
            for imt, mt in enumerate(mtbins):
                val_data, err_data = IntegralAndError(histos_count_Data[wpt][lepeta][chg][mt]["SR"])
                val_mc, err_mc = IntegralAndError(histos_count_MC[wpt][lepeta][chg][mt]["SR"])
                val_pred, err_pred = IntegralAndError(histos_count_QCD[wpt][lepeta][chg][mt]["SR"])
                hdata_mt.SetBinContent(imt + 1, val_data)
                hdata_mt.SetBinError(imt + 1, err_data)
                hmc_mt.SetBinContent(imt + 1, val_mc)
                hmc_mt.SetBinError(imt + 1, err_mc)
                hpred_mt.SetBinContent(imt + 1, val_pred)
                hpred_mt.SetBinError(imt + 1, err_pred)

            hdata_mt.SetLineColor(1)
            hdata_mt.SetMarkerColor(1)
            hmc_mt.SetLineColor(92)
            hmc_mt.SetFillColor(92)
            hmc_mt.SetMarkerColor(92)
            hpred_mt.SetLineColor(226)
            hpred_mt.SetMarkerColor(226)
            hpred_mt.SetFillColor(226)
            hsmc_mt = ROOT.THStack("hsmc_mt", "hsmc_mt")
            hsmc_mt.Add(hpred_mt)
            hsmc_mt.Add(hmc_mt) 
            ymax = hdata_mt.GetMaximum() * 1.3
            hratiopanel = GetRatioPanel(hsmc_mt) 
            DrawHistos([hdata_mt, hsmc_mt], ["Obs", "MC", "Pred QCD"], 0, 140, "m_{T}", 0., ymax, "Events / 0.1", f"{outdir}/histos_wjets_{lepname}_mt_CT_{chg}_{wpt}_{lepeta}_stacked", dology=False, showratio=True, yrmin = 0.89, yrmax = 1.11, ratiobase=1, redrawihist=0, hratiopanel=hratiopanel)
            
            histos_ToSave.append(hdata_mt)
            histos_ToSave.append(hmc_mt)
            histos_ToSave.append(hpred_mt)
            
ofile = ROOT.TFile.Open(f"{outdir}/histos_qcdFR.root", "RECREATE")
for h in histos_ToSave:
    h.SetDirectory(ofile)
    h.Write()
ofile.Close()
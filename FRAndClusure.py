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
            h1.SetBinContent(ix, iy, h1.GetBinContent(ix, iy)*h2.GetBinContent(ix, iy))
            
    return h1

def IncludeOverflow(h2, doUnderflow=False):
    nbinsX = h2.GetNbinsX()
    nbinsY = h2.GetNbinsY()
    
    for ix in range(1, nbinsX+1):
        h2.SetBinContent(ix, nbinsY, h2.GetBinContent(ix, nbinsY) + h2.GetBinContent(ix, nbinsY+1))
        h2.SetBinContent(ix, nbinsY+1, 0)
    for iy in range(1, nbinsY+1):
        h2.SetBinContent(nbinsX, iy, h2.GetBinContent(nbinsX, iy) + h2.GetBinContent(nbinsX+1, iy))
        h2.SetBinContent(nbinsX+1, iy, 0)
        
    if doUnderflow:
        for ix in range(1, nbinsX+1):
            h2.SetBinContent(ix, 1, h2.GetBinContent(ix, 1) + h2.GetBinContent(ix, 0))
            h2.SetBinContent(ix, 0, 0)
        for iy in range(1, nbinsY+1):
            h2.SetBinContent(1, iy, h2.GetBinContent(1, iy) + h2.GetBinContent(0, iy))
            h2.SetBinContent(0, iy, 0)

iname = "root/output_qcdLepPtVsEtaMean_munu_5TeV.root"

outdir = "plots/13TeV/FRAndClosure/"

f = ROOT.TFile.Open(iname) 

isoCuts = [0.0,0.001, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50,
           0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95]
isobins = []
for isobin in range(0, 20):
    isobins.append(f"iso{isobin}")
    
isogroups = OrderedDict()
isogroups["SR"] = ["iso0", "iso1", "iso2", "iso3"]
isogroups["antiIso1"] = ["iso5", "iso6", "iso7", "iso8", "iso9", "iso10", "iso11", "iso12", "iso13", "iso14", "iso15"]
    
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
histos_FR = OrderedDict()
histos_CT = OrderedDict()

for wpt in wptbins:
    histos_count[wpt] = OrderedDict()
    histos_FR[wpt] = OrderedDict()
    histos_CT[wpt] = OrderedDict()
    for lepeta in etabins:
        histos_count[wpt][lepeta] = OrderedDict()
        histos_FR[wpt][lepeta] = OrderedDict()
        histos_CT[wpt][lepeta] = OrderedDict()
        for chg in chgbins:
            histos_count[wpt][lepeta][chg] = OrderedDict()
            histos_FR[wpt][lepeta][chg] = OrderedDict()
            histos_CT[wpt][lepeta][chg] = OrderedDict()
            for mt in mtbins:
                histos_count[wpt][lepeta][chg][mt] = OrderedDict()
                histos_FR[wpt][lepeta][chg][mt] = OrderedDict()
                histos_CT[wpt][lepeta][chg][mt] = OrderedDict()
                
                for isogroup in isogroups:
                    h = None
                    for iso in isogroups[isogroup]:
                        strname = "weight_{}_{}_{}_{}_{}".format(chg, iso, wpt, lepeta, mt)
                        hname = f"histo_wjets_{lepname}_pt_vs_eta_{strname}"
                    
                        if not h:
                            h = f.Get(hname + f"_{data_suffix}")
                            h.Add(f.Get(hname + f"_{mc_suffix}"), -1)
                        
                        else:
                            h.Add(f.Get(hname + f"_{data_suffix}"))
                            h.Add(f.Get(hname + f"_{mc_suffix}"), -1)
                           
                        # inlcude overflow and underflow bins 
                        IncludeOverflow(h, True)
                
                    histos_count[wpt][lepeta][chg][mt][isogroup] = h
               
                for isogroup in isogroups: 
                    h = histos_count[wpt][lepeta][chg][mt]["SR"]
                    hnum = h.Clone(h.GetName() + "_FR")
                    
                    hden = histos_count[wpt][lepeta][chg][mt][isogroup] 
                    hnum.Divide(hden)
                    
                    histos_FR[wpt][lepeta][chg][mt][isogroup] = hnum 
                    
                    DrawHistos([hnum], [], 20, 70, "p^{l}_{T} [GeV]", -2.4, 2.4, "#eta ", f"{outdir}/histo_wjets_{lepname}_pt_vs_eta_FR_{chg}_{wpt}_{lepeta}_{mt}_{isogroup}", False, False, False, dologz=False, doth2=True, drawoptions="COLZ,text") 
            
            for mt in mtbins:
                hfr = histos_FR[wpt][lepeta][chg]["mt0"]["antiIso1"]
                hcounts_cr = histos_count[wpt][lepeta][chg][mt]["antiIso1"]
                
                h = hcounts_cr.Clone(hcounts_cr.GetName() + "_CT")
                
                hpred = MultiplyH2(h, hfr)
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
                DrawHistos([hcounts_sr_pt, hpred_lep_pt], ["Obs", "Pred"], 20, 70, "p^{l}_{T} [GeV]", 0, ymax, "Events / GeV", f"{outdir}/histos_wjets_{lepname}_lep_pt_CT_{chg}_{wpt}_{lepeta}_{mt}", dology=False, showratio=True, yrmin = 0.8, yrmax = 1.2)
                
                DrawHistos([hratio], [], 20, 70, "p^{l}_{T} [GeV]", -2.4, 2.4, "#eta ", f"{outdir}/histo_wjets_{lepname}_pt_vs_eta_CT_{chg}_{wpt}_{lepeta}_{mt}", False, False, False, dologz=False, doth2=True, drawoptions="COLZ,text")
                
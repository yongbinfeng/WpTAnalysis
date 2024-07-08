"""
script to measure FR and carry out closure test
"""
import ROOT
import numpy as np
import sys
import argparse
from collections import OrderedDict
from modules.Binnings import mass_bins_forqcd
from CMSPLOTS.myFunction import DrawHistos, MultiplyH2, PositiveProtection, IntegralAndError2D, CombineOneBin2D, GetRatioPanel, LHistos2Hist, IncludeOverflow2D, TH2ToTH1s, SymmetrizeHisto

def DrawDataMCStack(hdata, hmc, hpred = None, xmin = 0, xmax = 140, xlabel = "m_{T} GeV", outputname = "test", is5TeV = False):
    hdata.SetLineColor(1)
    hdata.SetMarkerColor(1)
    hmc.SetLineColor(92)
    hmc.SetFillColor(92)
    hmc.SetMarkerColor(92)
    hsmc = ROOT.THStack(f"hsmc_{outputname}", f"hsmc_{outputname}")
    labels = ["Obs", "MC"]
    if hpred:
        hpred.SetLineColor(226)
        hpred.SetMarkerColor(226)
        hpred.SetFillColor(226)
        hsmc.Add(hpred)
        labels.append("Pred QCD")
    # add hmc later than hpred, as first added will be drawn at the bottom
    hsmc.Add(hmc)
    hratiopanel = GetRatioPanel(hsmc) 
    DrawHistos([hdata, hsmc], labels, xmin, xmax, xlabel, 0., None, "Events", f"{outputname}", dology=False, showratio=True, yrmin = 0.89, yrmax = 1.11, ratiobase=1, redrawihist=0, hratiopanel=hratiopanel, is5TeV=is5TeV)
    
def GetFR(hnums_data, hnums_mc, hdens_data, hdens_mc, projX = False, projY = False, avergeEta = False, suffix = "Clone", scaleNumMC = 1.0, scaleDenMC = 1.0):
    """
    hFr = (hnum_data - hnum_mc * scaleNumMC) / (hden_data - hden_mc * scaleDenMC)
    """
    hnum_data = LHistos2Hist(hnums_data, f"{hnums_data[0].GetName()}_{suffix}")
    hnum_mc = LHistos2Hist(hnums_mc, f"{hnums_mc[0].GetName()}_{suffix}")
    hden_data = LHistos2Hist(hdens_data, f"{hdens_data[0].GetName()}_{suffix}")
    hden_mc = LHistos2Hist(hdens_mc, f"{hdens_mc[0].GetName()}_{suffix}")
    
    hnum = GetAbsYield(hnum_data, hnum_mc, sData = 1.0, sMC = scaleNumMC, suffix = f"num_FR_for_{suffix}")
    hden = GetAbsYield(hden_data, hden_mc, sData = 1.0, sMC = scaleDenMC, suffix = f"den_FR_for_{suffix}")
    
    if avergeEta:
        hnum = AverageEta2D(hnum)
        hden = AverageEta2D(hden)
        
    assert projX + projY < 2, "only one projection can be done at a time"
    hFR = hnum.Clone(f"{hnum.GetName()}_FR")
    
    if projX:
        # 1D histo along x axis
        hnumN = hnum.ProjectionX(f"{hnum.GetName()}_ProjX")
        hdenN = hden.ProjectionX(f"{hden.GetName()}_ProjX")
        PositiveProtection(hnumN)
        PositiveProtection(hdenN)
        hnumN.Divide(hdenN)
        for iy in range(1, hFR.GetNbinsY()+1):
            for ix in range(1, hFR.GetNbinsX()+1):
                hFR.SetBinContent(ix, iy, hnumN.GetBinContent(ix))
                hFR.SetBinError(ix, iy, hnumN.GetBinError(ix))
    elif projY:
        # 1D histo along y axis
        hnumN = hnum.ProjectionY(f"{hnum.GetName()}_ProjY")
        hdenN = hden.ProjectionY(f"{hden.GetName()}_ProjY")
        PositiveProtection(hnumN)
        PositiveProtection(hdenN)
        hnumN.Divide(hdenN)
        for ix in range(1, hFR.GetNbinsX()+1):
            for iy in range(1, hFR.GetNbinsY()+1):
                hFR.SetBinContent(ix, iy, hnumN.GetBinContent(iy))
                hFR.SetBinError(ix, iy, hnumN.GetBinError(iy))
    else:
        PositiveProtection(hnum)
        PositiveProtection(hden)
        hFR.Divide(hden)
    return hFR

def AverageEta2D(h2):
    """
    average eta+ and eta-, assuming eta is on the y axis
    """
    nbinsy = h2.GetNbinsY()
    for ix in range(1, h2.GetNbinsX()+1):
        for iy in range(1, int(nbinsy/2)+1):
            val = (h2.GetBinContent(ix, iy) + h2.GetBinContent(ix, nbinsy-iy+1)) * 0.5
            err = np.sqrt((h2.GetBinError(ix, iy)**2 + h2.GetBinError(ix, nbinsy-iy+1)**2)*0.5)
            h2.SetBinContent(ix, iy, val)
            h2.SetBinError(ix, iy, err)
            h2.SetBinContent(ix, nbinsy-iy+1, val)
            h2.SetBinError(ix, nbinsy-iy+1, err)
    return h2

def StatUnc2SysUnc(h1, prefix = "stat"):
    """
    prepare the stat varied histograms for combine
    """
    hs = []
    for ix in range(1, h1.GetNbinsX()+1):
        hup = h1.Clone(f"{h1.GetName()}_{prefix}_bin{ix}Up")
        hdn = h1.Clone(f"{h1.GetName()}_{prefix}_bin{ix}Down")
        val = h1.GetBinContent(ix)
        err = h1.GetBinError(ix)
        hup.SetBinContent(ix, val+err)
        hdn.SetBinContent(ix, val-err)
        hs.append(hup)
        hs.append(hdn)
    return hs

def GetAbsYield(hdata, hmc, sData = 1.0, sMC = 1.0, suffix=""):
    """
    hasb = (hdata * sData - hmc * sMC)
    """
    habs = hdata.Clone(f"{hdata.GetName()}_{suffix}")
    habs.Scale(sData)
    habs.Add(hmc, -1.0 * sMC)
    return habs
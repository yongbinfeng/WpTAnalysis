"""
script to make postfit comparisons
"""
import ROOT
import os,sys
import numpy as np

sys.path.append("/afs/cern.ch/work/y/yofeng/public/CMSPLOTS")
from myFunction import DrawHistos
from SampleManager import DrawConfig

ROOT.gROOT.SetBatch(True)

def doPlot(ifilename, channel, etabins, postfix):
    ifile = ROOT.TFile(ifilename)
    
    horgdata = ifile.Get("obs")
    
    ## read the postfit plots from input file
    hexpsig = ifile.Get("expproc_w_{}_sig_postfit".format(channel))
    # ewk bkg includes W->tau+nu, z->ll, and diboson process
    hexpewk = ifile.Get("expproc_taunu_postfit")
    hexpewk.Add( ifile.Get("expproc_zxx_postfit"))
    hexpewk.Add( ifile.Get("expproc_VV_postfit"))
    hexpttbar = ifile.Get("expproc_tt_postfit")

    # qcd process might have contributions from a few different etabins
    hexpqcd = ifile.Get("expproc_QCD_{}_{}_postfit".format(channel, etabins[0]))
    for etabin in etabins[1:]:
        hexpqcd.Add( ifile.Get("expproc_QCD_{}_{}_postfit".format(channel, etabin)) )

    # the combined prediction of all processes,
    # which should have included the correct total postfit uncertainties
    hexpfull = ifile.Get("expfull_postfit")

    # the histograms saved in the root file does not follow the original bining
    # recover the original binning
    # hard coded the binning for now, assuming 10GeV fixed bin width
    nbins = horgdata.GetNbinsX()
    xmin = 0.
    xmax = 10 * nbins
    binnings = (nbins, xmin, xmax)
    #binnings = (newbins.shape[0]-1, newbins)
    hdata   = ROOT.TH1D("hdata_{}_{}".format( channel, postfix),  "hdata_{}_{}".format( channel, postfix),  *binnings)
    hsig    = ROOT.TH1D("hsig_{}_{}".format(  channel, postfix),  "hsig_{}_{}".format(  channel, postfix),  *binnings)
    hewk    = ROOT.TH1D("hewk_{}_{}".format(  channel, postfix),  "hewk_{}_{}".format(  channel, postfix),  *binnings)
    httbar  = ROOT.TH1D("httbar_{}_{}".format(channel, postfix),  "httbar_{}_{}".format(channel, postfix),  *binnings)
    hqcd    = ROOT.TH1D("hqcd_{}_{}".format(  channel, postfix),  "hqcd_{}_{}".format(  channel, postfix),  *binnings)
    hratio  = ROOT.TH1D("hrato_{}_{}".format( channel, postfix),  "hratio_{}_{}".format(channel, postfix),  *binnings)
    for ibin in xrange(1, nbins+1):
        hdata.SetBinContent(ibin,   horgdata.GetBinContent(ibin))
        hdata.SetBinError(ibin,     horgdata.GetBinError(ibin))
        hsig.SetBinContent(ibin,    hexpsig.GetBinContent(ibin))
        hewk.SetBinContent(ibin,    hexpewk.GetBinContent(ibin))
        httbar.SetBinContent(ibin,  hexpttbar.GetBinContent(ibin))
        hqcd.SetBinContent(ibin,    hexpqcd.GetBinContent(ibin))
        hqcd.SetBinError(ibin, hexpqcd.GetBinError(ibin))

        hratio.SetBinContent(ibin, hexpfull.GetBinContent(ibin))
        hratio.SetBinError(ibin,   hexpfull.GetBinError(ibin))

    # deal with the uncertainty bar
    for ibin in xrange(1, hratio.GetNbinsX()+1):
        val = hratio.GetBinContent(ibin)
        err = hratio.GetBinError(ibin)
        hratio.SetBinContent(ibin, 1.0)
        if val!=0:
            hratio.SetBinError(ibin, err/val)
        else:
            hratio.SetBinError(ibin, 0.)

    hsig.SetFillColor(92)
    hsig.SetLineColor(92)

    hewk.SetFillColor(216)
    hewk.SetLineColor(216)

    httbar.SetFillColor(96)
    httbar.SetLineColor(96)

    hqcd.SetFillColor(226)
    hqcd.SetLineColor(226)

    hdata.SetMarkerStyle(20)
    hdata.SetMarkerSize(1)
    hdata.SetMarkerColor(1)
    hdata.SetLineColor(1)
    hdata.Scale(1.0, "width")
    hqcd.Scale(1.0, "width")
    httbar.Scale(1.0, "width")
    hewk.Scale(1.0, "width")
    hsig.Scale(1.0, "width")

    hs_gmc = ROOT.THStack("hs_stack_{}_{}".format(channel, postfix), "hs_stack")
    hs_gmc.Add(hqcd)
    hs_gmc.Add(httbar)
    hs_gmc.Add(hewk)
    hs_gmc.Add(hsig)

    ymaxs = {"pm": 3e4, "mm": 2e4}
    ymaxs = {"muplus": 3e4, "muminus": 2e4, "eplus": 1.5e4, "eminus": 1.0e4}
    "histo_wjets_eplus_mT_1_WpT_bin0_lepEta_bin1_data"
    drawconfigs = DrawConfig(xmin = xmin, xmax = xmax, xlabel = "m_{T} [GeV]", ymin = 0, ymax = ymaxs[channel], ylabel = "Events / GeV", outputname = "histo_wjets_{}_mT_PostFit_{}".format(channel, postfix), dology=False, addOverflow=False, addUnderflow=False, yrmin=0.95, yrmax=1.05)
    DrawHistos( [hdata, hs_gmc], ["Data", "Signal", "EWK", "ttbar", "QCD"], drawconfigs.xmin, drawconfigs.xmax, drawconfigs.xlabel, drawconfigs.ymin, drawconfigs.ymax, drawconfigs.ylabel, drawconfigs.outputname, dology=drawconfigs.dology, dologx=drawconfigs.dologx, showratio=drawconfigs.showratio, yrmax = drawconfigs.yrmax, yrmin = drawconfigs.yrmin, yrlabel = drawconfigs.yrlabel, donormalize=drawconfigs.donormalize, ratiobase=drawconfigs.ratiobase, legendPos = drawconfigs.legendPos, redrawihist = drawconfigs.redrawihist, extraText = drawconfigs.extraText, noCMS = drawconfigs.noCMS, addOverflow = drawconfigs.addOverflow, addUnderflow = drawconfigs.addUnderflow, nMaxDigits = drawconfigs.nMaxDigits, hratiopannel=hratio, drawoptions=['PE', 'HIST same'])

if __name__ == "__main__":
    doPlot("Cards/datacard_muplus_lepEta_bin0.root", "muplus", ["lepEta_bin0"], "")
    doPlot("Cards/datacard_eplus_lepEta.root",  "eplus", ["lepEta_bin1", "lepEta_bin2"], "")

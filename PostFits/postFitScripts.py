"""
script to make postfit comparisons
"""
import ROOT
import os,sys,math
from collections import OrderedDict
import json
import re

from CMSPLOTS.myFunction import DrawHistos
from SampleManager import DrawConfig

ROOT.gROOT.SetBatch(True)

def MakePostPlot(ifilename, channel, etabins, postfix, showpull=False):
    """
    compare the postfit of data and templates
    """
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
    hpull   = ROOT.TH1D("hpull_{}_{}".format( channel, postfix),  "hpull_{}_{}".format( channel, postfix),  *binnings)
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

        diff = horgdata.GetBinContent(ibin) - hexpfull.GetBinContent(ibin)
        # take the sigma as sqrt(data**2 + templates**2)
        # not 100% sure if this is the correct way to calculate pull
        sig = math.sqrt(horgdata.GetBinError(ibin)**2 + hexpfull.GetBinError(ibin)**2)
        hpull.SetBinContent(ibin, diff/(sig+1e-6))

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
    DrawHistos( [hdata, hs_gmc], ["Data", "Signal", "EWK", "ttbar", "QCD"], drawconfigs.xmin, drawconfigs.xmax, drawconfigs.xlabel, drawconfigs.ymin, drawconfigs.ymax, drawconfigs.ylabel, drawconfigs.outputname, dology=drawconfigs.dology, dologx=drawconfigs.dologx, showratio=drawconfigs.showratio, yrmax = drawconfigs.yrmax, yrmin = drawconfigs.yrmin, yrlabel = drawconfigs.yrlabel, donormalize=drawconfigs.donormalize, ratiobase=drawconfigs.ratiobase, legendPos = drawconfigs.legendPos, redrawihist = drawconfigs.redrawihist, extraText = drawconfigs.extraText, noCMS = drawconfigs.noCMS, addOverflow = drawconfigs.addOverflow, addUnderflow = drawconfigs.addUnderflow, nMaxDigits = drawconfigs.nMaxDigits, hratiopannel=hratio, drawoptions=['PE', 'HIST same'], showpull=showpull, hpulls=[hpull])


def result2json(ifilename, poiname, ofilename):
    """
    script to convert the postfit POI and impacts of nuisance parameters 
    to json file, which will be used to make impact plots later
    """
    nameMap = {
        "SysWeight1" : "mc",
        "SysWeight2" : "FSR",
        "SysWeight3" : "bkg",
        "SysWeight4" : "tagpt",
        "SysWeight6" : "Prefire",
        "SysRecoil2" : "recoil_eta",
        "SysRecoil3" : "recoil_keys",
        "SysRecoil6" : "recoil_stat0",
        "SysRecoil7" : "recoil_stat1",
        "SysRecoil8" : "recoil_stat2",
        "SysRecoil9" : "recoil_stat3",
        "SysRecoil10": "recoil_stat4",
        "SysRecoil11": "recoil_stat5",
        "SysRecoil12": "recoil_stat6",
        "SysRecoil13": "recoil_stat7",
        "SysRecoil14": "recoil_stat8",
        "SysRecoil15": "recoil_stat9",
    }

    def getNuisName(nuis):
        if nuis in nameMap.keys():
            return nameMap[nuis]
        elif bool(re.match(r"\w*bin\d+shape", nuis)):
            return "QCD_" + nuis
        else:
            return nuis

    ifile = ROOT.TFile(ifilename)
    himpact = ifile.Get("nuisance_impact_mu")
    himpact_grouped = ifile.Get("nuisance_group_impact_mu")
    tree = ifile.Get("fitresults")
    tree.GetEntry(0)

    # find the POI bin for poiname
    ibinX = -1
    for binX in range(1, himpact.GetNbinsX()+1):
        poi = himpact.GetXaxis().GetBinLabel(binX)
        if poi == poiname:
            ibinX = binX
            continue
    assert ibinX >=0, "Can not find the POI {} in the postfit file {}. Please check.".format(poiname, ifilename)

    results = OrderedDict()
    results['POIs'] = []
    val = getattr(tree, poiname)
    err = abs(getattr(tree, poiname+"_err"))
    poi = OrderedDict()
    poi['fit'] = [val-err, val, val+err]
    poi['name'] = poiname
    results['POIs'].append(poi)

    results['method'] = 'default'
    results['params'] = []

    # dump impacts
    impacts = OrderedDict()
    for ibinY in range(1, himpact.GetNbinsY()+1):
        nuis = himpact.GetYaxis().GetBinLabel(ibinY)
        impacts[nuis] = himpact.GetBinContent(ibinX, ibinY)

    # add the grouped QCD systematic
    for ibinY in range(1, himpact_grouped.GetNbinsY()+1):
        tmpY = himpact_grouped.GetYaxis().GetBinLabel(ibinY)
        if tmpY != 'QCD':
            continue
        impacts['QCD'] = himpact_grouped.GetBinContent(ibinX, ibinY)

    # sort impacts, descending
    impacts = OrderedDict(sorted(impacts.items(), key=lambda x: abs(x[1]), reverse=True))

    pulls = OrderedDict()
    for nuis in impacts.keys():
        if nuis != 'QCD':
            val = getattr(tree, nuis)
            err = getattr(tree, nuis+"_err")
            err = abs(err)
        else:
            # grouped QCD nuisance
            val = 0.
            err = 1.
        pulls[nuis] = [val - err, val, val + err]

    # save to results
    for nuis in impacts.keys():
        systematic = OrderedDict()
        systematic['fit'] = pulls[nuis]
        systematic['groups'] = []
        systematic['impact_' + poiname] = impacts[nuis]
        systematic['name'] = getNuisName(nuis)
        systematic['prefit'] = [-1.0, 0., 1.0]
        systematic[poiname] = [poi['fit'][1] - impacts[nuis], poi['fit'][1], poi['fit'][1] + impacts[nuis]]
        systematic['type'] = "Gaussian"
        print(getNuisName(nuis), pulls[nuis][1], pulls[nuis][1]-pulls[nuis][0], impacts[nuis])

        results['params'].append(systematic)

    with open(ofilename, 'w') as fp:
        json.dump(results, fp, indent=2)

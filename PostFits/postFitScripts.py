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

def MakePostPlot(ifilename, channel, postfix, showpull=False):
    """
    compare the unrolled postfit of data and templates
    """
    ifile = ROOT.TFile(ifilename)
    
    horgdata = ifile.Get("obs")

    # get the list of histograms saved in the file
    hkeys = ifile.GetListOfKeys()
    hnames_sig = []
    hnames_qcd = []
    for hkey in hkeys:
        if bool(re.match(r"expproc_w_"+channel+"_\w*sig_postfit", hkey.GetName())):
            hnames_sig.append( hkey.GetName() )
        elif bool(re.match(r"expproc_QCD_"+channel+"_\w*postfit", hkey.GetName())):
            hnames_qcd.append( hkey.GetName() )
    assert len(hnames_sig)>=1, "There should be at least one sig histogram in file: {}".format(ifilename)
    assert len(hnames_qcd)>=1, "There should be at least one QCD histogram in file: {}".format(ifilename)
    
    ## read the postfit plots from input file
    hexpsig = ifile.Get(hnames_sig[0])
    for hname_sig in hnames_sig[1:]:
        hexpsig.Add( ifile.Get(hname_sig) )
    # ewk bkg includes W->tau+nu, z->ll, and diboson process
    hexpewk = ifile.Get("expproc_taunu_postfit")
    hexpewk.Add( ifile.Get("expproc_zxx_postfit"))
    hexpewk.Add( ifile.Get("expproc_VV_postfit"))
    hexpttbar = ifile.Get("expproc_tt_postfit")

    # qcd process might have contributions from a few different etabins
    hexpqcd = ifile.Get(hnames_qcd[0])
    for hname_qcd in hnames_qcd[1:]:
        hexpqcd.Add( ifile.Get(hname_qcd) )

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

    nevts = OrderedDict()
    nevts['data'] = hdata.Integral()
    nevts['sig'] = hsig.Integral()
    nevts['ewk'] = hewk.Integral()
    nevts['ttbar'] = httbar.Integral()
    nevts['qcd'] = hqcd.Integral()

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

    ymaxs = {"muplus": 3e4, "muminus": 2e4, "eplus": 1.5e4, "eminus": 1.0e4}
    drawconfigs = DrawConfig(xmin = xmin, xmax = xmax, xlabel = "Unrolled m_{T} [GeV]", ymin = 0, ymax = ymaxs[channel] / (int(nbins/36)+1), ylabel = "Events / GeV", outputname = "histo_wjets_{}_mT_PostFit_{}".format(channel, postfix), dology=False, addOverflow=False, addUnderflow=False, yrmin=0.95, yrmax=1.05)
    DrawHistos( [hdata, hs_gmc], ["Data", "Signal", "EWK", "ttbar", "QCD"], drawconfigs.xmin, drawconfigs.xmax, drawconfigs.xlabel, drawconfigs.ymin, drawconfigs.ymax, drawconfigs.ylabel, drawconfigs.outputname, dology=drawconfigs.dology, dologx=drawconfigs.dologx, showratio=drawconfigs.showratio, yrmax = drawconfigs.yrmax, yrmin = drawconfigs.yrmin, yrlabel = drawconfigs.yrlabel, donormalize=drawconfigs.donormalize, ratiobase=drawconfigs.ratiobase, legendPos = drawconfigs.legendPos, redrawihist = drawconfigs.redrawihist, extraText = drawconfigs.extraText, noCMS = drawconfigs.noCMS, addOverflow = drawconfigs.addOverflow, addUnderflow = drawconfigs.addUnderflow, nMaxDigits = drawconfigs.nMaxDigits, hratiopannel=hratio, drawoptions=['PE', 'HIST same'], showpull=showpull, hpulls=[hpull], W_ref = 600 * int(nbins/36+1))

    return nevts


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

    # add the grouped QCD and Recoil systematic
    groupnames = []
    for ibinY in range(1, himpact_grouped.GetNbinsY()+1):
        tmpY = himpact_grouped.GetYaxis().GetBinLabel(ibinY)
        if tmpY == 'stat':
            continue
        impacts[tmpY] = himpact_grouped.GetBinContent(ibinX, ibinY)
        groupnames.append(tmpY)

    # sort impacts, descending
    impacts = OrderedDict(sorted(impacts.items(), key=lambda x: abs(x[1]), reverse=True))

    pulls = OrderedDict()
    for nuis in impacts.keys():
        if nuis not in groupnames:
            val = getattr(tree, nuis)
            err = getattr(tree, nuis+"_err")
            err = abs(err)
        else:
            # manually set the postfit of the grouped sys to [-1,1], and pulled at 0,
            # since only the impacts are useful to us
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


def MakeWpTPostFitPlots(jsonNames, postfix=""):
    """
    plot the signal strength and systematic unc impacts as a function of wpt
    """
    import numpy as np
    # hard code the w pt binning for now
    wptbins = np.array([0., 8.0, 16.0, 24.0, 32.0, 40.0, 50.0, 70.0, 100.0, 120.0])
    nbins = wptbins.shape[0] - 1

    hmu = ROOT.TH1F("hmu_wpt_{}".format(postfix), "hmu_wpt_{}".format(postfix), nbins, wptbins)

    sysUncs = ['lumi_13TeV', 'Recoil', 'QCD', 'Prefire', 'norm_tt', 'effstat', 'norm_taunu', 'FSR', 'tagpt', 'norm_z']
    colors = {"total": 1, 'lumi_13TeV': 14, "Recoil": 2, "QCD": 3, "Prefire": 4, "effstat":6, "norm_taunu": 7, "tagpt": 8, "norm_z": 9, "norm_tt": 28, "FSR": 38}

    himpacts = OrderedDict()
    himpacts['total'] = ROOT.TH1F("htotal_wpt_{}".format(postfix), "htotal_wpt_{}".format(postfix), nbins, wptbins)
    for unc in sysUncs: 
        himpacts[unc] = ROOT.TH1F("h{}_wpt_{}".format(unc, postfix), "h{}_wpt_{}".format(unc, postfix), nbins, wptbins)

    for ibin, jname in enumerate(jsonNames, start=1):
        with open(jname) as f:
            results = json.load(f)
    
        # assuming the 0th POI is the signal strength
        poiname = results['POIs'][0]['name']
        diff = results['POIs'][0]['fit'][1] - results['POIs'][0]['fit'][0]
        hmu.SetBinContent(ibin, results['POIs'][0]['fit'][1])
        hmu.SetBinError(ibin, diff)
        himpacts['total'].SetBinContent(ibin, abs(diff))

        # loop over the systematics, find the ones in the sysUncs list
        for systematic in results['params']:
            if systematic["name"] in sysUncs:
                diff = systematic["impact_"+poiname]
                himpacts[systematic["name"]].SetBinContent(ibin, abs(diff))

    for key, val in himpacts.iteritems():
        val.SetLineColor(colors[key])

    hmu.SetLineColor(1)
    DrawHistos([hmu], ['W^{+}#rightarrow#mu^{+}#nu'], 0, 120, "W p_{T} [GeV]", 0.5, 1.5, "#sigma_{Obs}/#sigma_{MC}", "hsignal_strength_"+postfix, dology=False, legendPos=[0.92, 0.88, 0.70, 0.80])

    DrawHistos(himpacts.values(), himpacts.keys(), 0, 120, "W p_{T} [GeV]", 0, 0.20, "Uncertainty", "himpacts_wpt_"+postfix, dology=False, legendPos=[0.92, 0.88, 0.70, 0.40], drawashist=True)


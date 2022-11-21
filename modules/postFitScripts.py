"""
script to make postfit comparisons
"""
import ROOT
import os,sys,math
from collections import OrderedDict
import json
import re
import numpy as np

from CMSPLOTS.myFunction import DrawHistos
from modules.SampleManager import DrawConfig

ROOT.gROOT.SetBatch(True)

def MakePostPlot(ifilename: str, channel: str, bins: np.array, suffix: str, showpull: bool = False, is5TeV: bool = False, startbin: int = 1, outdir: str = "plots"):
    """
    compare the unrolled postfit of data and templates
    """
    ifile = ROOT.TFile(ifilename)
    
    horgdata = ifile.Get("obs")

    # get the list of histograms saved in the file
    hkeys = ifile.GetListOfKeys()
    hkeys = [hkey.GetName() for hkey in hkeys]
    hnames_sig = []
    hnames_qcd = []
    hnames_ewks = []
    hnames_ttbar = []
    for hkey in hkeys:
        if is5TeV and "13TeV" in hkey:
            continue
        if not is5TeV and "5TeV" in hkey:
            continue
        # w signal
        if bool(re.match(r"expproc_\w*sig_postfit$", hkey)):
            hnames_sig.append( hkey )
        # z signal
        elif bool(re.match(r"expproc_\w*sig_postfit$", hkey)):
            hnames_sig.append( hkey )
        # qcd
        elif bool(re.match(r"expproc_QCD_\w*postfit$", hkey)):
            hnames_qcd.append( hkey )
        # ewks: DY, taunu, diboson, etc and EWK (for Z's)
        elif bool(re.match(r"expproc_taunu\w*postfit$", hkey)):
            hnames_ewks.append( hkey )
        elif bool(re.match(r"expproc_zxx\w*postfit$", hkey)):
            hnames_ewks.append( hkey )
        elif bool(re.match(r"expproc_VV\w*postfit$", hkey)):
            hnames_ewks.append( hkey )
        elif bool(re.match(r"expproc_EWK\w*postfit$", hkey)):
            hnames_ewks.append( hkey )
        # ttbar
        elif bool(re.match(r"expproc_tt\w*postfit$", hkey)):
            hnames_ttbar.append( hkey )
    assert len(hnames_sig)>=1, "There should be at least one sig histogram in file: {}".format(ifilename)
    print("sig ", hnames_sig)
    print("qcd ", hnames_qcd)
    print("ewk ", hnames_ewks)
    print("ttbar ", hnames_ttbar)

    ## read the postfit plots from input file
    hexpsig = None
    hexpewk = None
    hexpqcd = None
    hexpttbar = None

    for hkey in hkeys:
        if hkey in hnames_sig:
            if hexpsig is None:
                hexpsig = ifile.Get(hkey)
            else:
                hexpsig.Add( ifile.Get(hkey) )

        if hkey in hnames_ewks:
            if hexpewk is None:
                hexpewk = ifile.Get(hkey)
            else:
                hexpewk.Add( ifile.Get(hkey) )
        
        if hkey in hnames_ttbar:
            if hexpttbar is None:
                hexpttbar = ifile.Get(hkey)
            else:
                hexpttbar.Add( ifile.Get(hkey) )
        
        if hkey in hnames_qcd:
            if hexpqcd is None:
                hexpqcd = ifile.Get(hkey)
            else:
                hexpqcd.Add( ifile.Get(hkey) )
        
    # the combined prediction of all processes,
    # which should have included the correct total postfit uncertainties
    hexpfull = ifile.Get("expfull_postfit")

    # the histograms saved in the root file does not follow the original bining
    # recover the original binning
    nbins = len(bins) - 1
    binnings = (nbins, bins)
    
    #binnings = (newbins.shape[0]-1, newbins)
    print("channel" , channel)
    print("suffix" , suffix)
    hdata   = ROOT.TH1D("hdata_{}_{}".format( channel, suffix),  "hdata_{}_{}".format( channel, suffix),  *binnings)
    hsig    = ROOT.TH1D("hsig_{}_{}".format(  channel, suffix),  "hsig_{}_{}".format(  channel, suffix),  *binnings)
    hewk    = ROOT.TH1D("hewk_{}_{}".format(  channel, suffix),  "hewk_{}_{}".format(  channel, suffix),  *binnings)
    httbar  = ROOT.TH1D("httbar_{}_{}".format(channel, suffix),  "httbar_{}_{}".format(channel, suffix),  *binnings)
    hqcd    = ROOT.TH1D("hqcd_{}_{}".format(  channel, suffix),  "hqcd_{}_{}".format(  channel, suffix),  *binnings)
    hratio  = ROOT.TH1D("hrato_{}_{}".format( channel, suffix),  "hratio_{}_{}".format(channel, suffix),  *binnings)
    hpull   = ROOT.TH1D("hpull_{}_{}".format( channel, suffix),  "hpull_{}_{}".format( channel, suffix),  *binnings)
    for ibin in range(1, nbins + 1):
        hdata.SetBinContent(ibin,   horgdata.GetBinContent(ibin + startbin-1))
        hdata.SetBinError(ibin,     horgdata.GetBinError(ibin + startbin-1 ))
        if hexpsig:
            hsig.SetBinContent(ibin,    hexpsig.GetBinContent(ibin + startbin-1))
        if hexpewk:
            hewk.SetBinContent(ibin,    hexpewk.GetBinContent(ibin + startbin-1))
        if hexpttbar:
            httbar.SetBinContent(ibin,  hexpttbar.GetBinContent(ibin + startbin-1))
        if hexpqcd:
            hqcd.SetBinContent(ibin,    hexpqcd.GetBinContent(ibin + startbin-1))
            #hqcd.SetBinError(ibin,      hexpqcd.GetBinError(ibin + startbin))

        hratio.SetBinContent(ibin, hexpfull.GetBinContent(ibin + startbin - 1))
        hratio.SetBinError(ibin,   hexpfull.GetBinError(ibin + startbin - 1))

        diff = horgdata.GetBinContent(ibin + startbin - 1) - hexpfull.GetBinContent(ibin + startbin - 1)
        # take the sigma as sqrt(data**2 + templates**2)
        # not 100% sure if this is the correct way to calculate pull
        sig = math.sqrt(horgdata.GetBinError(ibin + startbin - 1)**2 + hexpfull.GetBinError(ibin + startbin - 1)**2)
        hpull.SetBinContent(ibin, diff/(sig+1e-6))

    # deal with the uncertainty bar
    for ibin in range(1, hratio.GetNbinsX()+1):
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

    hs_gmc = ROOT.THStack("hs_stack_{}_{}".format(channel, suffix), "hs_stack")
    hs_gmc.Add(hqcd)
    hs_gmc.Add(httbar)
    hs_gmc.Add(hewk)
    hs_gmc.Add(hsig)

    ymaxs = {"muplus": 3.5e4, "muminus": 2.5e4, "eplus": 2.5e4, "eminus": 1.8e4, "mumu": 2.0e4, "ee": 1.5e4}
    if is5TeV:
        for key, val in ymaxs.items():
            ymaxs[key] = val*0.7

    siglabels = {"muplus":  "W^{+}#rightarrow#mu^{+}#nu", 
                 "muminus": "W^{-}#rightarrow#mu^{-}#bar{#nu}", 
                 "eplus":   "W^{+}#rightarrow e^{+}#nu", 
                 "eminus":  "W^{-}#rightarrow e^{-}#bar{#nu}",
                 "ee":      "Z#rightarrow e^{+}e^{-}",
                 "mumu":    "Z#rightarrow #mu^{+}#mu^{-}"}

    if "ee" not in channel and "mumu" not in channel:
        # w's
        xlabel = "m_{T} [GeV]"
        outputname = f"histo_wjets_{channel}_{suffix}"
    else:
        # z's
        xlabel = "m_{ll} [GeV]"
        outputname = f"histo_zjets_{channel}_{suffix}"
    drawconfigs = DrawConfig(xmin = bins.min(), xmax = bins.max(), xlabel = xlabel, ymin = 0, ymax = ymaxs[channel] / (int(nbins/36)+1), ylabel = "Events / GeV", outputname = outdir + "/" + outputname, dology=False, addOverflow=False, addUnderflow=False, yrmin=0.95, yrmax=1.05, yrlabel = "Data / Pred")
    DrawHistos( [hdata, hs_gmc], ["Data", siglabels[channel], "EWK", "t#bar{t}", "QCD"], drawconfigs.xmin, drawconfigs.xmax, drawconfigs.xlabel, drawconfigs.ymin, drawconfigs.ymax, drawconfigs.ylabel, drawconfigs.outputname, dology=drawconfigs.dology, dologx=drawconfigs.dologx, showratio=drawconfigs.showratio, yrmax = drawconfigs.yrmax, yrmin = drawconfigs.yrmin, yrlabel = drawconfigs.yrlabel, donormalize=drawconfigs.donormalize, ratiobase=drawconfigs.ratiobase, legendPos = drawconfigs.legendPos, redrawihist = drawconfigs.redrawihist, extraText = drawconfigs.extraText, noCMS = drawconfigs.noCMS, addOverflow = drawconfigs.addOverflow, addUnderflow = drawconfigs.addUnderflow, nMaxDigits = drawconfigs.nMaxDigits, hratiopanel=hratio, drawoptions=['PE', 'HIST same'], showpull=showpull, hpulls=[hpull], W_ref = 600 * int(nbins/36+1), is5TeV = is5TeV) 

    return nevts

def GetPOIValue(ifilename, poiname = ""):
    """
    return the POI val and error given a postfit root file
    """
    f = ROOT.TFile(ifilename)
    tree = f.Get("fitresults")
    tree.GetEntry(0)
    val = getattr(tree, poiname)
    err = abs(getattr(tree, poiname+"_err"))
    return val, err

def ComparePOIs(vals_x: np.array, vals: list, errs: list, labels: list, colors: list, markers: list, output: str, is5TeV: bool):
    """
    compare the POI values with different selections
    """
    print(vals_x)
    graphs = []
    for idx in range(len(vals)):
        val = vals[idx]
        err = errs[idx]
        color = colors[idx]
        marker = markers[idx]
        g = ROOT.TGraphErrors(len(vals_x), vals_x+1.0*idx, val, np.zeros(len(vals_x)), err)
        g.SetLineColor(color)
        g.SetMarkerColor(color)
        g.SetMarkerStyle(markers[idx])
        graphs.append(g)
    ymin = 0.97
    ymax = 1.03
    DrawHistos(graphs, labels, 15.0, 60, "m_{T} [GeV]", ymin, ymax, "POI", output, dology=False, showratio=False, donormalize=False, drawoptions='EP', legendPos = [0.2, 0.8, 0.8, 0.9], noCMS = True, nMaxDigits = 3, legendNCols = 2, is5TeV = is5TeV, legendoptions=["LEP", "LEP", "LEP"])


def result2json(ifilename: str, poiname: str, ofilename: str, hname: str = "nuisance_impact_mu"):
    """
    script to convert the postfit POI and impacts of nuisance parameters 
    to json file, which will be used to make impact plots later
    """
    nameMap = {
        "SysWeight1" : "eff_MC",
        "SysWeight2" : "eff_FSR",
        "SysWeight3" : "eff_bkg",
        "SysWeight4" : "eff_tagpt",
        "SysWeight6" : "Prefire",
        "SysWeight8": "PrefireECAL",
        "SysWeight10": "PrefireMuon",
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
        "ScaledMCshape": "QCD_SigContam",
    }

    def getNuisName(nuis):
        result = nuis
        for key, val in nameMap.items():
            if nuis.endswith(key):
                #result = val
                result = nuis.replace(key, val)
                break
        if bool(re.match(r"\w*bin\d+shape", nuis)):
            result = ("QCD_" + nuis).replace("shape", "")
        return result.replace("lepEta_bin0_WpT_bin0_", "")

    ifile = ROOT.TFile(ifilename)
    himpact = ifile.Get(hname)
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

    # sort impacts, descending
    impacts = OrderedDict(sorted(list(impacts.items()), key=lambda x: abs(x[1]), reverse=True))

    pulls = OrderedDict()
    for nuis in list(impacts.keys()):
        val = getattr(tree, nuis)
        err = getattr(tree, nuis+"_err")
        err = abs(err)
        pulls[nuis] = [val - err, val, val + err]

    # save to results
    for nuis in list(impacts.keys()):
        systematic = OrderedDict()
        systematic['fit'] = pulls[nuis]
        systematic['groups'] = []
        systematic['impact_' + poiname] = impacts[nuis]
        systematic['name'] = getNuisName(nuis)
        systematic['prefit'] = [-1.0, 0., 1.0]
        systematic[poiname] = [poi['fit'][1] - impacts[nuis], poi['fit'][1], poi['fit'][1] + impacts[nuis]]
        systematic['type'] = "Gaussian"
        print((getNuisName(nuis), pulls[nuis][1], pulls[nuis][1]-pulls[nuis][0], impacts[nuis]))

        results['params'].append(systematic)

    dirpath = ofilename.rpartition('/')[0]
    if not os.path.exists(dirpath):
        print(f"Make the directory {dirpath}")
        os.makedirs(dirpath)
    with open(ofilename, 'w') as fp:
        json.dump(results, fp, indent=2)
    
def DumpGroupImpacts(ifilename: str, poiname: str, hname = "nuisance_group_impact_mu"):
    """
    print out the grouped impacts
    """
    val_poi = GetPOIValue(ifilename, poiname)[0]
    print("val_poi ", val_poi)
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

    impacts = OrderedDict()
    for ibinY in range(1, himpact_grouped.GetNbinsY()+1):
        nuis = himpact_grouped.GetYaxis().GetBinLabel(ibinY)
        impacts[nuis] = himpact_grouped.GetBinContent(ibinX, ibinY) * 100.0 / val_poi


    # sort impacts, descending
    impacts = OrderedDict(sorted(list(impacts.items()), key=lambda x: abs(x[1]), reverse=True))

    print(f"\nPrint grouped nuisance impacts for {poiname} in {ifilename}")
    for nuis in list(impacts.keys()):
        print(f"{nuis:20}: {impacts[nuis]:.3f}")
    print()

    return impacts


def MakeWpTPostFitPlots(jsonNames, suffix=""):
    """
    plot the signal strength and systematic unc impacts as a function of wpt
    """
    import numpy as np
    # hard code the w pt binning for now
    wptbins = np.array([0., 8.0, 16.0, 24.0, 32.0, 40.0, 50.0, 70.0, 100.0, 120.0])
    nbins = wptbins.shape[0] - 1

    hmu = ROOT.TH1F("hmu_wpt_{}".format(suffix), "hmu_wpt_{}".format(suffix), nbins, wptbins)

    sysUncs = ['lumi_13TeV', 'Recoil', 'QCD', 'Prefire', 'norm_tt', 'effstat', 'norm_taunu', 'FSR', 'tagpt', 'norm_z']
    colors = {"total": 1, 'lumi_13TeV': 14, "Recoil": 2, "QCD": 3, "Prefire": 4, "effstat":6, "norm_taunu": 7, "tagpt": 8, "norm_z": 9, "norm_tt": 28, "FSR": 38}

    himpacts = OrderedDict()
    himpacts['total'] = ROOT.TH1F("htotal_wpt_{}".format(suffix), "htotal_wpt_{}".format(suffix), nbins, wptbins)
    for unc in sysUncs: 
        himpacts[unc] = ROOT.TH1F("h{}_wpt_{}".format(unc, suffix), "h{}_wpt_{}".format(unc, suffix), nbins, wptbins)

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

    for key, val in himpacts.items():
        val.SetLineColor(colors[key])

    hmu.SetLineColor(1)
    DrawHistos([hmu], ['W^{+}#rightarrow#mu^{+}#nu'], 0, 120, "W p_{T} [GeV]", 0.5, 1.5, "#sigma_{Obs}/#sigma_{MC}", "hsignal_strength_"+suffix, dology=False, legendPos=[0.92, 0.88, 0.70, 0.80])

    DrawHistos(list(himpacts.values()), list(himpacts.keys()), 0, 120, "W p_{T} [GeV]", 0, 0.20, "Uncertainty", "himpacts_wpt_"+suffix, dology=False, legendPos=[0.92, 0.88, 0.70, 0.40], drawashist=True)

def WriteOutputToText(output, ofilename):
    """
    write the output to a tex file
    """
    outdir = ofilename.rpartition('/')[0]
    if not os.path.exists(outdir):
        print(f"Make the directory {outdir}")
        os.makedirs(outdir)

    with open(ofilename, 'w') as f:
        f.write(output)

import ROOT
import os,sys
import numpy as np
from CMSPLOTS.myFunction import DrawHistos
import json
import math

ROOT.gROOT.SetBatch(True)

def ExpltOneBin(isocenters, bincontents, binerrors, isoSR, mTmin, mTmax, suffix="", extraText="", bincontents_scaled = None, binerrors_scaled = None, showScaled = True, isMuon = True, is5TeV = False):
    """
    extrapolate the QCD shape from a set of control regions (isocenters) to the signal region (isoSR),
    using linear extrapolation and the 2nd order polynomial function
    """
    fitmax = 0.85 if isMuon else 0.65
    fitmin = 0.15 if isMuon else 0.15
    graph = ROOT.TGraphErrors(len(bincontents), np.array(isocenters), np.array(bincontents), np.zeros(len(bincontents)), np.array(binerrors))
    f1 = ROOT.TF1("pol1_"+suffix, "[0]*(x-{}) + [1]".format(str(isoSR)), -0.1, fitmax)
    #f1 = ROOT.TF1("pol2_"+suffix, "[0]*(x-{isoSR})*(x-{isoSR}) + [2]*(x-{isoSR}) + [1]".format(isoSR=str(isoSR)), -0.1, fitmax)
    #f1 = ROOT.TF1("pol0_"+suffix, "[1]", -0.1, fitmax)

    graph.Fit(f1, "R", "", fitmin, fitmax)
    #graph.Fit(f2, "R", "", fitmin, fitmax)
    #print("val at Signal region", f1.Eval(isoSR))

    val_pol1_par1 = f1.GetParameter(1)
    err_pol1_par1 = f1.GetParError(1)
    val_pol1_par0 = f1.GetParameter(0)
    err_pol1_par0 = f1.GetParError(0)
    #val_pol2_par2 = f2.GetParameter(2)
    #err_pol2_par2 = f2.GetParError(2)
    #print("val ", val, " error ", err)

    f1.SetLineStyle(2)
    f1.SetLineColor(46)
    #f2.SetLineStyle(2)
    #f2.SetLineColor(9)

    graph2 = ROOT.TGraphErrors(1, np.array([isoSR]), np.array([val_pol1_par1]), np.zeros(1), np.array([abs(err_pol1_par1)]))
    graph2.SetMarkerColor(46)
    graph2.SetMarkerStyle(47)
    graph2.SetMarkerSize(2)

    #graph3 = ROOT.TGraphErrors(1, np.array([isoSR]), np.array([val_pol2_par2]), np.zeros(1), np.array([abs(err_pol2_par2)]))
    #graph3.SetMarkerColor(9)
    #graph3.SetMarkerStyle(45)
    #graph3.SetMarkerSize(2)

    #h_todraws = [graph, f1, f2, graph2, graph3]
    h_todraws = [graph, f1, graph2]
    #labels = ["{} < mT < {}".format(mTmin, mTmax), "Pol1 Fit", "Pol2 Fit", "Pol1 Extrapolation", "Pol2 Extrapolation"]
    labels = ["{} < mT < {}".format(mTmin, mTmax), "Pol1 Fit", "Pol1 Extrapolation"]
    #drawoptions = ["P same", "L", "L", "P same", "P same"]
    drawoptions = ["P same", "L", "P same"]
    #legendoptions=["EP", "L", "L", "EP", "EP"]
    legendoptions = ["EP", "L", "EP"]

    val_scaled_pol1_par1 = None
    err_scaled_pol1_par1 = None
    if bincontents_scaled:
        # repeat the fitting precedure on the scaled templates
        graph_scaled = ROOT.TGraphErrors(len(bincontents), np.array(isocenters), np.array(bincontents_scaled), np.zeros(len(bincontents)), np.array(binerrors_scaled))
        f3 = ROOT.TF1("pol1_scaled"+suffix, "[0]*(x-{}) + [1]".format(str(isoSR)), -0.1, 0.95)
        graph_scaled.Fit(f3, "R", "", fitmin, fitmax)
        val_scaled_pol1_par1 = f3.GetParameter(1)
        err_scaled_pol1_par1 = f3.GetParError(1)

        f3.SetLineStyle(2)
        f3.SetLineColor(30)

        graph_scaled.SetMarkerColor(12)
        graph_scaled.SetMarkerStyle(31)
        graph4 = ROOT.TGraphErrors(1, np.array([isoSR]), np.array([val_scaled_pol1_par1]), np.zeros(1), np.array([abs(err_scaled_pol1_par1)]))
        graph4.SetMarkerColor(30)
        graph4.SetMarkerStyle(41)
        graph4.SetMarkerSize(2)

        if showScaled:
            h_todraws.append(graph_scaled)
            h_todraws.append(f3)
            h_todraws.append(graph4)

            labels.append("Templates with Scaled MC")
            labels.append("Pol1 Fit with Scaled MC")
            labels.append("Pol1 Extrapolation with Scaled MC")

            drawoptions.extend(["P same", "L", "P same"])
            legendoptions.extend(["PE", "L", "PE"])

    DrawHistos( h_todraws, labels, 0, fitmax, "Lepton Relative Isolation", 0.4*min(bincontents), 1.25*max(bincontents), "Bin Content", "QCDBinContentNorm_"+suffix, dology=False, drawoptions=drawoptions, legendoptions=legendoptions, nMaxDigits=3, legendPos=[0.55, 0.18, 0.88, 0.50], lheader=extraText, is5TeV = is5TeV, noCMS = True)

    #return (val_pol1_par1, err_pol1_par1), (val_pol1_par0, err_pol1_par0), (val_pol2_par2, err_pol2_par2), (val_scaled_pol1_par1, err_scaled_pol1_par1)
    return (val_pol1_par1, err_pol1_par1), (val_pol1_par0, err_pol1_par0), (val_scaled_pol1_par1, err_scaled_pol1_par1)


def ExtrapolateQCD(fname, oname, channels, wptbin, etabins, fname_scaled=None, is5TeV=False, rebinned = False, qcdnorms = {}, doTest = False):
    """
    run the QCd extrapolation in all mT bins,
    save the statistical and systematic variations for HComb, and
    make some shape comparison plots
    """
    etaLabels = {
        "lepEta_bin0": "Full Eta", 
        "lepEta_bin1": "Barrel", 
        "lepEta_bin2": "Endcap"
    }
    channelLabels = {
        "muplus":  "W^{+}#rightarrow #mu^{+}#nu",
        "muminus": "W^{-}#rightarrow #mu^{-}#bar{#nu}",
        "eplus":   "W^{+}#rightarrow e^{+}#nu",
        "eminus":  "W^{-}#rightarrow e^{-}#bar{#nu}",
    }
    
    print(f"read QCD histograms from {fname}")

    fqcd = ROOT.TFile.Open(fname)
    if fname_scaled:
        # qcd templates with scaled MC subtraction
        fqcd_scaled = ROOT.TFile.Open(fname_scaled)

    sqrtS = "5TeV" if is5TeV else "13TeV"

    outdir = oname.rpartition('/')[0]
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    ofile = ROOT.TFile.Open(oname, "recreate")
  
    if isinstance(channels, str):
        channels = [channels]

    isMuon = 0

    for channel in channels:
        for etabin in etabins:
            # hard code some isolation parameters
            if channel.startswith("mu"):
                isMuon = 1
                # muon channel
                isomin = 4
                isomax = 17
                isocuts = [0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95]
                isocenters = [0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.525, 0.575, 0.625, 0.675, 0.725, 0.775, 0.825, 0.875, 0.925]
                #isoSR = 0.025 # average isolation value in the signal region
                isoSR = 0.084 # average isolation value in the signal region
                if doTest:
                    # test closure
                    isomin = 0
                    isocuts = [0.0, 0.0001, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95]
                    isocenters = [0., 0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.525, 0.575, 0.625, 0.675, 0.725, 0.775, 0.825, 0.875, 0.925]
                    
            else:
                isomin = 6
                isomax = 13
                #isocuts = [0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.50, 0.65, 0.90]
                isocuts = [0.15, 0.17, 0.20, 0.225, 0.25, 0.30, 0.35, 0.40, 0.50, 0.65, 0.90]

                #isocenters = [0.173, 0.223, 0.273, 0.323, 0.373, 0.441, 0.556, 0.724]
                isocenters = [0.160, 0.184, 0.212, 0.237, 0.273, 0.322, 0.373, 0.440, 0.558, 0.721]
                if etabin == "lepEta_bin0":
                    isoSR = 0.044
                elif etabin == "lepEta_bin1":
                    isoSR = 0.1007
                else:
                    isoSR = 0.0779

            # before applying linear extrapolation,
            # scale the histograms in the control region 
            # to the same normalization first
            histos_norm = {}
            histos_scaled_norm = {}
            for iso in range(isomin, isomax):
                prefix = "histo_wjetsAntiIso_mtcorr_weight_"
                if rebinned:
                    prefix = "Rebinned_" + prefix
                hname = prefix + channel + "_iso" + str(iso) + "_" + wptbin + "_" + etabin
                print("hname: ", hname)
                h = fqcd.Get(hname)
                # set the overflow and underflow to zero
                h.SetBinContent(0, 0)
                h.SetBinContent(h.GetNbinsX()+1, 0)
                h_norm =  h.Clone(hname+"_Cloned")
                histos_norm[iso] = h_norm

                if fname_scaled:
                    h = fqcd_scaled.Get(hname)
                    h.SetBinContent(0, 0)
                    h.SetBinContent(h.GetNbinsX()+1, 0)
                    h_norm = h.Clone(hname+"_Scaled_Cloned")
                    histos_scaled_norm[iso] = h_norm

            href = histos_norm[isomin]
            #counts = href.Integral()
            counts = 1.0
            #counts = href.GetBinContent(1)
            for iso in range(isomin, isomax):
                histos_norm[iso].Scale(counts / (histos_norm[iso].Integral() + 1e-6))
                #histos_norm[iso].Scale(counts / (histos_norm[iso].GetBinContent(1) + 1e-6))
                if fname_scaled:
                    histos_scaled_norm[iso].Scale(counts / (histos_scaled_norm[iso].Integral() + 1e-6))
                    #histos_scaled_norm[iso].Scale(counts / (histos_scaled_norm[iso].GetBinContent(1) + 1e-6))

            # some histograms to save the trend of function parameter variation as a function of mT
            # mostly for plotting purpose
            h_pol1_par0 = href.Clone("h_pol1_par0_{}_{}_{}".format(channel, etabin, wptbin)) # slope
            h_pol1_par1 = href.Clone("h_pol1_par1_{}_{}_{}".format(channel, etabin, wptbin)) # interception
            #h_pol2_par2 = href.Clone("h_pol2_par2_{}_{}_{}".format(channel, etabin, wptbin)) # interception
            h_scaled_pol1_par1 = href.Clone("h_scaled_pol1_par1_{}_{}_{}".format(channel, etabin, wptbin)) # interception for the scaled templates

            # save the extrapolated shape for HComb
            hnew = href.Clone("h_QCD_Extrapolated_" + channel + "_" + etabin + "_" + wptbin)
            #hnew_pol2 = href.Clone("h_QCD_Extrapolated_" + channel + "_" + etabin + "_" + wptbin + "_Pol2shapeUp")
            hnew_scaled = href.Clone("h_QCD_Extrapolated_" + channel + "_" + etabin + "_" + wptbin + "_ScaledMCshape_" + sqrtS + "Up")

            vals_pol1_par1 = []
            #
            # run the linear extrapolation bin-by-bin
            #
            for ibin in range(1, histos_norm[isomin].GetNbinsX()+1):
                bincontents = []
                binerrors = []
                for iso, hist in list(histos_norm.items()):
                    bincontents.append( hist.GetBinContent(ibin) )
                    binerrors.append( hist.GetBinError(ibin) )

                bincontents_scaled = None
                binerrors_scaled = None
                if fname_scaled:
                    bincontents_scaled = []
                    binerrors_scaled = []
                    for iso, hist in list(histos_scaled_norm.items()):
                        bincontents_scaled.append( hist.GetBinContent(ibin) )
                        binerrors_scaled.append( hist.GetBinError(ibin) )

                mTmin = histos_norm[isomin].GetBinLowEdge(ibin)
                mTmax = histos_norm[isomin].GetBinLowEdge(ibin) + histos_norm[isomin].GetBinWidth(ibin)
                suffix = channel+"_bin_"+str(ibin)+"_"+ etabin+"_"+wptbin + "_" + sqrtS
                #extraText = channelLabels[channel] +" "+ etaLabels[etabin]
                extraText = channelLabels[channel]
                #results_pol1_par1, results_pol1_par0, results_pol2_par2, results_scaled_pol1_par1 = ExpltOneBin(isocenters, bincontents, binerrors, isoSR, mTmin, mTmax, suffix=suffix, extraText=extraText, bincontents_scaled = bincontents_scaled, binerrors_scaled = binerrors_scaled)
                results_pol1_par1, results_pol1_par0, results_scaled_pol1_par1 = ExpltOneBin(isocenters, bincontents, binerrors, isoSR, mTmin, mTmax, suffix=suffix, extraText=extraText, bincontents_scaled = bincontents_scaled, binerrors_scaled = binerrors_scaled, showScaled=False, isMuon = isMuon, is5TeV=is5TeV)
                
                hnew.SetBinContent(ibin, max(results_pol1_par1[0], 0))
                hnew.SetBinError(ibin, 0.)
                #hnew_pol2.SetBinContent(ibin, max(results_pol2_par2[0], 0))
                #hnew_pol2.SetBinError(ibin, 0.)
                if fname_scaled:
                    hnew_scaled.SetBinContent(ibin, max(results_scaled_pol1_par1[0], 0))
                    hnew_scaled.SetBinError(ibin, 0.)

                h_pol1_par1.SetBinContent(ibin, results_pol1_par1[0])
                h_pol1_par1.SetBinError(ibin,   results_pol1_par1[1])
                h_pol1_par0.SetBinContent(ibin, results_pol1_par0[0])
                h_pol1_par0.SetBinError(ibin,   results_pol1_par0[1])
                #h_pol2_par2.SetBinContent(ibin, results_pol2_par2[0])
                #h_pol2_par2.SetBinError(ibin,   results_pol2_par2[1])
                if fname_scaled:
                    h_scaled_pol1_par1.SetBinContent(ibin, results_scaled_pol1_par1[0])
                    h_scaled_pol1_par1.SetBinError(ibin, results_scaled_pol1_par1[1])

                vals_pol1_par1.append(results_pol1_par1)

            # set the bin-by-bin shape variation (stat.) for HComb
            hnew_ups = []
            hnew_downs = []
            for ibin in range(1, histos_norm[isomin].GetNbinsX()+1):
                val = max(vals_pol1_par1[ibin-1][0], 0.)
                #val = href.GetBinContent(ibin)
                err = vals_pol1_par1[ibin-1][1]
                hnew_up   = hnew.Clone("h_QCD_Extrapolated_"+channel+"_"+etabin+"_"+wptbin+"_bin{}shape_{}Up".format(str(ibin), sqrtS))
                hnew_down = hnew.Clone("h_QCD_Extrapolated_"+channel+"_"+etabin+"_"+wptbin+"_bin{}shape_{}Down".format(str(ibin), sqrtS))
                hnew_up.SetBinContent(ibin, val+err)
                hnew_up.SetBinError(ibin, 0.)
                hnew_down.SetBinContent(ibin, max(val-err, 0.))
                hnew_down.SetBinError(ibin, 0.)

                hnew_ups.append(hnew_up)
                hnew_downs.append(hnew_down)

            # pol2 as another systematic
            #hnew_pol2.Scale(hnew.Integral() / hnew_pol2.Integral())
            #hnew_pol2Dn = hnew_pol2.Clone("h_QCD_Extrapolated_" + channel + "_" + etabin + "_" + wptbin + "_Pol2shapeDown")
            #for ibin in range(1, hnew.GetNbinsX()+1):
            #    hnew_pol2Dn.SetBinContent(ibin, 2*hnew.GetBinContent(ibin) - hnew_pol2.GetBinContent(ibin))
            #hnew_ups.append(hnew_pol2)
            #hnew_downs.append(hnew_pol2Dn)

            # scaled MC as another systematic
            if fname_scaled:
                #hnew_scaled.Scale(hnew.Integral() / (hnew_scaled.Integral() + 1e-6))
                hnew_scaledDn = hnew_scaled.Clone("h_QCD_Extrapolated_" + channel + "_" + etabin + "_" + wptbin + "_ScaledMCshape_" + sqrtS + "Down")
                for ibin in range(1, hnew.GetNbinsX()+1):
                    hnew_scaledDn.SetBinContent(ibin, 2*hnew.GetBinContent(ibin) - hnew_scaled.GetBinContent(ibin))
                hnew_ups.append(hnew_scaled)
                hnew_downs.append(hnew_scaledDn)

            h_pol1_par1.SetLineColor(46)
            h_pol1_par1.SetMarkerColor(46)
            #h_pol2_par2.Scale(h_pol1_par1.Integral() / h_pol2_par2.Integral())
            #h_pol2_par2.SetLineColor(9)
            #h_pol2_par2.SetMarkerColor(9)
            #h_todraws = [h_pol1_par1, h_pol2_par2]

            iso_to_show = 7
            isobin = iso_to_show - isomin
            h_todraws = [h_pol1_par1, histos_norm[iso_to_show]]
            #labels = ["Pol1 Extrapolation", "Pol2 Extrapolation"]
            labels = ["Extrapolated", f"{isocuts[isobin]} < I < {isocuts[isobin+1]}"]
            if fname_scaled:
                h_scaled_pol1_par1.SetLineColor(30)
                h_scaled_pol1_par1.SetMarkerColor(30)
                h_todraws.append(h_scaled_pol1_par1)
                labels.append("With Scaled MC")
            DrawHistos( h_todraws, labels, 0, 140, "m_{T} [GeV]", 0., 0.3, "A.U.", "QCDShapeCompare_"+channel+"_"+etabin+"_"+wptbin+ "_" + sqrtS, dology=False, nMaxDigits=3, legendPos=[0.60, 0.65, 0.88, 0.83], lheader=extraText, is5TeV = is5TeV, noCMS = True, donormalize = True, addOverflow = True, showratio = True, yrmin = 0.81, yrmax = 1.19)

            # 
            # draw the original (normalized) mT distribution in different anti-isolated regions
            #
            h_todraws = []
            labels = []
            suffix = channel+"_"+ etabin+"_"+wptbin + "_"+sqrtS
            extraText = channelLabels[channel] 
            colors = [1, 2, 3, 4, 6, 7, 8, 9, 28]
            icolor = 0
            for idx, (iso, hmt) in enumerate(histos_norm.items()):
                #if idx != 0 and idx % 3 == 0:
                #    # remove some histograms to make the plot more clean
                #    continue
                if idx >= 9:
                    continue
                hmt.SetLineColor(colors[icolor])
                icolor += 1
                h_todraws.append(hmt)

                labels.append(f"{isocuts[idx]} < I < {isocuts[idx+1]}")
            DrawHistos( h_todraws, labels, 0, 140, "m_{T} [GeV]", 0., 0.3, "A.U.", "mTShape_QCD_"+suffix, dology=False, nMaxDigits=3, legendPos=[0.65, 0.50, 0.88, 0.85], lheader=extraText, noCMS = True, donormalize=True, drawashist=True, showratio=True, yrmin=0.81, addOverflow = True, yrmax=1.39, is5TeV = is5TeV)

            norm = 1.0
            if qcdnorms:   
                norm = GetQCDNorm(qcdnorms, channel, is5TeV)

            #
            # write the variations to the output
            #
            hnew.Scale(norm)
            hnew.SetDirectory(ofile)
            hnew.Write()
            #for h in h_todraws + [h_pol1_par0]:
            #    h.SetDirectory(ofile)
            #    h.Write()
            for hnew_up in hnew_ups:
                hnew_up.Scale(norm)
                hnew_up.SetDirectory(ofile)
                hnew_up.Write()
            for hnew_down in hnew_downs:
                hnew_down.Scale(norm)
                hnew_down.SetDirectory(ofile)
                hnew_down.Write()

    ofile.Close()
    
def InclusiveQCD(fname, oname, channels, wptbin, etabins, fname_scaled=None, is5TeV=False, rebinned = False, qcdnorms = {}):
    """
    test using directly the template from the anti-isolated region
    """
    etaLabels = {
        "lepEta_bin0": "Full Eta", 
        "lepEta_bin1": "Barrel", 
        "lepEta_bin2": "Endcap"
    }
    channelLabels = {
        "muplus":  "W^{+}#rightarrow #mu^{+}#nu",
        "muminus": "W^{-}#rightarrow #mu^{-}#bar{#nu}",
        "eplus":   "W^{+}#rightarrow e^{+}#nu",
        "eminus":  "W^{-}#rightarrow e^{-}#bar{#nu}",
    }
    
    print(f"read QCD histograms from {fname}")

    fqcd = ROOT.TFile.Open(fname)

    sqrtS = "5TeV" if is5TeV else "13TeV"

    outdir = oname.rpartition('/')[0]
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    ofile = ROOT.TFile.Open(oname, "recreate")
  
    if isinstance(channels, str):
        channels = [channels]

    for channel in channels:
        for etabin in etabins:
            # hard code some isolation parameters
            if channel.startswith("mu"):
                isogroups = [[7,8,9,10,11], [12,13,14,15,16]]
            else:
                isogroups = [[6,7,8,9],[10,11,12,13]]

            # read the histograms
            histos_norm = {}
            isobin = 0
            for isogroup in isogroups:
                h_norm = None
                for iso in isogroup:
                    prefix = "histo_wjetsAntiIso_mtcorr_weight_"
                    if rebinned:
                        prefix = "Rebinned_" + prefix
                    hname = prefix + channel + "_iso" + str(iso) + "_" + wptbin + "_" + etabin
                    print("hname: ", hname)
                    h = fqcd.Get(hname)
                    if not h_norm:
                        h_norm = h.Clone(hname+"_Cloned")
                    else:
                        h_norm.Add(h)
                histos_norm[isobin] = h_norm
                isobin += 1

            href = histos_norm[0]
            counts = href.Integral()
            #counts = href.GetBinContent(1)
            for iso in histos_norm.keys():
                histos_norm[iso].Scale(counts / (histos_norm[iso].Integral() + 1e-6))
                #histos_norm[iso].Scale(counts / (histos_norm[iso].GetBinContent(1) + 1e-6))

            # save the extrapolated shape for HComb
            hnew = href.Clone("h_QCD_Extrapolated_" + channel + "_" + etabin + "_" + wptbin)

            # set the bin-by-bin shape variation (stat.) for HComb
            hnew_ups = []
            hnew_downs = []
            for ibin in range(1, histos_norm[0].GetNbinsX()+1):
                val = href.GetBinContent(ibin)
                err = math.sqrt((href.GetBinContent(ibin) - histos_norm[1].GetBinContent(ibin))**2 + href.GetBinError(ibin)**2)
                hnew_up   = hnew.Clone("h_QCD_Extrapolated_"+channel+"_"+etabin+"_"+wptbin+"_bin{}shape_{}Up".format(str(ibin), sqrtS))
                hnew_down = hnew.Clone("h_QCD_Extrapolated_"+channel+"_"+etabin+"_"+wptbin+"_bin{}shape_{}Down".format(str(ibin), sqrtS))
                hnew_up.SetBinContent(ibin, val+err)
                hnew_up.SetBinError(ibin, 0.)
                hnew_down.SetBinContent(ibin, max(val-err, 0.))
                hnew_down.SetBinError(ibin, 0.)

                hnew_ups.append(hnew_up)
                hnew_downs.append(hnew_down)

            norm = 1.0
            if qcdnorms:   
                norm = GetQCDNorm(qcdnorms, channel, is5TeV)

            #
            # write the variations to the output
            #
            hnew.Scale(norm)
            hnew.SetDirectory(ofile)
            hnew.Write()
            for hnew_up in hnew_ups:
                hnew_up.Scale(norm)
                hnew_up.SetDirectory(ofile)
                hnew_up.Write()
            for hnew_down in hnew_downs:
                hnew_down.Scale(norm)
                hnew_down.SetDirectory(ofile)
                hnew_down.Write()

    ofile.Close()

def LoadQCDNorms(fname):
    """
    read the QCD normalization from the json file
    """
    with open(fname, "r") as f:
        norms = json.load(f)
    return norms

def GetQCDNorm(norms, channel, is5TeV = False):
    """
    get the QCD normalization for a given channel
    """
    suffix = "5TeV" if is5TeV else "13TeV"
    return norms[f"QCD_{channel}_lepEta_bin0_WpT_bin0_{suffix}"]

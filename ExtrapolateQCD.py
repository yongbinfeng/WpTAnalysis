import ROOT
import os,sys
import numpy as np
from CMSPLOTS.myFunction import DrawHistos

doMuon = True
doElectron = True
doWpT = False

ROOT.gROOT.SetBatch(True)

etaLabels = {
    "lepEta_bin0": "Full Eta", 
    "lepEta_bin1": "Barrel", 
    "lepEta_bin2": "Endcap"
}
channelLabels = {
    "muplus":  "W^{+}#rightarrow #mu^{+}#nu",
    "muminus": "W^{-}#rightarrow #mu^{-}#nu",
    "eplus":   "W^{+}#rightarrow e^{+}#nu",
    "eminus":  "W^{-}#rightarrow e^{-}#nu",
}

if doWpT:
    wptbins = ["WpT_bin1", "WpT_bin2", "WpT_bin3", "WpT_bin4", "WpT_bin5", "WpT_bin6", "WpT_bin7", "WpT_bin8", "WpT_bin9"]
else:
    wptbins = ["WpT_bin0"]


def ExpltOneBin(isocenters, bincontents, binerrors, isoSR, mTmin, mTmax, suffix="", extraText="", bincontents_scaled = None, binerrors_scaled = None):
    """
    extrapolate the QCD shape from a set of control regions (isocenters) to the signal region (isoSR),
    using linear extrapolation and the 2nd order polynomial function
    """
    graph = ROOT.TGraphErrors(len(bincontents), np.array(isocenters), np.array(bincontents), np.zeros(len(bincontents)), np.array(binerrors))
    f1 = ROOT.TF1("pol1_"+suffix, "[0]*(x-{}) + [1]".format(str(isoSR)), -0.1, 0.60)
    f2 = ROOT.TF1("pol2_"+suffix, "[0]*(x-{isoSR})*(x-{isoSR}) + [1]*(x-{isoSR}) + [2]".format(isoSR=str(isoSR)), -0.1, 0.60)
    # fit range
    fitmin = 0.25
    fitmax = 0.60
    graph.Fit(f1, "R", "", fitmin, fitmax)
    graph.Fit(f2, "R", "", fitmin, fitmax)
    #print("val at Signal region", f1.Eval(isoSR))

    val_pol1_par1 = f1.GetParameter(1)
    err_pol1_par1 = f1.GetParError(1)
    val_pol1_par0 = f1.GetParameter(0)
    err_pol1_par0 = f1.GetParError(0)
    val_pol2_par2 = f2.GetParameter(2)
    err_pol2_par2 = f2.GetParError(2)
    #print("val ", val, " error ", err)

    f1.SetLineStyle(2)
    f1.SetLineColor(46)
    f2.SetLineStyle(2)
    f2.SetLineColor(9)

    graph2 = ROOT.TGraphErrors(1, np.array([isoSR]), np.array([val_pol1_par1]), np.zeros(1), np.array([abs(err_pol1_par1)]))
    graph2.SetMarkerColor(46)
    graph2.SetMarkerStyle(47)
    graph2.SetMarkerSize(2)

    graph3 = ROOT.TGraphErrors(1, np.array([isoSR]), np.array([val_pol2_par2]), np.zeros(1), np.array([abs(err_pol2_par2)]))
    graph3.SetMarkerColor(9)
    graph3.SetMarkerStyle(45)
    graph3.SetMarkerSize(2)

    h_todraws = [graph, f1, f2, graph2, graph3]
    labels = ["{} < mT < {}".format(mTmin, mTmax), "Pol1 Fit", "Pol2 Fit", "Pol1 Extrapolation", "Pol2 Extrapolation"]
    drawoptions = ["P same", "L", "L", "P same", "P same"]
    legendoptions=["EP", "L", "L", "EP", "EP"]

    val_scaled_pol1_par1 = None
    err_scaled_pol1_par1 = None
    if bincontents_scaled:
        # repeat the fitting precedure on the scaled templates
        graph_scaled = ROOT.TGraphErrors(len(bincontents), np.array(isocenters), np.array(bincontents_scaled), np.zeros(len(bincontents)), np.array(binerrors_scaled))
        f3 = ROOT.TF1("pol1_scaled"+suffix, "[0]*(x-{}) + [1]".format(str(isoSR)), -0.1, 0.60)
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

        h_todraws.append(graph_scaled)
        h_todraws.append(f3)
        h_todraws.append(graph4)

        labels.append("Templates with Scaled MC")
        labels.append("Pol1 Fit with Scaled MC")
        labels.append("Pol1 Extrapolation with Scaled MC")
        
        drawoptions.extend(["P same", "L", "P same"])
        legendoptions.extend(["PE", "L", "PE"])

    DrawHistos( h_todraws, labels, 0, 0.6, "Lepton Relative Isolation", 0.5*min(bincontents), 1.25*max(bincontents), "Bin Content", "QCDBinContentNorm_"+suffix, dology=False, drawoptions=drawoptions, legendoptions=legendoptions, nMaxDigits=3, legendPos=[0.65, 0.18, 0.88, 0.58], lheader=extraText)

    return (val_pol1_par1, err_pol1_par1), (val_pol1_par0, err_pol1_par0), (val_pol2_par2, err_pol2_par2), (val_scaled_pol1_par1, err_scaled_pol1_par1)


def ExtrapolateQCD(fname, oname, channel, wptbin, etabins, fname_scaled=None):
    """
    run the QCd extrapolation in all mT bins,
    save the statistical and systematic variations for HComb, and
    make some shape comparison plots
    """
    fqcd = ROOT.TFile(fname)
    if fname_scaled:
        # qcd templates with scaled MC subtraction
        fqcd_scaled = ROOT.TFile(fname_scaled)

    if not os.path.exists("root/QCD"):
        os.makedirs("root/QCD")
    ofile = ROOT.TFile("root/QCD/"+oname, "recreate")
    
    for etabin in etabins:
        # hard code some isolation parameters
        if channel.startswith("mu"):
            # muon channel
            isomin = 5
            isomax = 13
            isocuts = [0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65]
            isocenters = [0.225, 0.275, 0.320, 0.375, 0.425, 0.475, 0.525, 0.575]
            isoSR = 0.025 # average isolation value in the signal region
        else:
            isomin = 5
            isomax = 11
            isocuts = [0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.55, 0.70]
            isocenters = [0.2217, 0.2721, 0.3221, 0.372, 0.4219, 0.4902, 0.603]
            if etabin == "lepEta_bin0":
                isoSR = 0.09
            elif etabin == "lepEta_bin1":
                isoSR = 0.1007
            else:
                isoSR = 0.0779

        # before applying linear extrapolation,
        # scale the histograms in the control region 
        # to the same normalization first
        histos_norm = {}
        histos_scaled_norm = {}
        for iso in xrange(isomin, isomax):
            hname = wptbin + "/histo_wjetsAntiIso_fullrange_mtcorr_weight_" + channel + "_iso" + str(iso) + "_" + wptbin + "_" + etabin
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
        counts = href.Integral()
        for iso in xrange(isomin,isomax):
            histos_norm[iso].Scale(counts / histos_norm[iso].Integral())
            if fname_scaled:
                histos_scaled_norm[iso].Scale(counts / histos_scaled_norm[iso].Integral())

        # some histograms to save the trend of function parameter variation as a function of mT
        # mostly for plotting purpose
        h_pol1_par1 = href.Clone("h_pol1_par1_{}_{}_{}".format(channel, etabin, wptbin)) # interception
        h_pol1_par0 = href.Clone("h_pol1_par0_{}_{}_{}".format(channel, etabin, wptbin)) # slope
        h_pol2_par2 = href.Clone("h_pol2_par2_{}_{}_{}".format(channel, etabin, wptbin)) # interception
        h_scaled_pol1_par1 = href.Clone("h_scaled_pol1_par1_{}_{}_{}".format(channel, etabin, wptbin)) # interception for the scaled templates

        # save the extrapolated shape for HComb
        hnew = href.Clone("h_QCD_Extrapolated_" + channel + "_" + etabin + "_" + wptbin)
        hnew_pol2 = href.Clone("h_QCD_Extrapolated_" + channel + "_" + etabin + "_" + wptbin + "_Pol2shapeUp")
        hnew_scaled = href.Clone("h_QCD_Extrapolated_" + channel + "_" + etabin + "_" + wptbin + "_ScaledMCshapeUp")

        vals_pol1_par1 = []
        #
        # run the linear extrapolation bin-by-bin
        #
        for ibin in xrange(1, histos_norm[isomin].GetNbinsX()+1):
            bincontents = []
            binerrors = []
            for iso, hist in histos_norm.iteritems():
                bincontents.append( hist.GetBinContent(ibin) )
                binerrors.append( hist.GetBinError(ibin) )

            bincontents_scaled = None
            binerrors_scaled = None
            if fname_scaled:
                bincontents_scaled = []
                binerrors_scaled = []
                for iso, hist in histos_scaled_norm.iteritems():
                    bincontents_scaled.append( hist.GetBinContent(ibin) )
                    binerrors_scaled.append( hist.GetBinError(ibin) )

            mTmin = histos_norm[isomin].GetBinLowEdge(ibin)
            mTmax = histos_norm[isomin].GetBinLowEdge(ibin) + histos_norm[isomin].GetBinWidth(ibin)
            suffix = channel+"_bin_"+str(ibin)+"_"+ etabin+"_"+wptbin
            extraText = channelLabels[channel] +" "+ etaLabels[etabin]
            results_pol1_par1, results_pol1_par0, results_pol2_par2, results_scaled_pol1_par1 = ExpltOneBin(isocenters, bincontents, binerrors, isoSR, mTmin, mTmax, suffix=suffix, extraText=extraText, bincontents_scaled = bincontents_scaled, binerrors_scaled = binerrors_scaled)

            hnew.SetBinContent(ibin, max(results_pol1_par1[0], 0))
            hnew.SetBinError(ibin, 0.)
            hnew_pol2.SetBinContent(ibin, max(results_pol2_par2[0], 0))
            hnew_pol2.SetBinError(ibin, 0.)
            if fname_scaled:
                hnew_scaled.SetBinContent(ibin, max(results_scaled_pol1_par1[0], 0))
                hnew_scaled.SetBinError(ibin, 0.)

            h_pol1_par1.SetBinContent(ibin, results_pol1_par1[0])
            h_pol1_par1.SetBinError(ibin,   results_pol1_par1[1])
            h_pol1_par0.SetBinContent(ibin, results_pol1_par0[0])
            h_pol1_par0.SetBinError(ibin,   results_pol1_par0[1])
            h_pol2_par2.SetBinContent(ibin, results_pol2_par2[0])
            h_pol2_par2.SetBinError(ibin,   results_pol2_par2[1])
            if fname_scaled:
                h_scaled_pol1_par1.SetBinContent(ibin, results_scaled_pol1_par1[0])
                h_scaled_pol1_par1.SetBinError(ibin, results_scaled_pol1_par1[1])

            vals_pol1_par1.append(results_pol1_par1)

        # set the bin-by-bin shape variation (stat.) for HComb
        hnew_ups = []
        hnew_downs = []
        for ibin in xrange(1, histos_norm[isomin].GetNbinsX()+1):
            val = max(vals_pol1_par1[ibin-1][0], 0.)
            err = vals_pol1_par1[ibin-1][1]
            hnew_up   = hnew.Clone("h_QCD_Extrapolated_"+channel+"_"+etabin+"_"+wptbin+"_bin{}shapeUp".format(str(ibin)))
            hnew_down = hnew.Clone("h_QCD_Extrapolated_"+channel+"_"+etabin+"_"+wptbin+"_bin{}shapeDown".format(str(ibin)))
            hnew_up.SetBinContent(ibin, val+err)
            hnew_up.SetBinError(ibin, 0.)
            hnew_down.SetBinContent(ibin, max(val-err, 0.))
            hnew_down.SetBinError(ibin, 0.)

            hnew_ups.append(hnew_up)
            hnew_downs.append(hnew_down)

        # pol2 as another systematic
        hnew_pol2.Scale(hnew.Integral() / hnew_pol2.Integral())
        hnew_pol2Dn = hnew_pol2.Clone("h_QCD_Extrapolated_" + channel + "_" + etabin + "_" + wptbin + "_Pol2shapeDown")
        for ibin in xrange(1, hnew.GetNbinsX()+1):
            hnew_pol2Dn.SetBinContent(ibin, 2*hnew.GetBinContent(ibin) - hnew_pol2.GetBinContent(ibin))
        hnew_ups.append(hnew_pol2)
        hnew_downs.append(hnew_pol2Dn)

        # scaled MC as another systematic
        if fname_scaled:
            hnew_scaled.Scale(hnew.Integral() / hnew_scaled.Integral())
            hnew_scaledDn = hnew_scaled.Clone("h_QCD_Extrapolated_" + channel + "_" + etabin + "_" + wptbin + "_ScaledMCshapeDown")
            for ibin in xrange(1, hnew.GetNbinsX()+1):
                hnew_scaledDn.SetBinContent(ibin, 2*hnew.GetBinContent(ibin) - hnew_scaled.GetBinContent(ibin))
            hnew_ups.append(hnew_scaled)
            hnew_downs.append(hnew_scaledDn)

        h_pol1_par1.SetLineColor(46)
        h_pol1_par1.SetMarkerColor(46)
        h_pol2_par2.Scale(h_pol1_par1.Integral() / h_pol2_par2.Integral())
        h_pol2_par2.SetLineColor(9)
        h_pol2_par2.SetMarkerColor(9)
        h_todraws = [h_pol1_par1, h_pol2_par2]
        labels = ["Pol1 Extrapolation", "Pol2 Extrapolation"]
        if fname_scaled:
            h_scaled_pol1_par1.SetLineColor(30)
            h_scaled_pol1_par1.SetMarkerColor(30)
            h_todraws.append(h_scaled_pol1_par1)
            labels.append("Scaled MC")
        DrawHistos( h_todraws, labels, 0, 120, "m_{T} [GeV]", 0., 1.25*h_pol1_par1.GetMaximum(), "A.U.", "QCDShapeCompare_"+channel+"_"+etabin+"_"+wptbin, dology=False, nMaxDigits=3, legendPos=[0.60, 0.72, 0.88, 0.88], lheader=extraText)

        #
        # write the variations to the output
        #
        hnew.SetDirectory(ofile)
        hnew.Write()
        for h in h_todraws + [h_pol1_par0]:
            h.SetDirectory(ofile)
            h.Write()
        for hnew_up in hnew_ups:
            hnew_up.SetDirectory(ofile)
            hnew_up.Write()
        for hnew_down in hnew_downs:
            hnew_down.SetDirectory(ofile)
            hnew_down.Write()

    ofile.Close()


if __name__ == "__main__":
    if doMuon:
        fname = "root/output_qcdshape_fullrange_munu.root"
        oname = "qcdshape_extrapolated_mu"
        for wptbin in wptbins:
            ExtrapolateQCD(fname, oname+"_"+wptbin+".root", "muplus", wptbin, ["lepEta_bin0"], fname_scaled = fname.replace(".root", "_applyScaling.root"))
    if doElectron:
        # on electron channel
        # currently separate barrel and endcap for electrons
        fname = "root/output_qcdshape_fullrange_enu.root"
        oname = "qcdshape_extrapolated_e"
        for wptbin in ["WpT_bin0"]:
            # for electrons, only implement the inclusive version for now
            ExtrapolateQCD(fname, oname+"_"+wptbin+".root", "eplus", wptbin, ["lepEta_bin1", "lepEta_bin2"], fname_scaled = fname.replace(".root", "_applyScaling.root"))


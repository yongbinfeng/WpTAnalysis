import ROOT
import os,sys
import numpy as np
from CMSPLOTS.myFunction import DrawHistos

doPol2 = False # use Pol1 if doPol2 set to false
doMuon = False
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
    wptbins = ["WpT_bin0", "WpT_bin1", "WpT_bin2", "WpT_bin3", "WpT_bin4", "WpT_bin5", "WpT_bin6", "WpT_bin7", "WpT_bin8", "WpT_bin9"]
else:
    wptbins = ["WpT_bin0"]

def ExtrapolateQCD(fname, oname, channel, wptbin, etabins):
    fqcd = ROOT.TFile(fname)
    if not os.path.exists("QCD"):
        os.makedirs("QCD")
    ofile = ROOT.TFile("QCD/"+oname, "recreate")
    
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
        for iso in xrange(isomin, isomax):
            #hname = "histo_wjetsAntiIso_fullrange_mtcorr_weight_" + channel + "_iso" + str(iso) + "_" + wptbin + "_" + etabin
            # for QCD bkg template try not binning in W pT for now. test if it works
            hname = "histo_wjetsAntiIso_fullrange_mtcorr_weight_" + channel + "_iso" + str(iso) + "_" + "WpT_bin0" + "_" + etabin
            h = fqcd.Get(hname)
            # set the overflow and underflow to zero
            h.SetBinContent(0, 0)
            h.SetBinContent(h.GetNbinsX()+1, 0)
            h_norm =  h.Clone(hname+"_Cloned")
            histos_norm[iso] = h_norm

        href = histos_norm[isomin]
        counts = href.Integral()
        for iso in xrange(isomin,isomax):
            histos_norm[iso].Scale(counts / histos_norm[iso].Integral())

        # define two histograms to save the trend of p0 and p1 as a function of mT
        h_p0 = href.Clone("h_p0_{}_{}_{}".format(channel, etabin, wptbin))
        h_p1 = href.Clone("h_p1_{}_{}_{}".format(channel, etabin, wptbin))

        # save the extrapolated shape
        hnew = href.Clone("h_QCD_Extrapolated_" + channel + "_" + etabin + "_" + wptbin)

        vals = []
        uncs = []
        # run the linear extrapolation bin-by-bin
        for ibin in xrange(1, histos_norm[isomin].GetNbinsX()+1):
            bincontents = []
            binerrors = []
            for iso, hist in histos_norm.iteritems():
                bincontents.append( hist.GetBinContent(ibin) )
                binerrors.append( hist.GetBinError(ibin) )

            extraText = channelLabels[channel] +" "+ etaLabels[etabin]

            graph = ROOT.TGraphErrors(len(bincontents), np.array(isocenters), np.array(bincontents), np.zeros(len(bincontents)), np.array(binerrors))
            if not doPol2:
                f1 = ROOT.TF1("pol1_"+channel+"_iso"+str(iso)+"_"+ etabin+"_"+wptbin, "[0]*(x-{}) + [1]".format(str(isoSR)), -0.1, 0.60)
            else:
                f1 = ROOT.TF1("pol2_"+channel+"_iso"+str(iso)+"_"+ etabin+"_"+wptbin, "[0]*(x-{isoSR})*(x-{isoSR}) + [1]*(x-{isoSR}) + [2]".format(isoSR=str(isoSR)), -0.1, 0.60)
            # fit range
            graph.Fit(f1, "R", "", 0.21, 0.60)
            #print("val at Signal region", f1.Eval(isoSR))

            if not doPol2:
                val = f1.GetParameter(1)
                err = f1.GetParError(1)
            else:   
                val = f1.GetParameter(2)
                err = f1.GetParError(2)
            #print("val ", val, " error ", err)
            val = max(val, 0.)
            hnew.SetBinContent(ibin, val)
            hnew.SetBinError(ibin, 0.)

            h_p0.SetBinContent(ibin, val)
            h_p0.SetBinError(ibin,   err)
            h_p1.SetBinContent(ibin, f1.GetParameter(0))
            h_p1.SetBinError(ibin,   f1.GetParError(0))

            vals.append(val)
            uncs.append(err)

            f1.SetLineStyle(2)
            f1.SetLineColor(2)

            graph2 = ROOT.TGraphErrors(1, np.array([isoSR]), np.array([val]), np.zeros(1), np.array([abs(err)]))
            graph2.SetMarkerColor(3)
            graph2.SetMarkerStyle(47)
            graph2.SetMarkerSize(2)

            label = "Pol2 Fit" if doPol2 else "Pol1 Fit"
            mTmin = histos_norm[isomin].GetBinLowEdge(ibin)
            mTmax = histos_norm[isomin].GetBinLowEdge(ibin) + histos_norm[isomin].GetBinWidth(ibin)
            DrawHistos( [graph, f1, graph2], ["{} < mT < {}".format(mTmin, mTmax), label, "Extrapolation"], 0, 0.6, "Lepton Relative Isolation", 0.7*min(bincontents), 1.25*max(bincontents), "Bin Content", "QCDBinContentNorm_"+channel+"_bin_"+str(ibin)+"_"+etabin+"_"+wptbin, dology=False, drawoptions=["P same", "L", "P same"], legendoptions=["P", "L", "P"], nMaxDigits=3, legendPos=[0.65, 0.18, 0.88, 0.48], lheader=extraText)

        # set the bin-by-bin shape variation
        hnew_ups = []
        hnew_downs = []
        for ibin in xrange(1, histos_norm[isomin].GetNbinsX()+1):
            val = vals[ibin-1]
            err = abs(uncs[ibin-1])
            hnew_up   = hnew.Clone("h_QCD_Extrapolated_"+channel+"_"+etabin+"_"+wptbin+"_bin{}shapeUp".format(str(ibin)))
            hnew_down = hnew.Clone("h_QCD_Extrapolated_"+channel+"_"+etabin+"_"+wptbin+"_bin{}shapeDown".format(str(ibin)))
            hnew_up.SetBinContent(ibin, val+err)
            hnew_up.SetBinError(ibin, 0.)
            hnew_down.SetBinContent(ibin, max(val-err, 0.))
            hnew_down.SetBinError(ibin, 0.)

            hnew_ups.append(hnew_up)
            hnew_downs.append(hnew_down)

        #
        # write the variations to the output
        #
        hnew.SetDirectory(ofile)
        hnew.Write()
        h_p1.SetDirectory(ofile)
        h_p1.Write()
        h_p0.SetDirectory(ofile)
        h_p0.Write()
        for hnew_up in hnew_ups:
            hnew_up.SetDirectory(ofile)
            hnew_up.Write()
        for hnew_down in hnew_downs:
            hnew_down.SetDirectory(ofile)
            hnew_down.Write()

    ofile.Close()


if __name__ == "__main__":
    if doMuon:
        fname = "output_qcdshape_fullrange_munu.root"
        oname = "pe_extrapolated_mu_pol2" if doPol2 else "qcdshape_extrapolated_mu"
        for wptbin in wptbins:
            ExtrapolateQCD(fname, oname+"_"+wptbin+".root", "muplus", wptbin, ["lepEta_bin0"])
    else:
        # on electron channel
        # currently separate barrel and endcap for electrons
        fname = "output_qcdshape_fullrange_enu.root"
        oname = "pe_extrapolated_e_pol2" if doPol2 else "qcdshape_extrapolated_e"
        for wptbin in ["WpT_bin0"]:
            # for electrons, only implement the inclusive version for now
            ExtrapolateQCD(fname, oname+"_"+wptbin+".root", "eplus", wptbin, ["lepEta_bin1", "lepEta_bin2"])


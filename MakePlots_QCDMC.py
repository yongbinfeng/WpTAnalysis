"""
make some plots on the QCD MC samples
in order to compare with the extrapolated version
"""

import ROOT
import numpy as np
from collections import OrderedDict
import sys
from CMSPLOTS.myFunction import DrawHistos, THStack2TH1, AddOverflowsTH1, RebinHisto
from modules.histProcessor import ProcessHists, DoRebin
from modules.qcdExtrapolater import ExtrapolateQCD, LoadQCDNorms
from modules.Binnings import mass_bins_test
import CMSPLOTS.CMS_lumi
import pickle
from modules.SampleManager import DrawConfig, Sample, SampleManager
import argparse

ROOT.gROOT.SetBatch(True)
ROOT.ROOT.EnableImplicitMT(10)

def main():
    doReading = True
    doFit = False
    if doReading:
        input_data = "inputs/QCDMC/inputs_qcdmc.txt"
    
        DataSamp = Sample(input_data, isMC=False, legend="Data", name="Data", doTheoryVariation=False, isZSR=False)
    
        sampMan = SampleManager(DataSamp)
        #sampMan.ApplyCutAll("PV_npvs < 15")
        sampMan.ApplyCutAll("nMuon > 0 && Muon_pt[0] > 25 && Muon_tightId[0] == 1 && Muon_eta[0] > -2.4 && Muon_eta[0] < 2.4")
    
        outdir = "plots/QCDMC/"
        sampMan.outdir = outdir
    
        # copy the same as MakePlotsIso_Wlnu.py
        # such that the qcd extrapolation script can be used directly
        lepname = "mu"
        sampMan.DefineAll("q", "Muon_charge[0]")
        sampMan.DefineAll(lepname+"plus",  "q > 0")
        sampMan.DefineAll(lepname+"minus", "q < 0")
        chgbins = [lepname+"plus", lepname + "minus"]
    
        etabins = ["lepEta_bin0"]
        sampMan.DefineAll("lepEta_bin0", "1.0")
    
        wptbins = ["WpT_bin0"]
        sampMan.DefineAll("WpT_bin0",  "1.0")
    
        sampMan.DefineAll("relIso", "Muon_pfRelIso04_all[0]")
        sampMan.DefineAll("mtCorr", "sqrt(2*Muon_pt[0]*MET_pt*(1-cos(MET_phi-Muon_phi[0])))")
        sampMan.DefineAll("Lep_pt", "Muon_pt[0]")
        sampMan.DefineAll("Lep_eta", "Muon_eta[0]")
        sampMan.DefineAll("met_pt", "MET_pt")
    
        isoCuts = [0.0, 0.001, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50,
                       0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95]
        isobins = []
        sampMan.DefineAll("RelIso", "relIso")
        for isobin in range(0, 20):
            sampMan.DefineAll(
                f"w_iso{isobin}", f"(relIso >= {isoCuts[isobin]} && relIso < {isoCuts[isobin+1]})")
            isobins.append(f"iso{isobin}")
            
        # define a few signal and anti-isolated regions for comparisons
        sampMan.DefineAll("w_isosig", "relIso < 0.15")
        sampMan.DefineAll("w_isoAnti1", "relIso > 0.20 && relIso < 0.60")
        sampMan.DefineAll("w_isoAnti2", "relIso > 0.30 && relIso < 0.40")
        sampMan.DefineAll("w_isoAnti3", "relIso > 0.40 && relIso < 0.50")
        
        isobins.append("isosig")
        isobins.append("isoAnti1")
        isobins.append("isoAnti2")
        isobins.append("isoAnti3")

        sampMan.cacheDraw("RelIso", "histo_qcd_RelIso", 75, 0, 0.75, DrawConfig(xmin=0, xmax=0.75, xlabel="Relative Isolation", ylabel="Events / 0.01", dology=True, ymax=1e5, donormalizebin=False, addOverflow=True, addUnderflow=True, showratio=False, legendPos=[0.94, 0.88, 0.70, 0.68]))
    
        nbins = 12
        xmin = 0
        xmax = 140
        mass_bins = np.array([0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0,
                             70.0, 75.0, 80.0, 85.0, 90.0, 95.0, 100.0, 105.0, 110.0, 115.0, 120.0, 125.0, 130.0, 135.0, 140.])
    
        for wpt in wptbins:
            for lepeta in etabins:
                for chg in chgbins:
                    idx = 0
                    for iso in isobins:
                        strname = "weight_{}_{}_{}_{}".format(
                            chg, iso, wpt, lepeta)

                        sampMan.DefineAll(
                            strname, "w_{} * {} * {} * {}".format(iso, wpt, lepeta, chg))
                        
                        if iso == "isosig":
                            lheader = "I < 0.15"
                        elif iso == "isoAnti1":
                            lheader = "0.20 < I < 0.60"
                        elif iso == "isoAnti2":
                            lheader = "0.30 < I < 0.40"
                        elif iso == "isoAnti3":
                            lheader = "0.40 < I < 0.50"
                        else:
                            lheader = f"{isoCuts[idx]} < I < {isoCuts[idx + 1]}"
                        idx += 1

                        outputname = "histo_wjetsAntiIso_mtcorr_" + strname 
                        sampMan.cacheDraw("mtCorr", outputname, mass_bins, DrawConfig(xmin=xmin, xmax=xmax, xlabel="m_{T} [GeV]", ylabel=f"Events / {(xmax-xmin)/nbins:.0f} GeV", dology=True, ymax=6e5, donormalizebin=False, addOverflow=True, addUnderflow=True, showratio=False, lheader=lheader, legendPos=[0.94, 0.88, 0.70, 0.68]), weightname=strname)
                        
                        outputname = "histo_wjetsAntiIso_Lep_pt_" + strname
                        sampMan.cacheDraw("Lep_pt", outputname, 50, 20, 70, DrawConfig(xmin=20, xmax=70, xlabel="Lepton p_{T} [GeV]", ylabel="Events / 1 GeV", dology=True, ymax=1e5, donormalizebin=False, addOverflow=True, addUnderflow=True, showratio=False, lheader=lheader, legendPos=[0.94, 0.88, 0.70, 0.68]), weightname=strname)
                        
                        outputname = "histo_wjetsAntiIso_Lep_eta_" + strname
                        sampMan.cacheDraw("Lep_eta", outputname, 50, -2.5, 2.5, DrawConfig(xmin=-2.5, xmax=2.5, xlabel="Lepton #eta", ylabel="Events / 0.1", dology=True, ymax=1e5, donormalizebin=False, addOverflow=True, addUnderflow=True, showratio=False, lheader=lheader, legendPos=[0.94, 0.88, 0.70, 0.68]), weightname=strname)
                        
                        outputname = "histo_wjetsAntiIso_met_pt_" + strname
                        sampMan.cacheDraw("met_pt", outputname, 50, 0, 100, DrawConfig(xmin=0, xmax=100, xlabel="MET [GeV]", ylabel="Events / 2 GeV", dology=True, ymax=1e5, donormalizebin=False, addOverflow=True, addUnderflow=True, showratio=False, lheader=lheader, legendPos=[0.94, 0.88, 0.70, 0.68]), weightname=strname)
    
        sampMan.launchDraw()
        
        postfix = lepname + "nu.root"
        outfile =  "root/output_qcdMCshape_"+postfix
        outfile = ROOT.TFile.Open(outfile, "recreate")

        for wpt in wptbins:
            for lepeta in etabins:
                for chg in chgbins:
                    histos_mt = OrderedDict()
                    histos_leppt = OrderedDict()
                    histos_lepeta = OrderedDict()
                    histos_metpt = OrderedDict()
                    idx = 1
                    for iso in isobins:
                        strname = "weight_{}_{}_{}_{}".format(
                            chg, iso, wpt, lepeta)

                        # for mT
                        outputname = "histo_wjetsAntiIso_mtcorr_" + strname
                        hdata = sampMan.hdatas[outputname]
                        
                        if "isoAnti" in iso or "isosig" in iso:
                            histos_mt[iso] = hdata.Clone("histo_qcdAntiIso_mtcorr_" + strname + "_clone")
                            histos_mt[iso] = DoRebin(histos_mt[iso], mass_bins_test[0])
                            histos_mt[iso].SetMarkerColor(idx)
                            histos_mt[iso].SetLineColor(idx)
                        hdata.SetName(outputname)
                        hdata.SetDirectory(outfile)

                        hdata.Write()
                        
                        outputname = "histo_wjetsAntiIso_Lep_pt_" + strname
                        hdata = sampMan.hdatas[outputname]
                        if "isoAnti" in iso or "isosig" in iso:
                            histos_leppt[iso] = hdata.Clone("histo_qcdAntiIso_Lep_pt_" + strname + "_clone")
                            histos_leppt[iso].SetMarkerColor(idx)
                            histos_leppt[iso].SetLineColor(idx)
                            
                        outputname = "histo_wjetsAntiIso_Lep_eta_" + strname
                        hdata = sampMan.hdatas[outputname]
                        if "isoAnti" in iso or "isosig" in iso:
                            histos_lepeta[iso] = hdata.Clone("histo_qcdAntiIso_Lep_eta_" + strname + "_clone")
                            histos_lepeta[iso].SetMarkerColor(idx)
                            histos_lepeta[iso].SetLineColor(idx)
                        
                        outputname = "histo_wjetsAntiIso_met_pt_" + strname
                        hdata = sampMan.hdatas[outputname]
                        if "isoAnti" in iso or "isosig" in iso:
                            histos_metpt[iso] = hdata.Clone("histo_qcdAntiIso_met_pt_" + strname + "_clone")
                            histos_metpt[iso].SetMarkerColor(idx)
                            histos_metpt[iso].SetLineColor(idx)
                            
                            idx += 1
                    
                    channelLabels = {
                        "muplus":  "W^{+}#rightarrow #mu^{+}#nu",
                        "muminus": "W^{-}#rightarrow #mu^{-}#bar{#nu}",
                        "eplus":   "W^{+}#rightarrow e^{+}#nu",
                        "eminus":  "W^{-}#rightarrow e^{-}#bar{#nu}",
                    }
                    
                    labels = {}
                    labels["isosig"] = "Iso < 0.15"
                    labels["isoAnti1"] = "0.20 < Iso < 0.60"
                    labels["isoAnti2"] = "0.30 < Iso < 0.40"
                    labels["isoAnti3"] = "0.40 < Iso < 0.50"
                    
                    legends = [labels[isobin] for isobin in list(histos_mt.keys())]
                        
                    DrawHistos(list(histos_mt.values()), legends, 0, 140, "m_{T} [GeV]", 0., 0.3, "A.U.", "QCDShapeCompare_mT_" + chg, dology=False, nMaxDigits=3, legendPos=[0.60, 0.60, 0.88, 0.85], lheader=channelLabels[chg], noCMS = True, donormalize = True, addOverflow = True, showratio = True, yrmin = 0.96, yrmax = 1.04, outdir=outdir)
                    
                    DrawHistos(list(histos_leppt.values()), legends, 20, 70, "Lepton p_{T} [GeV]", 0., 0.21, "A.U.", "QCDShapeCompare_lepPt_" + chg, dology=False, nMaxDigits=3, legendPos=[0.60, 0.60, 0.88, 0.85], lheader=channelLabels[chg], noCMS = True, donormalize = True, addOverflow = True, showratio = True, yrmin = 0.7, yrmax = 1.3, outdir=outdir)
                    
                    DrawHistos(list(histos_lepeta.values()), legends, -2.5, 2.5, "Lepton #eta", 0., 0.06, "A.U.", "QCDShapeCompare_lepEta_" + chg, dology=False, nMaxDigits=3, legendPos=[0.60, 0.60, 0.88, 0.85], lheader=channelLabels[chg], noCMS = True, donormalize = True, addOverflow = True, showratio = True, yrmin = 0.5, yrmax = 1.5, outdir=outdir)
                    
                    DrawHistos(list(histos_metpt.values()), legends, 0, 100, "MET [GeV]", 0., 0.1, "A.U.", "QCDShapeCompare_metpt_" + chg, dology=False, nMaxDigits=3, legendPos=[0.60, 0.60, 0.88, 0.85], lheader=channelLabels[chg], noCMS = True, donormalize = True, addOverflow = True, showratio = True, yrmin = 0.81, yrmax = 1.19, outdir=outdir)
    
    
        outfile.Write()
        outfile.Close()

        sampMan.dumpCounts()
        
    if doFit:
        #
        # rebin and extrapolate qcd template
        #
        lepname = "mu"
        postfix = lepname + "nu.root"
        outfile =  "root/output_qcdMCshape_"+postfix
        fqcd_rebin = "root/output_qcdMCshape_rebin"+postfix
        fqcd_output = "root/output_qcdMCshape_extrap"+postfix
        includeUnderflow = True
        includeOverflow = True
        ProcessHists(outfile, fqcd_rebin, mass_bins_test[0], includeUnderflow, includeOverflow)

        # extrapolate the QCD template from anti-isolated region to isolated region
        # manually set the QCD normalization factors for the prefit impacts
        # very hacky, not ideal at all; only temporary solution
        fqcdnorm = "data/QCDNorms_inc.json"
        qcdnorms = LoadQCDNorms(fqcdnorm)
        print("QCD Norms: ", qcdnorms)
        # qcdnorms = {}
        ExtrapolateQCD(fqcd_rebin, fqcd_output, [lepname+"plus", lepname+"minus"], "WpT_bin0", [
                       "lepEta_bin0"], fname_scaled=None, rebinned=True, is5TeV=False, qcdnorms=qcdnorms, doTest=True)

    print("Program end...")
   
    
if __name__ == "__main__":
    main()

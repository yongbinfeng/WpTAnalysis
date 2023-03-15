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
    doGenMET = True
    
    if doReading:
        input_data = "inputs/QCDMC/inputs_qcdmc.txt"
    
        DataSamp = Sample(input_data, isMC=False, legend="Data", name="Data", doTheoryVariation=False, isZSR=False)
    
        sampMan = SampleManager(DataSamp)
        #sampMan.ApplyCutAll("PV_npvs < 15")
        sampMan.ApplyCutAll("nMuon > 0 && Muon_pt[0] > 25 && Muon_tightId[0] == 1 && Muon_eta[0] > -2.4 && Muon_eta[0] < 2.4")
    
        outdir = "plots/QCDMC/"
        sampMan.outdir = outdir
    
        # copy the same as MakePlotsIso_Wlnu.py
        # such that the qcd extrapolation and FR scripts can be used directly
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
        if doGenMET:
            sampMan.DefineAll("met_pt", "GenMET_pt")
            sampMan.DefineAll("met_phi", "GenMET_phi")
        else:
            sampMan.DefineAll("met_pt", "MET_pt")
            sampMan.DefineAll("met_phi", "MET_phi")
        sampMan.DefineAll("mtCorr", "sqrt(2*Muon_pt[0]*met_pt*(1-cos(met_phi-Muon_phi[0])))")
        sampMan.DefineAll("Lep_pt", "Muon_pt[0]")
        sampMan.DefineAll("Lep_eta", "Muon_eta[0]")
        sampMan.DefineAll("deltaPhi", "abs(TVector2::Phi_mpi_pi(met_phi - Muon_phi[0]))")
    
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
        
        mass_bins = mass_bins_test[0]
        mtbins = []
        for imt in range(len(mass_bins)-1):
            sampMan.DefineAll(
                f"mt{imt}", f"mtCorr >= {mass_bins[imt]} && mtCorr < {mass_bins[imt+1]}")
            mtbins.append(f"mt{imt}")
            
        def makeDraws(strname, ymax):
            sampMan.cacheDraw("RelIso", f"histo_wjets_{lepname}_RelIso_{strname}", 72, 0, 0.72, DrawConfig(xmin=0., xmax=0.75, xlabel="Relative Isolation", ylabel=f"Events / 0.01", dology=True, ymax=ymax, donormalizebin=False, addOverflow=True, addUnderflow=True, showratio=False, legendPos=[0.94, 0.88, 0.70, 0.68]), weightname=strname)
        
            sampMan.cacheDraw("Lep_pt", f"histo_wjets_{lepname}_Lep_pt_{strname}", 50, 20, 70, DrawConfig(xmin=20, xmax=70, xlabel="Lepton p_{T} [GeV]", ylabel=f"Events / GeV", dology=True, ymax=ymax, donormalizebin=False, addOverflow=True, addUnderflow=True, showratio=False, legendPos=[0.94, 0.88, 0.70, 0.68]), weightname=strname)
        
            sampMan.cacheDraw("Lep_eta", f"histo_wjets_{lepname}_Lep_eta_{strname}", 50, -2.5, 2.5, DrawConfig(xmin=-2.5, xmax=2.5, xlabel="Lepton #eta", ylabel=f"Events / {(2.5+2.5)/50:.2f}", dology=True, ymax=ymax, donormalizebin=False, addOverflow=True, addUnderflow=True, showratio=False, legendPos=[0.94, 0.88, 0.70, 0.68]), weightname=strname)
        
            sampMan.cacheDraw("met_pt", f"histo_wjets_{lepname}_met_pt_{strname}", 70, 0, 70, DrawConfig(xmin=0, xmax=70, xlabel="MET [GeV]", ylabel=f"Events / GeV", dology=True, ymax=ymax, donormalizebin=False, addOverflow=True, addUnderflow=True, showratio=False, legendPos=[0.94, 0.88, 0.70, 0.68]), weightname=strname)
        
            sampMan.cacheDraw("deltaPhi", f"histo_wjets_{lepname}_deltaPhi_{strname}", 30, 0, ROOT.TMath.Pi(), DrawConfig(xmin=0, xmax=ROOT.TMath.Pi(), xlabel="#Delta#phi(#mu, MET)", ylabel=f"Events / {(2*ROOT.TMath.Pi())/30:.2f}", dology=True, ymax=ymax, donormalizebin=False, addOverflow=True, addUnderflow=True, showratio=False, legendPos=[0.94, 0.88, 0.70, 0.68]), weightname=strname)
        
            sampMan.cacheDraw2D("Lep_pt", "Lep_eta", f"histo_wjets_{lepname}_pt_vs_eta_{strname}", 9, 25, 70, 8, -2.4, 2.4, DrawConfig(xmin=25, xmax=70, ymin=-2.4, ymax=2.4, xlabel="Lepton p_{T} [GeV]", ylabel="Lepton #eta", zlabel="Events / 10 GeV", dology=False, donormalizebin=False, addOverflow=True, addUnderflow=True, showratio=False, legendPos=[0.94, 0.88, 0.70, 0.68]), weightname=strname)
            sampMan.cacheDraw2D("Lep_pt", "met_pt", f"histo_wjets_{lepname}_pt_vs_met_{strname}", 9, 25, 70, 10, 0, 70, DrawConfig(xmin=25, xmax=70, ymin=0, ymax=70, xlabel="Lepton p_{T} [GeV]", ylabel="MET [GeV]", zlabel="Events / 10 GeV", dology=False, donormalizebin=False, addOverflow=True, addUnderflow=True, showratio=False, legendPos=[0.94, 0.88, 0.70, 0.68]), weightname=strname)
            sampMan.cacheDraw2D("Lep_pt", "deltaPhi", f"histo_wjets_{lepname}_pt_vs_deltaPhi_{strname}", 9, 25, 70, 10, 0., ROOT.TMath.Pi(), DrawConfig(xmin=25, xmax=70, ymin=0, ymax=ROOT.TMath.Pi(), xlabel="Lepton p_{T} [GeV]", ylabel="#Delta#phi(#mu, MET)", zlabel="Events / 10 GeV", dology=False, donormalizebin=False, addOverflow=True, addUnderflow=True, showratio=False, legendPos=[0.94, 0.88, 0.70, 0.68]), weightname=strname)

            sampMan.cacheDraw("RelIso", "histo_qcd_RelIso", 75, 0, 0.75, DrawConfig(xmin=0, xmax=0.75, xlabel="Relative Isolation", ylabel="Events / 0.01", dology=True, ymax=1e5, donormalizebin=False, addOverflow=True, addUnderflow=True, showratio=False, legendPos=[0.94, 0.88, 0.70, 0.68]))
    
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
                        
                        ymax = 1e5
                        makeDraws(strname, ymax)
                        
                        for mt in mtbins:
                            strname = "weight_{}_{}_{}_{}_{}".format(chg, iso, wpt, lepeta, mt)
                            sampMan.DefineAll(strname, "w_{} * {} * {} * {} * {}".format(iso, wpt, lepeta, chg, mt))
                            makeDraws(strname, ymax)

        sampMan.launchDraw()
        
        hIsos = OrderedDict()
        hIsos_mt = OrderedDict()
    
        hleppts = OrderedDict()
        hleppts_mt = OrderedDict()
        hlepetas = OrderedDict()
        hlepetas_mt = OrderedDict()
        hmetpts = OrderedDict()
        hmetpts_mt = OrderedDict()
        hdeltaphis = OrderedDict()
        hdeltaphis_mt = OrderedDict()
    
        hpt_vs_etas_data = OrderedDict()
        hpt_vs_etas_mc = OrderedDict()
        hpt_vs_etas_data_mt = OrderedDict()
        hpt_vs_etas_mc_mt = OrderedDict()
    
        hpt_vs_mets_data = OrderedDict()
        hpt_vs_mets_mc = OrderedDict()
        hpt_vs_mets_data_mt = OrderedDict()
        hpt_vs_mets_mc_mt = OrderedDict()
    
        hpt_vs_deltaphis_data = OrderedDict()
        hpt_vs_deltaphis_mc = OrderedDict()
        hpt_vs_deltaphis_data_mt = OrderedDict()
        hpt_vs_deltaphis_mc_mt = OrderedDict()
    
        def getHisto(name):
            hdata = sampMan.hdatas[name]
            #hstacked = THStack2TH1(sampMan.hsmcs[name])
            #for ibin in range(hstacked.GetNbinsX()+1):
            #    hstacked.SetBinContent(ibin, max(hstacked.GetBinContent(ibin), 0))
            #hdata.Add(hstacked, -1.0)
            hdata.SetName(name)
            return hdata
    
        def getHisto2D(name):
            """
            return both data and merged mc
            """
            hdata = sampMan.hdatas2D[name]
            #hmc = sampMan.hmcs2D[name]
            hdata.SetName(name + "_Data")
            #hmc.SetName(name + "_MC")
            # dumb histogram as a placeholder
            hmc = hdata.Clone(name + "_MC")
            for ibinx in range(0, hmc.GetNbinsX()+2):
                for ibiny in range(0, hmc.GetNbinsY()+2):
                    hmc.SetBinContent(ibinx, ibiny, 0.)
                    hmc.SetBinError(ibinx, ibiny, 0.)
            return hdata, hmc
    
        def dumpHistos(hlist, outputname, doPrint = False):
            outfile = ROOT.TFile.Open(f"root/{outputname}", "recreate")
            if not isinstance(hlist, list):
                hlist = [hlist]
            for hs in hlist:
                for isobin, h in hs.items():
                    if doPrint:
                        print(f"{isobin} mean: {h.GetMean():.3f}")
                    h.SetDirectory(outfile)
                    h.Write()
            outfile.Close()

        # hetas_mtCut_comp = OrderedDict()
        for iso in isobins:
            for wpt in wptbins:
                for lepeta in etabins:
                    for chg in chgbins:
                        strname = "weight_{}_{}_{}_{}".format(
                            chg, iso, wpt, lepeta)

                        # for isolation
                        outputname = f"histo_wjets_{lepname}_RelIso_{strname}"
                        hIsos[strname] = getHisto(outputname)

                        outputname = f"histo_wjets_{lepname}_Lep_pt_{strname}"
                        hleppts[strname] = getHisto(outputname)

                        outputname = f"histo_wjets_{lepname}_Lep_eta_{strname}"
                        hlepetas[strname] = getHisto(outputname)

                        outputname = f"histo_wjets_{lepname}_met_pt_{strname}"
                        hmetpts[strname] = getHisto(outputname)

                        outputname = f"histo_wjets_{lepname}_deltaPhi_{strname}"
                        hdeltaphis[strname] = getHisto(outputname)

                        outputname = f"histo_wjets_{lepname}_pt_vs_eta_{strname}"
                        hpt_vs_etas_data[strname], hpt_vs_etas_mc[strname] = getHisto2D(outputname)

                        outputname = f"histo_wjets_{lepname}_pt_vs_met_{strname}"
                        hpt_vs_mets_data[strname], hpt_vs_mets_mc[strname] = getHisto2D(outputname)

                        outputname = f"histo_wjets_{lepname}_pt_vs_deltaPhi_{strname}"
                        hpt_vs_deltaphis_data[strname], hpt_vs_deltaphis_mc[strname] = getHisto2D(outputname)

                        for mt in mtbins:
                            strname = "weight_{}_{}_{}_{}_{}".format(
                                chg, iso, wpt, lepeta, mt)

                            # for isolation
                            outputname = f"histo_wjets_{lepname}_RelIso_{strname}"
                            hIsos_mt[strname] = getHisto(outputname)

                            outputname = f"histo_wjets_{lepname}_Lep_pt_{strname}"
                            hleppts_mt[strname] = getHisto(outputname)

                            outputname = f"histo_wjets_{lepname}_Lep_eta_{strname}"
                            hlepetas_mt[strname] = getHisto(outputname)

                            outputname = f"histo_wjets_{lepname}_met_pt_{strname}"
                            hmetpts_mt[strname] = getHisto(outputname)

                            outputname = f"histo_wjets_{lepname}_deltaPhi_{strname}"
                            hdeltaphis_mt[strname] = getHisto(outputname)

                            outputname = f"histo_wjets_{lepname}_pt_vs_eta_{strname}"
                            hpt_vs_etas_data_mt[strname], hpt_vs_etas_mc_mt[strname] = getHisto2D(outputname)

                            outputname = f"histo_wjets_{lepname}_pt_vs_met_{strname}"
                            hpt_vs_mets_data_mt[strname], hpt_vs_mets_mc_mt[strname] = getHisto2D(outputname)

                            outputname = f"histo_wjets_{lepname}_pt_vs_deltaPhi_{strname}"
                            hpt_vs_deltaphis_data_mt[strname], hpt_vs_deltaphis_mc_mt[strname] = getHisto2D(outputname)

        suffix = "relIso_" + lepname + "nu_QCDMC"
        if doGenMET:
            suffix += "_GenMET"
        else:
            suffix += "_PFMET"
        suffix += ".root"

        dumpHistos([hIsos, hIsos_mt], f"output_qcdIsoMean_{suffix}", doPrint=True)
        dumpHistos([hleppts, hleppts_mt], f"output_qcdLepPtMean_{suffix}", doPrint=False)
        dumpHistos([hlepetas, hlepetas_mt], f"output_qcdLepEtaMean_{suffix}", doPrint=False)
        dumpHistos([hmetpts, hmetpts_mt], f"output_qcdMetPtMean_{suffix}", doPrint=False)
        dumpHistos([hdeltaphis, hdeltaphis_mt], f"output_qcdDeltaPhiMean_{suffix}", doPrint=False)
    
        dumpHistos([hpt_vs_etas_data, hpt_vs_etas_mc, hpt_vs_etas_data_mt, hpt_vs_etas_mc_mt], f"output_qcdLepPtVsEtaMean_{suffix}", doPrint=False)
    
        dumpHistos([hpt_vs_mets_data, hpt_vs_mets_mc, hpt_vs_mets_data_mt, hpt_vs_mets_mc_mt], f"output_qcdLepPtVsMetMean_{suffix}", doPrint=False)
    
        dumpHistos([hpt_vs_deltaphis_data, hpt_vs_deltaphis_mc, hpt_vs_deltaphis_data_mt, hpt_vs_deltaphis_mc_mt], f"output_qcdLepPtVsDeltaPhiMean_{suffix}", doPrint=False) 
        
        return
    
        outfile =  "root/output_qcdMCshape_"+ suffix
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

from ast import parse
import ROOT
import numpy as np
from collections import OrderedDict
from modules.SampleManager import DrawConfig, Sample, SampleManager
from modules.Binnings import Vptbins
from modules.theoryUncIndices import theoryUnc_13TeV_wjets, theoryUnc_5TeV_wjets
from CMSPLOTS.myFunction import DrawHistos
import argparse

ROOT.gROOT.SetBatch(True)
ROOT.ROOT.EnableImplicitMT(10)

def main():
    print("Program start...")

    parser = argparse.ArgumentParser(description="Make plots for Wlnu analysis")
    parser.add_argument("--doTest", action="store_true", dest="doTest", help="Run on a subset of samples for debugging; false runs on all dataset.")
    parser.add_argument("--doWpT", action="store_true", dest="doWpT", help="Bin in different W pt bins; false runs inclusively")
    parser.add_argument("--is5TeV", action="store_true", dest="is5TeV", help="Analyze the 5TeV data; false runs on 13TeV data")
    parser.add_argument("--doElectron", action="store_true", dest="doElectron", help="Analyze the electron channel; false runs the muon channel")
    args = parser.parse_args()

    doMuon = not args.doElectron
    doTest = args.doTest
    doWpT  = args.doWpT
    is5TeV = args.is5TeV
    doTheoryNorm = True

    print("doMuon: ", doMuon)
    print("doTest:", doTest)
    print("doWpT:", doWpT)
    print("is5TeV:", is5TeV)
    print("doTheoryNorm", doTheoryNorm)

    #ROOT.gROOT.ProcessLine('TFile* f_zpt = TFile::Open("results/zpt_weight.root")')
    #ROOT.gROOT.ProcessLine('TH1D* h_zpt_ratio  = (TH1D*)f_zpt->Get("h_zpt_ratio")')

    if not is5TeV:
        if doMuon:
            input_data    =    "inputs/wmunu/input_data.txt"
            input_wl0     =    "inputs/wmunu/input_wm0.txt"
            input_wl1     =    "inputs/wmunu/input_wm1.txt"
            input_wl2     =    "inputs/wmunu/input_wm2.txt"
            input_ttbar   =    "inputs/wmunu/input_ttbar_dilepton.txt"
            input_ttbar_1lep = "inputs/wmunu/input_ttbar_singlelepton.txt"
            input_ttbar_0lep = "inputs/wmunu/input_ttbar_hadronic.txt"
            input_ww      =    "inputs/wmunu/input_ww.txt"
            input_wz      =    "inputs/wmunu/input_wz.txt"
            input_zz      =    "inputs/wmunu/input_zz.txt"
            input_zxx     =    "inputs/wmunu/input_zxx.txt"
            input_wx0     =    "inputs/wmunu/input_wx0.txt"
            input_wx1     =    "inputs/wmunu/input_wx1.txt"
            input_wx2     =    "inputs/wmunu/input_wx2.txt"
        else:
            input_data    = "inputs/wenu/input_data.txt"
            input_wl0     = "inputs/wenu/input_we0.txt"
            input_wl1     = "inputs/wenu/input_we1.txt"
            input_wl2     = "inputs/wenu/input_we2.txt"
            input_ttbar   = "inputs/wenu/input_ttbar_dilepton.txt"
            input_ttbar_1lep = "inputs/wenu/input_ttbar_singlelepton.txt"
            input_ttbar_0lep = "inputs/wenu/input_ttbar_hadronic.txt"
            input_ww      = "inputs/wenu/input_ww.txt"
            input_wz      = "inputs/wenu/input_wz.txt"
            input_zz      = "inputs/wenu/input_zz.txt"
            input_zxx     = "inputs/wenu/input_zxx.txt"
            input_wx0     = "inputs/wenu/input_wx0.txt"
            input_wx1     = "inputs/wenu/input_wx1.txt"
            input_wx2     = "inputs/wenu/input_wx2.txt"
    else:
        if doMuon:
            input_data  = "inputs_5TeV/wmunu/input_data.txt"
            input_wl    = "inputs_5TeV/wmunu/input_wm.txt"
            input_ttbar = "inputs_5TeV/wmunu/input_ttbar.txt"
            input_ww    = "inputs_5TeV/wmunu/input_ww.txt"
            input_wz    = "inputs_5TeV/wmunu/input_wz.txt"
            input_zz2l  = "inputs_5TeV/wmunu/input_zz2l.txt"
            input_zz4l  = "inputs_5TeV/wmunu/input_zz4l.txt" 
            input_zxx   = "inputs_5TeV/wmunu/input_zxx.txt"
            input_wx    = "inputs_5TeV/wmunu/input_wx.txt"
        else:
            input_data  = "inputs_5TeV/wenu/input_data.txt"
            input_wl    = "inputs_5TeV/wenu/input_we.txt"
            input_ttbar = "inputs_5TeV/wenu/input_ttbar.txt"
            input_ww    = "inputs_5TeV/wenu/input_ww.txt"
            input_wz    = "inputs_5TeV/wenu/input_wz.txt"
            input_zz2l  = "inputs_5TeV/wenu/input_zz2l.txt"
            input_zz4l  = "inputs_5TeV/wenu/input_zz4l.txt"
            input_zxx   = "inputs_5TeV/wenu/input_zxx.txt"
            input_wx    = "inputs_5TeV/wenu/input_wx.txt"

    DataSamp  = Sample(input_data, isMC=False, legend="Data", name="Data", isWSR=True, doTheoryVariation=False)
    signame = "wlnu"
    if not is5TeV:
        # W -> munu
        #wjetsnorm = 1.06
        wjetsnorm = 1.0
        Wl0Samp   = Sample(input_wl0, isMC=True, name = "wl0", isWSR=True, additionalnorm = wjetsnorm)
        Wl1Samp   = Sample(input_wl1, isMC=True, name = "wl1", isWSR=True, additionalnorm = wjetsnorm)
        Wl2Samp   = Sample(input_wl2, isMC=True, name = "wl2", isWSR=True, additionalnorm = wjetsnorm)
        # ttbar
        TTbarSamp  = Sample(input_ttbar, isMC=True, name = "ttbar_dilepton", isWSR=True)
        TT1LepSamp = Sample(input_ttbar_1lep, isMC=True, name = "ttbar_1lepton", isWSR=True)
        TT0LepSamp = Sample(input_ttbar_0lep, isMC=True, name = "ttbar_0lepton", isWSR=True)
        ## dibosons
        WWSamp = Sample(input_ww, isMC=True, name = "WW", isWSR=True)
        WZSamp = Sample(input_wz, isMC=True, name = "WZ", isWSR=True)
        ZZSamp = Sample(input_zz, isMC=True, name = "ZZ", isWSR=True)
        # tau
        ZXXSamp = Sample(input_zxx, isMC=True, name = "ZXX", isWSR=True)
        Wx0Samp = Sample(input_wx0, isMC=True, name = "wx0", isWSR=True)
        Wx1Samp = Sample(input_wx1, isMC=True, name = "wx1", isWSR=True)
        Wx2Samp = Sample(input_wx2, isMC=True, name = "wx2", isWSR=True)

        if not doTest:
            sampMan = SampleManager(DataSamp, [Wl0Samp, Wl1Samp, Wl2Samp, TTbarSamp, TT1LepSamp, TT0LepSamp, WWSamp, WZSamp, ZZSamp, ZXXSamp, Wx0Samp, Wx1Samp, Wx2Samp])
        else:
            sampMan = SampleManager(DataSamp, [Wl1Samp, Wl2Samp])
        wxname = "wx"
        sampMan.groupMCs(["wx0", "wx1", "wx2"], wxname, 216, 'wx')
        sampMan.groupMCs(["WW", "WZ", "ZZ"], 'vv', 216, 'vv')
        sampMan.groupMCs(['ZXX'], 'zxx', 216, 'zxx')
        sampMan.groupMCs(["ttbar_dilepton", "ttbar_1lepton", "ttbar_0lepton"], "ttbar", 96, "t#bar{t}")
        label = "W#rightarrow#mu#nu" if doMuon else "W#rightarrow e#nu"
        sampMan.groupMCs(['wl0', 'wl1', 'wl2'], signame, 92,"W#rightarrow#mu#nu")

        # for the signal samples
        if not doTest:
            signalSamps = [Wl0Samp, Wl1Samp, Wl2Samp]
        else:
            signalSamps = [Wl1Samp, Wl2Samp]
        signalSampnames = [samp.name for samp in signalSamps]

    else:
        # for the 5TeV files
        # W -> munu
        wjetsnorm = 1.0
        WlSamp     = Sample(input_wl, isMC=True, name = signame, color = 92, legend = "W#rightarrow#mu#nu", isWSR=True, additionalnorm = wjetsnorm, is5TeV = True)
        # ttbar
        TTbarSamp  = Sample(input_ttbar, isMC=True, name = "ttbar", color = 86, legend = "t#bar{t}", isWSR=True, is5TeV = True)
        ## dibosons
        WWSamp = Sample(input_ww, isMC=True, name = "WW", isWSR=True, is5TeV = True)
        WZSamp = Sample(input_wz, isMC=True, name = "WZ", isWSR=True, is5TeV = True)
        ZZ2LSamp = Sample(input_zz2l, isMC=True, name = "ZZ2L", isWSR=True, is5TeV = True)
        ZZ4LSamp = Sample(input_zz4l, isMC=True, name = "ZZ4L", isWSR=True, is5TeV = True)
        # tau
        ZXXSamp = Sample(input_zxx, isMC=True, name = "ZXX", isWSR=True, is5TeV = True)
        WxSamp = Sample(input_wx, isMC=True, name = "wx", color = 216, legend = "wx", isWSR=True, is5TeV = True)

        sampMan = SampleManager(DataSamp, [WlSamp, TTbarSamp, WWSamp, WZSamp, ZZ2LSamp, ZZ4LSamp, ZXXSamp, WxSamp], is5TeV = True)
        sampMan.groupMCs(["WW", "WZ", "ZZ2L", "ZZ4L"], 'vv', 75, 'vv')
        sampMan.groupMCs(['ZXX'], 'zxx', 84, 'zxx')
        label = "W#rightarrow#mu#nu" if doMuon else "W#rightarrow e#nu"

        # for the signal samples
        signalSamps = [WlSamp]
        signalSampnames = [samp.name for samp in signalSamps]
    print(("signal sample names: ", signalSampnames))
    h_sigs = OrderedDict()

    # define variables and weights
    sampMan.DefineAll("Lep_pt",  "lep.Pt()")
    sampMan.DefineAll("Lep_eta", "lep.Eta()")
    sampMan.DefineAll("Lep_phi", "lep.Phi()")
    for i in [1, 2, 3, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]:
        # variations on the recoil correction
        sampMan.DefineMC("met_{}_pt".format(str(i)),  "metVars[{}]".format(str(i))) 
        sampMan.DefineMC("met_{}_phi".format(str(i)), "metVarsPhi[{}]".format(str(i)))
        # for data, always use the uncorrected
        DataSamp.Define("met_{}_pt".format(str(i)),    "metVars[0]")
        DataSamp.Define("met_{}_phi".format(str(i)),   "metVarsPhi[0]")

        # mT
        sampMan.DefineAll("mT_{}".format(str(i)), "TMath::Sqrt(2 * Lep_pt * met_{}_pt * (1- TMath::Cos(Lep_phi - met_{}_phi)))".format(str(i), str(i)))

    # deltaPhi
    sampMan.DefineAll("deltaPhi", "TVector2::Phi_0_2pi(lep.Phi() - met_1_phi)")

    # recoil variables
    sampMan.DefineAll("V2W", "UVec(lep.Pt(), lep.Phi(), met_1_pt, met_1_phi)")
    sampMan.DefineAll("WpT", "V2W.Mod()")
    sampMan.DefineAll("Wphi", "TVector2::Phi_mpi_pi(V2W.Phi())")

    # charge
    lepname = "mu" if doMuon else "e"
    sampMan.DefineAll(lepname+"plus",  "q > 0")
    sampMan.DefineAll(lepname+"minus", "q < 0")
    chgbins = [lepname+"plus", lepname + "minus"]

    # WpT bins
    # bin0 is the inclusive bin
    sampMan.DefineAll("WpT_bin0",  "WpT>=0.")
    wptbins = ["WpT_bin0"]
    if doWpT:
        sampMan.DefineAll("WpT_bin1",  "WpT>=0. && WpT<8.0")
        sampMan.DefineAll("WpT_bin2",  "WpT>=8.0 &&  WpT<16.0")
        sampMan.DefineAll("WpT_bin3",  "WpT>=16.0 && WpT<24.0")
        sampMan.DefineAll("WpT_bin4",  "WpT>=24.0 && WpT<32.0")
        sampMan.DefineAll("WpT_bin5",  "WpT>=32.0 && WpT<40.0")
        sampMan.DefineAll("WpT_bin6",  "WpT>=40.0 && WpT<50.0")
        sampMan.DefineAll("WpT_bin7",  "WpT>=50.0 && WpT<70.0")
        sampMan.DefineAll("WpT_bin8",  "WpT>=70.0 && WpT<100.0")
        sampMan.DefineAll("WpT_bin9",  "WpT>=100.0")
        wptbins = ["WpT_bin0", "WpT_bin1", "WpT_bin2", "WpT_bin3", "WpT_bin4", "WpT_bin5", "WpT_bin6", "WpT_bin7", "WpT_bin8", "WpT_bin9"]
        if doTest:
            wptbins = ["WpT_bin0", "WpT_bin1"]

    # truth level
    wpttruthbins = []
    if doWpT:
        sampMan.DefineSpecificMCs("WpT_truth", "genV.Pt()", sampnames=signalSampnames)
        sampMan.DefineSpecificMCs("WpT_truth_bin1",  "WpT_truth>=0.   &&  WpT_truth<8.0",  sampnames=signalSampnames)
        sampMan.DefineSpecificMCs("WpT_truth_bin2",  "WpT_truth>=8.0  &&  WpT_truth<16.0", sampnames=signalSampnames)
        sampMan.DefineSpecificMCs("WpT_truth_bin3",  "WpT_truth>=16.0 &&  WpT_truth<24.0", sampnames=signalSampnames)
        sampMan.DefineSpecificMCs("WpT_truth_bin4",  "WpT_truth>=24.0 &&  WpT_truth<32.0", sampnames=signalSampnames)
        sampMan.DefineSpecificMCs("WpT_truth_bin5",  "WpT_truth>=32.0 &&  WpT_truth<40.0", sampnames=signalSampnames)
        sampMan.DefineSpecificMCs("WpT_truth_bin6",  "WpT_truth>=40.0 &&  WpT_truth<50.0", sampnames=signalSampnames)
        sampMan.DefineSpecificMCs("WpT_truth_bin7",  "WpT_truth>=50.0 &&  WpT_truth<70.0", sampnames=signalSampnames)
        sampMan.DefineSpecificMCs("WpT_truth_bin8",  "WpT_truth>=70.0 &&  WpT_truth<100.0",sampnames=signalSampnames)
        sampMan.DefineSpecificMCs("WpT_truth_bin9",  "WpT_truth>=100.0",                   sampnames=signalSampnames)
        wpttruthbins  = ["WpT_truth_bin1", "WpT_truth_bin2", "WpT_truth_bin3", "WpT_truth_bin4", "WpT_truth_bin5", "WpT_truth_bin6", "WpT_truth_bin7", "WpT_truth_bin8", "WpT_truth_bin9"]
        if doTest:
            wpttruthbins = ["WpT_truth_bin1", "WpT_truth_bin2"]


    # eta bins for electrons: barral and endcap
    sampMan.DefineAll("lepEta_bin0", "1.0")
    sampMan.DefineAll("lepEta_bin1", "abs(lep.Eta()) <= 1.4442") 
    sampMan.DefineAll("lepEta_bin2", "abs(lep.Eta()) > 1.4442")

    #if doMuon:
    #    etabins = ["lepEta_bin0"]
    #else:
    #    etabins = ["lepEta_bin0", "lepEta_bin1", "lepEta_bin2"]
    etabins = ["lepEta_bin0"]
    
    # sample weight
    # 0 is the central one with all corrections
    for i in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]:
        for wpt in wptbins:
            for lepeta in etabins:
                #if i!= 5:
                #    sampMan.DefineMC("weight_{}_{}_{}".format(str(i), wpt, lepeta), "evtWeight[{}] * self.fnorm * {} * {}".format(str(i), wpt, lepeta))
                #else:
                #    sampMan.DefineMC("weight_{}_{}_{}".format(str(i), wpt, lepeta), "(evtWeight[0] + TMath::Sqrt(TMath::Abs(evtWeight[{}]))) * self.fnorm * {} * {}".format(str(i), wpt, lepeta))
                sampMan.DefineMC("weight_{}_{}_{}".format(str(i), wpt, lepeta), "(TMath::IsNaN(evtWeight[{}]) ? 0.: evtWeight[{}]) * self.fnorm * {} * {}".format(str(i), str(i), wpt, lepeta))
                DataSamp.Define("weight_{}_{}_{}".format(str(i), wpt, lepeta),  "1.0 * {} * {}".format(wpt, lepeta))

                for chg in chgbins:
                    sampMan.DefineAll("weight_{}_{}_{}_{}".format(chg, str(i),  wpt, lepeta), "weight_{}_{}_{} * {}".format(str(i), wpt, lepeta, chg))

                    # for signal MC, define the truth WpT bins
                    for wpttruth in wpttruthbins:
                        for samp in signalSamps:
                            samp.Define("weight_{}_{}_{}_{}_{}".format(chg, str(i), wpt, lepeta, wpttruth), "weight_{}_{}_{}_{} * {}".format(chg, str(i),  wpt, lepeta, wpttruth))

    # scale the theory variations back
    # some theory variations are far from 1.0
    # use this to scale them back
    # this should be equavalent as using fnorm_varied
    for i in range(111):
        sampMan.DefineMC(f"lheweightS_{i}", f"lheweight[{i}] * self.nmcevt / self.nmcevts_varied[{i}]")

    vsamples = ["wl0", "wl1", "wl2", "wx0", "wx1", "wx2"]
    sampMan.DefineSpecificMCs("VpT", "genV.Pt()", sampnames=vsamples)

    theoryMaps = theoryUnc_13TeV_wjets if not is5TeV else theoryUnc_5TeV_wjets
    theoryVariations = OrderedDict()

    for wpt in range(len(Vptbins)-1):
        ptmin = Vptbins[wpt]
        ptmax = Vptbins[wpt+1]
        for var in ["MuRVPTUp", "MuRVPTDown", "MuFVPTUp", "MuFVPTDown", "MuFMuRVPTUp", "MuFMuRVPTDown"]:
            idx = theoryMaps.getIndex(var.replace("VPT", ""))
            sampMan.DefineSpecificMCs(f"{var}".replace("VPT", str(wpt)), f"(VpT >= {ptmin} && VpT < {ptmax}) ? lheweightS_{idx}: 1.0", sampnames=vsamples)
            sampMan.DefineAll(f"{var}".replace("VPT", str(wpt)), "1.0", excludes=vsamples)

            theoryVariations[f"{var}".replace("VPT", str(wpt))] = f"{var}".replace("VPT", str(wpt))

    # pdf
    pdfStart = theoryMaps.getIndex("PDFStart")
    for i in range(1, 101):
        theoryVariations[f"PDF{i}Up"] = f"lheweightS_{i + pdfStart - 1}"
    # alphaS
    alphaSUp = theoryMaps.getIndex("alphaSUp")
    theoryVariations["alphaSUp"] = f"lheweightS_{alphaSUp}"

    # weights with theory variations
    for wpt in wptbins:
        for lepeta in etabins:
            for chg in chgbins:
                for var in theoryVariations.values():
                    sampMan.DefineMC(f"weight_{chg}_{var}_{wpt}_{lepeta}", f"(TMath::IsNaN(evtWeight[0])) ? 0. : evtWeight[0] * self.fnorm * {wpt} * {lepeta} * {var} * {chg}")
                    DataSamp.Define( f"weight_{chg}_{var}_{wpt}_{lepeta}", f"1.0 * {wpt} * {lepeta} * {chg}")

    #
    # Some histograms for simple data-MC comparison
    #
    met_pt_bins = np.array([0., 2.0, 4., 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 33, 36, 39, 42, 45, 48, 51, 55, 60, 65, 70, 75, 80])
    mt_bins = np.concatenate((np.linspace(40, 90, 51), np.linspace(92, 100, 5), np.linspace(104, 120, 5)))
    phimin = -ROOT.TMath.Pi()
    phimax = ROOT.TMath.Pi()
    dphi_bins = np.concatenate((np.linspace(0,1, 3), np.linspace(1.2, 2.4, 7),np.linspace(2.5, 3.5, 11), np.linspace(3.6, 4.8, 7), np.linspace(5.0, 6.0, 3))) * ROOT.TMath.Pi() / 3

    eta_bins = np.concatenate((np.linspace(-2.4, -2.0, 3),np.linspace(-1.8,1.8,37), np.linspace(2.0, 2.4, 3)))
    pt_bins = np.concatenate((np.linspace(25, 35, 6),np.linspace(36,55,20), np.linspace(57, 63, 4), np.linspace(66,70,2)))
    #mass_bins = np.concatenate((np.linspace(70, 76, 3),np.linspace(78,86,5), np.linspace(87, 96, 10), np.linspace(98,104,5), np.linspace(106, 110, 2)))
    # do the fine binning first; then rebin in the processHists
    mass_bins = np.array([0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0, 95.0, 100.0, 105.0, 110.0, 115.0, 120.0, 125.0, 130.0, 135.0, 140.])

    # nPV
    sampMan.cacheDraw("npv", "histo_nPV_"+lepname, 10, 0, 10, DrawConfig(xmin=0, xmax=10, xlabel='# PV', ylabel='Events / 1', dology=False, ymax = 9e5), weightname = "weight_0_WpT_bin0_lepEta_bin0")
    # lepton pt
    sampMan.cacheDraw("Lep_pt", "histo_Lep_pt_"+lepname, pt_bins, DrawConfig(xmin=20, xmax=70, xlabel='p_{T}(#mu) [GeV]', dology=False, ymax=1.05e5), weightname="weight_0_WpT_bin0_lepEta_bin0")
    sampMan.cacheDraw("Lep_eta", "histo_Lep_eta_"+lepname, eta_bins, DrawConfig(xmin=-2.6, xmax=2.6, xlabel='#eta (#mu) [GeV]', ylabel='Events / 1', dology=False, ymax=6e5), weightname = "weight_0_WpT_bin0_lepEta_bin0")
    sampMan.cacheDraw("Lep_pt", "histo_Lep_pt_"+lepname+"plus", pt_bins, DrawConfig(xmin=20, xmax=70, xlabel='p_{T}(#mu) [GeV]', dology=False, ymax=6.05e4), weightname = "weight_"+lepname+"plus_0_WpT_bin0_lepEta_bin0")
    sampMan.cacheDraw("Lep_pt", "histo_Lep_pt_"+lepname+"minus", pt_bins, DrawConfig(xmin=20, xmax=70, xlabel='p_{T}(#mu) [GeV]', dology=False, ymax=5.05e4), weightname = "weight_"+lepname+"minus_0_WpT_bin0_lepEta_bin0")

    # met and raw_met
    sampMan.cacheDraw("met_1_pt", "histo_wjets_"+lepname+"_pfmet_pt", met_pt_bins, DrawConfig(xmin=0, xmax=80, xlabel='PF MET [GeV]', dology=False, ymax=8e4), weightname="weight_0_WpT_bin0_lepEta_bin0")
    sampMan.cacheDraw("met_1_phi", "histo_wjets_"+lepname+"_pfmet_phi", 30, phimin, phimax, DrawConfig(xmin=phimin, xmax=phimax, xlabel='PF MET #phi', dology=False, ymax=6e5), weightname = "weight_0_WpT_bin0_lepEta_bin0")

    # deltaPhi
    sampMan.cacheDraw("deltaPhi", "histo_wjets_"+lepname+"_deltaphi",      dphi_bins,  DrawConfig(xmin=0., xmax=2*phimax, xlabel='#Delta #phi(#mu, MET)', dology=False, ymax=2e6), weightname="weight_0_WpT_bin0_lepEta_bin0")

    # mt
    sampMan.cacheDraw("mtCorr", "histo_wjets_"+lepname+"_mtcorr",      mt_bins, DrawConfig(xmin=0, xmax=120, xlabel="m_{T} [GeV]", dology=False, ymax=7e4), weightname="weight_0_WpT_bin0_lepEta_bin0")
    sampMan.cacheDraw("mtCorr", "histo_wjets_"+lepname+"plus_mtcorr",  mt_bins, DrawConfig(xmin=0, xmax=120, xlabel="m_{T} [GeV]", dology=False, ymax=4e4), weightname = "weight_"+lepname+"plus_0_WpT_bin0_lepEta_bin0")
    sampMan.cacheDraw("mtCorr", "histo_wjets_"+lepname+"minus_mtcorr", mt_bins, DrawConfig(xmin=0, xmax=120, xlabel="m_{T} [GeV]", dology=False, ymax=3.5e4), weightname = "weight_"+lepname+"minus_0_WpT_bin0_lepEta_bin0")

    # wpt
    sampMan.cacheDraw("WpT", "histo_wjets_"+lepname+"_wpt",        300, 0,  300, DrawConfig(xmin=0, xmax=300, xlabel="p^{W}_{T} [GeV]", dology=True, ymax=1e6, ymin=1e1), weightname="weight_0_WpT_bin0_lepEta_bin0")
    sampMan.cacheDraw("WpT", "histo_wjets_"+lepname+"plus_wpt",    300, 0,  300, DrawConfig(xmin=0, xmax=300, xlabel="p^{W}_{T} [GeV]", dology=True, ymax=1e6, ymin=1e1), weightname="weight_"+lepname+"plus_0_WpT_bin0_lepEta_bin0")
    sampMan.cacheDraw("WpT", "histo_wjets_"+lepname+"minus_wpt",   300, 0,  300, DrawConfig(xmin=0, xmax=300, xlabel="p^{W}_{T} [GeV]", dology=True, ymax=1e6, ymin=1e1), weightname="weight_"+lepname+"minus_0_WpT_bin0_lepEta_bin0")


    ymaxs = {lepname+"plus": 4e4, lepname+"minus": 3.5e4}
    nbins = 12
    xmin = 0
    xmax = 120

    #
    # templates and variations for combine
    #
    for chg in chgbins:
        # variation on recoil corections
        for i in [1, 2, 3, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]:
            for wpt in wptbins:
                for lepeta in etabins:
                    sampMan.cacheDraw("mT_{}".format(str(i)), "histo_wjets_{}_mT_{}_{}_{}".format(chg, str(i), wpt, lepeta),  mass_bins, DrawConfig(xmin=xmin, xmax=xmax, xlabel="m_{T} [GeV]", dology=False, ymax=ymaxs[chg], donormalizebin=False, addOverflow=False, addUnderflow=False), weightname = "weight_{}_0_{}_{}".format(chg, wpt, lepeta))

                    # for signal MC, draw the histogram in different truth WpT bins
                    for wpttruth in wpttruthbins:
                        h_list = []
                        hname = "histo_wjets_{}_mT_{}_{}_{}_{}_signalMC".format(chg, str(i), wpt, lepeta, wpttruth)
                        for isamp, samp in enumerate(signalSamps):    
                            h_list.append( samp.rdf.Histo1D((hname+str(isamp), hname, mass_bins), "mT_{}".format(str(i)), "weight_{}_0_{}_{}_{}".format(chg, wpt, lepeta, wpttruth)) )
                        h_sigs[hname] = h_list
                            

        # variation on the scale factors and sample weights
        for i in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]:
            for wpt in wptbins:
                for lepeta in etabins:
                    sampMan.cacheDraw("mT_1", f"histo_wjets_{chg}_mtcorr_weight_{i}_{wpt}_{lepeta}", mass_bins, DrawConfig(xmin=xmin, xmax=xmax, xlabel="m_{T} [GeV]", dology=False, ymax=ymaxs[chg], donormalizebin=False, addOverflow=False, addUnderflow=False), weightname = f"weight_{chg}_{i}_{wpt}_{lepeta}")

                    # for signal MC, draw the histogram in different truth WpT bins
                    for wpttruth in wpttruthbins:
                        h_list = []
                        hname = "histo_wjets_{}_mtcorr_weight_{}_{}_{}_{}_signalMC".format(chg, str(i), wpt, lepeta, wpttruth)
                        for isamp, samp in enumerate(signalSamps):
                            h_list.append( samp.rdf.Histo1D((hname+str(isamp), hname, mass_bins), "mT_1", "weight_{}_{}_{}_{}_{}".format(chg, str(i), wpt, lepeta, wpttruth)))
                        h_sigs[hname] = h_list
        
        # variations on the theory uncertainties
        for wpt in wptbins:
            for lepeta in etabins:
                for var in theoryVariations.values():
                    sampMan.cacheDraw("mT_1", f"histo_wjets_{chg}_mtcorr_weight_{var}_{wpt}_{lepeta}", mass_bins, DrawConfig(xmin=xmin, xmax=xmax, xlabel="m_{T} [GeV]", dology=False, ymax=ymaxs[chg], donormalizebin=False, addOverflow=False, addUnderflow=False), weightname = f"weight_{chg}_{var}_{wpt}_{lepeta}")

    # Draw all these histograms
    sampMan.launchDraw()

    #sampMan.cacheDraw("WpT", "histo_wjets_"+lepname+"minus_wpt_test",   300, 0,  300, DrawConfig(xmin=0, xmax=300, xlabel="p^{W}_{T} [GeV]", dology=True, ymax=1e6, ymin=1e1), weightname="weight_"+lepname+"minus_0_WpT_bin0_lepEta_bin0")
    #sampMan.launchDraw()

    #
    # draw the mT distributions of different processes for comparisons
    #
    for chg in chgbins:
        for wpt in wptbins:
            for lepeta in etabins:
                labels = []
                hists = []
                hmcs = sampMan.hmcs["histo_wjets_{}_mT_1_{}_{}".format(chg, wpt, lepeta)]
                icolor = 1
                for key, val in hmcs.items():
                    if any(x in key for x in ["WW", "WZ", "ZZ", "ttbar_0lepton", "ttbar_dilepton"]):
                        continue
                    labels.append(key)
                    hists.append(val)
                    val.SetLineColor(icolor)
                    val.SetMarkerColor(icolor)
                    icolor += 1
                DrawHistos(hists, labels, xmin, xmax, "m_{T} [GeV]", 0., 0.24, "A.U.", "Comparison_mTShape_wjets_{}_mT_1_{}_{}".format(chg, wpt, lepeta), dology=False, addOverflow=True, addUnderflow=True, donormalize = True, is5TeV = is5TeV, legendPos=[0.30, 0.92, 0.90, 0.72], legendNCols=4, legendoptions="LE", drawoptions="LE") 

    # 
    # merge the signal MC histograms in different samples
    # mostly for WpT usage
    #
    def addSignals(hname):
        h_list = h_sigs[hname]
        hadded = h_list[0].Clone(hname)
        # merge the signal MCs
        for idx in range(1, len(h_list)):
            hadded.Add(h_list[idx].GetValue())
        return hadded

    h_sigsMerged = OrderedDict()
    for hname in list(h_sigs.keys()):
        hadded = addSignals(hname)
        h_sigsMerged[hname] = hadded

    #
    # preparing the outputs
    #
    sqrtS = "5TeV" if is5TeV else "13TeV"
    output_suffix = lepname + "nu"
    if doWpT:
        output_suffix += "_WpT"
    output_suffix += f"_{sqrtS}.root"

    #
    # write out mT histograms for combine
    #
    outfile = ROOT.TFile.Open("root/output_shapes_"+output_suffix, "recreate")

    # Data
    #odir = outfile.mkdir("Data")
    #outfile.cd("Data")
    for chg in chgbins:
        for wpt in wptbins:
            for lepeta in etabins:
                hdata = sampMan.hdatas["histo_wjets_{}_mT_1_{}_{}".format(chg, wpt, lepeta)]
                #hdata.SetDirectory(odir)
                hdata.SetDirectory(outfile)
                hdata.Write()

    # MC
    for wpttruth in ["MCTemplates"] + wpttruthbins:
        # the 1st wpttruth bin is 'fake', for all MC templates
        # the other wpttruth bins are for signal MC only
        #odir = outfile.mkdir(wpttruth)
        #outfile.cd(wpttruth)
        for chg in chgbins:
            for wpt in wptbins:
                for lepeta in etabins:

                    if wpttruth == "MCTemplates":
                        hsmcs_central = sampMan.hsmcs["histo_wjets_{}_mT_1_{}_{}".format(chg, wpt, lepeta)]
                        hlists_central = list(hsmcs_central.GetHists())
                    else:
                        # signal MC
                        hname = "histo_wjets_{}_mT_1_{}_{}_{}_signalMC".format(chg, wpt, lepeta, wpttruth)
                        hlists_central = [h_sigsMerged[hname]]

                    # central values for MC
                    for ih in range(len(hlists_central)):
                        hcen = hlists_central[ih]
                        #hcen.SetDirectory(odir)
                        hcen.SetDirectory(outfile)
                        hcen.Write()

                    # recoil systematics
                    for i in [2, 3, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]:
                        if wpttruth == "MCTemplates":
                            hsmcs_up = sampMan.hsmcs["histo_wjets_{}_mT_{}_{}_{}".format(chg, str(i), wpt, lepeta)]
                            hlists_up = list(hsmcs_up.GetHists())
                        else:
                            # signal MC
                            hname = "histo_wjets_{}_mT_{}_{}_{}_{}_signalMC".format(chg, str(i), wpt, lepeta, wpttruth)
                            hlists_up = [h_sigsMerged[hname]]

                        for ih in range(len(hlists_up)):
                            # loop over different simulated processes
                            hcen = hlists_central[ih]
                            hup  = hlists_up[ih]
                            #if not doMuon and i !=2 and i!=3:
                            #if i!=2 and i!=3:
                            #    suffix = lepname + "_" + lepeta + "_" + wpt + "_SysRecoil" + str(i)
                            #else:
                            #    suffix = lepname + "_" + "SysRecoil" + str(i)
                            suffix = lepname + "_SysRecoil" + str(i) + "_" + sqrtS

                            hup.SetName("{}_{}Up".format(hcen.GetName(), suffix))

                            hdn  = hcen.Clone("{}_{}Down".format(hcen.GetName(), suffix))
                            for ibin in range(1, hup.GetNbinsX()+1):
                                hdn.SetBinContent(ibin, 2*hcen.GetBinContent(ibin) - hup.GetBinContent(ibin))

                            #hup.SetDirectory(odir)
                            hup.SetDirectory(outfile)
                            hup.Write()
                            #hdn.SetDirectory(odir)
                            hdn.SetDirectory(outfile)
                            hdn.Write()

                    # weights/corrections
                    for i in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]:
                        if wpttruth == "MCTemplates":
                            hsmcs_up = sampMan.hsmcs["histo_wjets_{}_mtcorr_weight_{}_{}_{}".format(chg, str(i), wpt, lepeta)]
                            hlists_up = list(hsmcs_up.GetHists())
                        else:
                            # signal MC
                            hname = "histo_wjets_{}_mtcorr_weight_{}_{}_{}_{}_signalMC".format(chg, str(i), wpt, lepeta, wpttruth)
                            hlists_up = [h_sigsMerged[hname]]
                        for ih in range(len(hlists_up)):
                            # loop over different processes
                            hcen = hlists_central[ih]
                            hup  = hlists_up[ih]
                            prefix = ""
                            if i == 2 or i==3 or i==4 or i==12 or i==13:
                                # the systematics that should be decorrelated between lepton flavors
                                prefix = lepname + "_"
                            if i==5:
                                # stat unc, decorrelate the plus and minus
                                prefix = chg + "_"
                            if i<=6 or i==8 or i==10 or i==12 or i==13:
                                suffix = prefix + "SysWeight" + str(i) + "_" + sqrtS
                                #hup.SetName("{}_SysWeight{}Up".format(hcen.GetName(), str(i)))
                                hup.SetName("{}_{}Up".format(hcen.GetName(), suffix))
                            else:
                                suffix = prefix + "SysWeight" + str(i-1) + "_" + sqrtS
                                #hup.SetName("{}_SysWeight{}Down".format(hcen.GetName(), str(i-1)))
                                hup.SetName("{}_{}Down".format(hcen.GetName(), suffix))

                            if i<6 or i==12 or i==13:
                                # for prefire, the up and down weights are calculated seperately
                                suffix = prefix + "SysWeight" + str(i) + "_" + sqrtS
                                #hdn  = hcen.Clone("{}_SysWeight{}Down".format(hcen.GetName(), str(i)))
                                hdn  = hcen.Clone("{}_{}Down".format(hcen.GetName(), suffix))
                                for ibin in range(1, hup.GetNbinsX()+1):
                                    hdn.SetBinContent(ibin, 2*hcen.GetBinContent(ibin) - hup.GetBinContent(ibin))
                            
                            #hup.SetDirectory(odir)
                            hup.SetDirectory(outfile)
                            hup.Write()
                            if i<6 or i==12 or i==13:
                                #hdn.SetDirectory(odir)
                                hdn.SetDirectory(outfile)
                                hdn.Write()

                    # pdf, alphaS, and QCD scale uncertainties
                    for sysname, var in theoryVariations.items():
                        if wpttruth == "MCTemplates":
                            hsmcs_up = sampMan.hsmcs[f"histo_wjets_{chg}_mtcorr_weight_{var}_{wpt}_{lepeta}"]
                            hlists_up = list(hsmcs_up.GetHists())
                        else:
                            continue
                        for ih in range(len(hlists_up)):
                            # loop over different processes
                            hcen = hlists_central[ih]
                            hup  = hlists_up[ih]
                            
                            suffix = sqrtS + "_" + sysname
                            hup.SetName("{}_{}".format(hcen.GetName(), suffix))

                            if doTheoryNorm:
                                # normalize the number of events to the central values
                                # if the theoretical variations are normalized
                                # ONLY for the signals
                                if signame in hcen.GetName():
                                    hup.Scale( hcen.Integral(0, 1000000) / hup.Integral(0, 1000000) )

                            if "PDF" or "alphaS" in sysname:
                                # prepare the down variation
                                # alphaS down variation is not saved in W + jets samples
                                # use the symmetrized up variation instead
                                suffix = sqrtS + "_" + sysname.replace("Up", "Down")
                                hdn  = hcen.Clone("{}_{}".format(hcen.GetName(), suffix))
                                for ibin in range(1, hup.GetNbinsX()+1):
                                    hdn.SetBinContent(ibin, 2*hcen.GetBinContent(ibin) - hup.GetBinContent(ibin))
                                hdn.SetDirectory(outfile)
                                hdn.Write()
                            
                            hup.SetDirectory(outfile)
                            hup.Write()

    outfile.Close()

    #
    # save all the unstacked mc histograms for understanding and debugging
    #
    outfile = ROOT.TFile.Open("root/output_sigshapes_"+output_suffix, "recreate")
    for chg in chgbins:
        for wpt in wptbins:
            for lepeta in etabins:
                hmcs = sampMan.hmcs["histo_wjets_{}_mT_1_{}_{}".format(chg, wpt, lepeta)]
                for key, val in hmcs.items():
                    val.SetDirectory(outfile)
                    val.Write()
    outfile.Close()

    sampMan.dumpCounts()

    print("Program end...")

    #input()
    
    return 

if __name__ == "__main__":
   main()

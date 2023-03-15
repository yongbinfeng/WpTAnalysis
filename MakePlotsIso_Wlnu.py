"""
plot the isolation distributions in different mt bins for different processes
"""
import ROOT
import numpy as np
from collections import OrderedDict
import sys
import argparse
from CMSPLOTS.myFunction import THStack2TH1

from modules.SampleManager import DrawConfig, Sample, SampleManager
from modules.Binnings import mass_bins_test

ROOT.gROOT.SetBatch(True)
ROOT.ROOT.EnableImplicitMT(16)


def main():
    print("Program start...")

    parser = argparse.ArgumentParser(
        description="Make plots for Wlnu analysis in the lepton anti-isolated region")
    parser.add_argument("--doTest", action="store_true", dest="doTest",
                        help="Run on a subset of DataSamples for debugging; false runs on all dataset.")
    parser.add_argument("--doWpT", action="store_true", dest="doWpT",
                        help="Bin in different W pt bins; false runs inclusively")
    parser.add_argument("--is5TeV", action="store_true", dest="is5TeV",
                        help="Analyze the 5TeV data; false runs on 13TeV data")
    parser.add_argument("--doElectron", action="store_true", dest="doElectron",
                        help="Analyze the electron channel; false runs the muon channel")
    args = parser.parse_args()

    doMuon = not args.doElectron
    doTest = args.doTest
    doWpT = args.doWpT
    is5TeV = args.is5TeV
    
    reweightZpt = True

    print("doMuon: ", doMuon)
    print("doTest:", doTest)
    print("doWpT:", doWpT)
    print("is5TeV:", is5TeV)
  
    useChgIso = True
    if useChgIso: 
        isovar = "pfChIso / lep.Pt()"
        suffix = "pfChIso"
    else:
        isovar = "relIso"
        suffix = "relIso"
    
    if not is5TeV:
        if doMuon:
            ROOT.gROOT.ProcessLine(
                'TFile* f_zpt = TFile::Open("data/output_Zpt_mumu_13TeV.root")')
            ROOT.gROOT.ProcessLine(
                'TH1D* h_zpt_ratio = (TH1D*)f_zpt->Get("histo_zjets_zpt_mumu_Ratio")')
        else:
            ROOT.gROOT.ProcessLine(
                'TFile* f_zpt = TFile::Open("data/output_Zpt_ee_13TeV.root")')
            ROOT.gROOT.ProcessLine(
                'TH1D* h_zpt_ratio = (TH1D*)f_zpt->Get("histo_zjets_zpt_ee_Ratio")')
    else:
        if doMuon:
            ROOT.gROOT.ProcessLine(
                'TFile* f_zpt = TFile::Open("data/output_Zpt_mumu_5TeV.root")')
            ROOT.gROOT.ProcessLine(
                'TH1D* h_zpt_ratio = (TH1D*)f_zpt->Get("histo_zjets_zpt_mumu_Ratio")')
        else:
            ROOT.gROOT.ProcessLine(
                'TFile* f_zpt = TFile::Open("data/output_Zpt_ee_5TeV.root")')
            ROOT.gROOT.ProcessLine(
                'TH1D* h_zpt_ratio = (TH1D*)f_zpt->Get("histo_zjets_zpt_ee_Ratio")')

    outdir = "plots/"
    outdir += "5TeV/" if is5TeV else "13TeV/"
    outdir += "iso_mu/" if doMuon else "iso_e/"

    if not is5TeV:
        if doMuon:
            input_antiiso_data = "inputs/awmunu/input_data.txt"
            input_antiiso_wl0 = "inputs/awmunu/input_wm0.txt"
            input_antiiso_wl1 = "inputs/awmunu/input_wm1.txt"
            input_antiiso_wl2 = "inputs/awmunu/input_wm2.txt"
            input_antiiso_ttbar = "inputs/awmunu/input_ttbar_dilepton.txt"
            input_antiiso_ttbar_1lep = "inputs/awmunu/input_ttbar_singlelepton.txt"
            input_antiiso_ttbar_0lep = "inputs/awmunu/input_ttbar_hadronic.txt"
            input_antiiso_ww = "inputs/awmunu/input_ww.txt"
            input_antiiso_wz = "inputs/awmunu/input_wz.txt"
            input_antiiso_zz = "inputs/awmunu/input_zz.txt"
            input_antiiso_zxx = "inputs/awmunu/input_zxx.txt"
            input_antiiso_wx0 = "inputs/awmunu/input_wx0.txt"
            input_antiiso_wx1 = "inputs/awmunu/input_wx1.txt"
            input_antiiso_wx2 = "inputs/awmunu/input_wx2.txt"
            input_data = "inputs/wmunu/input_data.txt"
            input_wl0 = "inputs/wmunu/input_wm0.txt"
            input_wl1 = "inputs/wmunu/input_wm1.txt"
            input_wl2 = "inputs/wmunu/input_wm2.txt"
            input_ttbar = "inputs/wmunu/input_ttbar_dilepton.txt"
            input_ttbar_1lep = "inputs/wmunu/input_ttbar_singlelepton.txt"
            input_ttbar_0lep = "inputs/wmunu/input_ttbar_hadronic.txt"
            input_ww = "inputs/wmunu/input_ww.txt"
            input_wz = "inputs/wmunu/input_wz.txt"
            input_zz = "inputs/wmunu/input_zz.txt"
            input_zxx = "inputs/wmunu/input_zxx.txt"
            input_wx0 = "inputs/wmunu/input_wx0.txt"
            input_wx1 = "inputs/wmunu/input_wx1.txt"
            input_wx2 = "inputs/wmunu/input_wx2.txt"
        else:
            input_antiiso_data = "inputs/awenu/input_data.txt"
            input_antiiso_wl0 = "inputs/awenu/input_we0.txt"
            input_antiiso_wl1 = "inputs/awenu/input_we1.txt"
            input_antiiso_wl2 = "inputs/awenu/input_we2.txt"
            input_antiiso_ttbar = "inputs/awenu/input_ttbar_dilepton.txt"
            input_antiiso_ttbar_1lep = "inputs/awenu/input_ttbar_singlelepton.txt"
            input_antiiso_ttbar_0lep = "inputs/awenu/input_ttbar_hadronic.txt"
            input_antiiso_ww = "inputs/awenu/input_ww.txt"
            input_antiiso_wz = "inputs/awenu/input_wz.txt"
            input_antiiso_zz = "inputs/awenu/input_zz.txt"
            input_antiiso_zxx = "inputs/awenu/input_zxx.txt"
            input_antiiso_wx0 = "inputs/awenu/input_wx0.txt"
            input_antiiso_wx1 = "inputs/awenu/input_wx1.txt"
            input_antiiso_wx2 = "inputs/awenu/input_wx2.txt"
            input_data = "inputs/wenu/input_data.txt"
            input_wl0 = "inputs/wenu/input_we0.txt"
            input_wl1 = "inputs/wenu/input_we1.txt"
            input_wl2 = "inputs/wenu/input_we2.txt"
            input_ttbar = "inputs/wenu/input_ttbar_dilepton.txt"
            input_ttbar_1lep = "inputs/wenu/input_ttbar_singlelepton.txt"
            input_ttbar_0lep = "inputs/wenu/input_ttbar_hadronic.txt"
            input_ww = "inputs/wenu/input_ww.txt"
            input_wz = "inputs/wenu/input_wz.txt"
            input_zz = "inputs/wenu/input_zz.txt"
            input_zxx = "inputs/wenu/input_zxx.txt"
            input_wx0 = "inputs/wenu/input_wx0.txt"
            input_wx1 = "inputs/wenu/input_wx1.txt"
            input_wx2 = "inputs/wenu/input_wx2.txt"

    else:
        # analyze 5TeV data
        if doMuon:
            input_antiiso_data = "inputs_5TeV/awmunu/input_data.txt"
            input_antiiso_wl = "inputs_5TeV/awmunu/input_wm.txt"
            input_antiiso_ttbar = "inputs_5TeV/awmunu/input_ttbar.txt"
            input_antiiso_ww = "inputs_5TeV/awmunu/input_ww.txt"
            input_antiiso_wz = "inputs_5TeV/awmunu/input_wz.txt"
            input_antiiso_zz2l = "inputs_5TeV/awmunu/input_zz2l.txt"
            input_antiiso_zz4l = "inputs_5TeV/awmunu/input_zz4l.txt"
            input_antiiso_zxx = "inputs_5TeV/awmunu/input_zxx.txt"
            input_antiiso_wx = "inputs_5TeV/awmunu/input_wx.txt"
            input_data = "inputs_5TeV/wmunu/input_data.txt"
            input_wl = "inputs_5TeV/wmunu/input_wm.txt"
            input_ttbar = "inputs_5TeV/wmunu/input_ttbar.txt"
            input_ww = "inputs_5TeV/wmunu/input_ww.txt"
            input_wz = "inputs_5TeV/wmunu/input_wz.txt"
            input_zz2l = "inputs_5TeV/wmunu/input_zz2l.txt"
            input_zz4l = "inputs_5TeV/wmunu/input_zz4l.txt"
            input_zxx = "inputs_5TeV/wmunu/input_zxx.txt"
            input_zxx2 = "inputs_5TeV/wmunu/input_zxx2.txt"
            input_wx = "inputs_5TeV/wmunu/input_wx.txt"
        else:
            input_antiiso_data = "inputs_5TeV/awenu/input_data.txt"
            input_antiiso_wl = "inputs_5TeV/awenu/input_we.txt"
            input_antiiso_ttbar = "inputs_5TeV/awenu/input_ttbar.txt"
            input_antiiso_ww = "inputs_5TeV/awenu/input_ww.txt"
            input_antiiso_wz = "inputs_5TeV/awenu/input_wz.txt"
            input_antiiso_zz2l = "inputs_5TeV/awenu/input_zz2l.txt"
            input_antiiso_zz4l = "inputs_5TeV/awenu/input_zz4l.txt"
            input_antiiso_zxx = "inputs_5TeV/awenu/input_zxx.txt"
            input_antiiso_wx = "inputs_5TeV/awenu/input_wx.txt"
            input_data = "inputs_5TeV/wenu/input_data.txt"
            input_wl = "inputs_5TeV/wenu/input_we.txt"
            input_ttbar = "inputs_5TeV/wenu/input_ttbar.txt"
            input_ww = "inputs_5TeV/wenu/input_ww.txt"
            input_wz = "inputs_5TeV/wenu/input_wz.txt"
            input_zz2l = "inputs_5TeV/wenu/input_zz2l.txt"
            input_zz4l = "inputs_5TeV/wenu/input_zz4l.txt"
            input_zxx = "inputs_5TeV/wenu/input_zxx.txt"
            input_wx = "inputs_5TeV/wenu/input_wx.txt"

    inputs_data = []
    inputs_data.append(input_antiiso_data)
    if not doTest:
        # test only runs on the antiiso data for quick testing of code
        inputs_data.append(input_data)

    DataSamp = Sample(inputs_data, isMC=False, name="Data",
                      isWSR=True, legend='QCD', color='226')
    if not is5TeV:
        qcdnorm = 1.0
        mcscale = 1.0
        Wl0AisoSamp = Sample(input_antiiso_wl0, isMC=True, name="wl0_aiso",
                             isWSR=True, additionalnorm=qcdnorm * mcscale, doTheoryVariation=False)
        Wl1AisoSamp = Sample(input_antiiso_wl1, isMC=True, name="wl1_aiso",
                             isWSR=True, additionalnorm=qcdnorm * mcscale, doTheoryVariation=False)
        Wl2AisoSamp = Sample(input_antiiso_wl2, isMC=True, name="wl2_aiso",
                             isWSR=True, additionalnorm=qcdnorm * mcscale, doTheoryVariation=False)
        # ttbar
        TTbarAisoSamp = Sample(input_antiiso_ttbar, isMC=True, name="ttbar_dilepton_aiso",
                               isWSR=True, additionalnorm=qcdnorm * mcscale, doTheoryVariation=False)
        TT1LepAisoSamp = Sample(input_antiiso_ttbar_1lep, isMC=True, name="ttbar_1lepton_aiso",
                                isWSR=True, additionalnorm=qcdnorm * mcscale, doTheoryVariation=False)
        TT0LepAisoSamp = Sample(input_antiiso_ttbar_0lep, isMC=True, name="ttbar_0lepton_aiso",
                                isWSR=True, additionalnorm=qcdnorm * mcscale, doTheoryVariation=False)
        # dibosons
        WWAisoSamp = Sample(input_antiiso_ww, isMC=True, name="WW_aiso",
                            isWSR=True, additionalnorm=qcdnorm * mcscale, doTheoryVariation=False)
        WZAisoSamp = Sample(input_antiiso_wz, isMC=True, name="WZ_aiso",
                            isWSR=True, additionalnorm=qcdnorm * mcscale, doTheoryVariation=False)
        ZZAisoSamp = Sample(input_antiiso_zz, isMC=True, name="ZZ_aiso",
                            isWSR=True, additionalnorm=qcdnorm * mcscale, doTheoryVariation=False)
        # tau
        ZXXAisoSamp = Sample(input_antiiso_zxx, isMC=True, name="ZXX_aiso",
                             isWSR=True, additionalnorm=qcdnorm * mcscale, doTheoryVariation=False)
        Wx0AisoSamp = Sample(input_antiiso_wx0, isMC=True, name="wx0_aiso",
                             isWSR=True, additionalnorm=qcdnorm * mcscale, doTheoryVariation=False)
        Wx1AisoSamp = Sample(input_antiiso_wx1, isMC=True, name="wx1_aiso",
                             isWSR=True, additionalnorm=qcdnorm * mcscale, doTheoryVariation=False)
        Wx2AisoSamp = Sample(input_antiiso_wx2, isMC=True, name="wx2_aiso",
                             isWSR=True, additionalnorm=qcdnorm * mcscale, doTheoryVariation=False)
        # W -> munu
        # wjetsnorm = 1.06
        wjetsnorm = 1.0
        Wl0Samp = Sample(input_wl0, isMC=True, name="wl0",
                         isWSR=True, additionalnorm=wjetsnorm, reweightZpt=reweightZpt)
        Wl1Samp = Sample(input_wl1, isMC=True, name="wl1",
                         isWSR=True, additionalnorm=wjetsnorm, reweightZpt=reweightZpt)
        Wl2Samp = Sample(input_wl2, isMC=True, name="wl2",
                         isWSR=True, additionalnorm=wjetsnorm, reweightZpt=reweightZpt)
        # ttbar
        TTbarSamp = Sample(input_ttbar, isMC=True,
                           name="ttbar_dilepton", isWSR=True)
        TT1LepSamp = Sample(input_ttbar_1lep, isMC=True,
                            name="ttbar_1lepton", isWSR=True)
        TT0LepSamp = Sample(input_ttbar_0lep, isMC=True,
                            name="ttbar_0lepton", isWSR=True)
        # dibosons
        WWSamp = Sample(input_ww, isMC=True, name="WW", isWSR=True)
        WZSamp = Sample(input_wz, isMC=True, name="WZ", isWSR=True)
        ZZSamp = Sample(input_zz, isMC=True, name="ZZ", isWSR=True)
        # tau
        ZXXSamp = Sample(input_zxx, isMC=True, name="ZXX",
                         isWSR=True, reweightZpt=reweightZpt)
        Wx0Samp = Sample(input_wx0, isMC=True, name="wx0",
                         isWSR=True, reweightZpt=reweightZpt)
        Wx1Samp = Sample(input_wx1, isMC=True, name="wx1",
                         isWSR=True, reweightZpt=reweightZpt)
        Wx2Samp = Sample(input_wx2, isMC=True, name="wx2",
                         isWSR=True, reweightZpt=reweightZpt)

        if not doTest:
            sampMan = SampleManager(DataSamp, [Wl0Samp, Wl1Samp, Wl2Samp, TTbarSamp, TT1LepSamp,
                                               TT0LepSamp, WWSamp, WZSamp, ZZSamp, ZXXSamp, Wx0Samp, Wx1Samp, Wx2Samp, Wl0AisoSamp, Wl1AisoSamp, Wl2AisoSamp, TTbarAisoSamp, TT1LepAisoSamp,
                                               TT0LepAisoSamp, WWAisoSamp, WZAisoSamp, ZZAisoSamp, ZXXAisoSamp, Wx0AisoSamp, Wx1AisoSamp, Wx2AisoSamp])
        else:
            sampMan = SampleManager(DataSamp, [TTbarSamp, TT1LepSamp])
        sampMan.groupMCs(["wx0", "wx1", "wx2", "wx0_aiso",
                         "wx1_aiso", "wx2_aiso"], "wx", 216, 'wx')
        sampMan.groupMCs(['ZXX', 'ZXX_aiso'], 'zxx', 216, 'zxx')
        sampMan.groupMCs(["WW_aiso", "WZ_aiso", "ZZ_aiso", "WW", "WZ", "ZZ"], "VV", 216, "VV")
        sampMan.groupMCs(["ttbar_dilepton_aiso", "ttbar_1lepton_aiso",
                         "ttbar_0lepton_aiso", "ttbar_dilepton", "ttbar_1lepton",
                          "ttbar_0lepton"], "ttbar", 96, "t#bar{t}")
        label = "W#rightarrow#mu#nu" if doMuon else "W#rightarrow e#nu"
        sampMan.groupMCs(
            ['wl0_aiso', 'wl1_aiso', 'wl2_aiso', 'wl0', 'wl1', 'wl2'], "wlnu", 92, label)

    else:
        qcdnorm = 1.0
        mcscale = 1.0
        # 5 TeV
        # W -> lnu
        label = "W#rightarrow#mu#nu" if doMuon else "W#rightarrow e#nu"
        WlAisoSamp = Sample(input_antiiso_wl, isMC=True, name="wlnu_aiso", isWSR=True, additionalnorm=qcdnorm *
                            mcscale, is5TeV=True, color=92, legend=label, doTheoryVariation=False)
        # ttbar
        TTbarAisoSamp = Sample(input_antiiso_ttbar, isMC=True, name="ttbar_aiso",     isWSR=True,
                               additionalnorm=qcdnorm * mcscale, is5TeV=True, color=96, legend="t#bar{t}", doTheoryVariation=False)
        # dibosons
        WWAisoSamp = Sample(input_antiiso_ww, isMC=True, name="WW_aiso", isWSR=True,
                            additionalnorm=qcdnorm * mcscale, is5TeV=True, doTheoryVariation=False)
        WZAisoSamp = Sample(input_antiiso_wz, isMC=True, name="WZ_aiso", isWSR=True,
                            additionalnorm=qcdnorm * mcscale, is5TeV=True, doTheoryVariation=False)
        ZZ2LAisoSamp = Sample(input_antiiso_zz2l, isMC=True, name="ZZ2L_aiso", isWSR=True,
                              additionalnorm=qcdnorm * mcscale, is5TeV=True, doTheoryVariation=False)
        ZZ4LAisoSamp = Sample(input_antiiso_zz4l, isMC=True, name="ZZ4L_aiso", isWSR=True,
                              additionalnorm=qcdnorm * mcscale, is5TeV=True, doTheoryVariation=False)
        # tau
        ZXXAisoSamp = Sample(input_antiiso_zxx, isMC=True, name="ZXX_aiso", isWSR=True,
                             additionalnorm=qcdnorm * mcscale, is5TeV=True, doTheoryVariation=False)
        WxAisoSamp = Sample(input_antiiso_wx,  isMC=True, name="wx_aiso",  isWSR=True,
                            additionalnorm=qcdnorm * mcscale, is5TeV=True, doTheoryVariation=False)

        wjetsnorm = 1.0
        WlSamp = Sample(input_wl, isMC=True, name="wlnu", color=92,
                        legend="W#rightarrow#mu#nu", isWSR=True, additionalnorm=wjetsnorm, is5TeV=True, reweightZpt=reweightZpt)
        # ttbar
        TTbarSamp = Sample(input_ttbar, isMC=True, name="ttbar",
                           color=86, legend="t#bar{t}", isWSR=True, is5TeV=True)
        # dibosons
        WWSamp = Sample(input_ww, isMC=True, name="WW",
                        isWSR=True, is5TeV=True)
        WZSamp = Sample(input_wz, isMC=True, name="WZ",
                        isWSR=True, is5TeV=True)
        ZZ2LSamp = Sample(input_zz2l, isMC=True, name="ZZ2L",
                          isWSR=True, is5TeV=True)
        ZZ4LSamp = Sample(input_zz4l, isMC=True, name="ZZ4L",
                          isWSR=True, is5TeV=True)
        # tau
        ZXXSamp = Sample(input_zxx, isMC=True, name="ZXX",
                         isWSR=True, is5TeV=True, reweightZpt=reweightZpt)
        #ZXX2Samp = Sample(input_zxx2, isMC=True, name="ZXX2",
        #                  isWSR=True, is5TeV=True, reweightZpt=reweightZpt)
        WxSamp = Sample(input_wx, isMC=True, name="WX",
                        color=216, legend="WX", isWSR=True, is5TeV=True, reweightZpt=reweightZpt)

        sampMan = SampleManager(DataSamp, [WlAisoSamp, TTbarAisoSamp, WWAisoSamp, WZAisoSamp, ZZ2LAisoSamp, ZZ4LAisoSamp,
                                ZXXAisoSamp, WxAisoSamp, WlSamp, TTbarSamp, WWSamp, WZSamp, ZZ2LSamp, ZZ4LSamp, ZXXSamp, WxSamp], is5TeV=True)
        sampMan.groupMCs(["WW", "WZ", "ZZ2L", "ZZ4L", "WW_aiso",
                         "WZ_aiso", "ZZ2L_aiso", "ZZ4L_aiso"], 'VV', 216, 'VV')
        sampMan.groupMCs(['ZXX', 'ZXX2', 'ZXX_aiso'], 'zxx', 216, 'zxx')
        sampMan.groupMCs(['WX', 'wx_aiso'], 'wx', 216, 'wx')
        sampMan.groupMCs(['ttbar', 'ttbar_aiso'], 'ttbar', 86, 't#bar{t}')
        label = "W#rightarrow#mu#nu" if doMuon else "W#rightarrow e#nu"
        sampMan.groupMCs(['wlnu_aiso', 'wlnu'], "wl", 92, label)

    sampMan.outdir = outdir

    sampMan.DefineAll("Lep_pt", "lep.Pt()")
    sampMan.ApplyCutAll("Lep_pt > 25.0")
    sampMan.DefineAll("Lep_eta", "lep.Eta()")
    sampMan.ApplyCutAll("fabs(Lep_eta) < 2.4")

    # muon and electron isolation distributions are different
    # more coarse binning for electrons to make sure enough statistics
    if doMuon:
        #isoCuts = [0.0,0.001, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50,
        #           0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95]
        isoCuts = [0.0,0.001, 0.15, 0.25, 0.50, 0.75, 1.0]
        isobins = []
        sampMan.DefineAll("RelIso", isovar)
        for isobin in range(0, 6):
            sampMan.DefineAll(
                f"w_iso{isobin}", f"(RelIso >= {isoCuts[isobin]} && RelIso < {isoCuts[isobin+1]})")
            isobins.append(f"iso{isobin}")

    else:
        isoCuts = [0.0, 0.10, 0.15, 0.20, 0.25, 0.40, 0.60, 0.90]
        isobins = []
        sampMan.DefineAll("isEB",   "fabs(Lep_eta) <= 1.4442")
        # sampMan.DefineAll("RelIso", "isEB ? (relIso + 0.0287 - 0.0478) : (relIso + 0.0445 - 0.0658)")
        #sampMan.DefineAll("RelIso", "(pfCombIso/lep.Pt())")
        sampMan.DefineAll("RelIso", isovar)
        for isobin in range(3, 10):
            sampMan.DefineAll(
                f"w_iso{isobin}", f"(RelIso > {isoCuts[isobin-3]} && RelIso < {isoCuts[isobin-2]})")
            isobins.append(f"iso{isobin}")

    sampMan.DefineMC("met_pt", "metVars[1]")
    sampMan.DefineMC("met_phi", "TVector2::Phi_mpi_pi(metVarsPhi[1])")
    DataSamp.Define("met_pt", "metVars[0]")
    DataSamp.Define("met_phi", "TVector2::Phi_mpi_pi(metVarsPhi[0])")

    # recoil variables
    sampMan.DefineAll("V2W", "UVec(lep.Pt(), lep.Phi(), met_pt, met_phi)")
    sampMan.DefineAll("WpT", "V2W.Mod()")
    sampMan.DefineAll("Wphi", "TVector2::Phi_mpi_pi(V2W.Phi())")
    
    sampMan.DefineAll("deltaPhi", "abs(TVector2::Phi_mpi_pi(met_phi - lep.Phi()))")

    # charge
    lepname = "mu" if doMuon else "e"
    sampMan.DefineAll(lepname+"plus",  "q > 0")
    sampMan.DefineAll(lepname+"minus", "q < 0")
    chgbins = [lepname+"plus", lepname + "minus"]

    # WpT bins
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
        wptbins = ["WpT_bin0", "WpT_bin1", "WpT_bin2", "WpT_bin3", "WpT_bin4",
                   "WpT_bin5", "WpT_bin6", "WpT_bin7", "WpT_bin8", "WpT_bin9"]

    # eta bins for electrons: barral and endcap
    sampMan.DefineAll("lepEta_bin0", "1.0")
    sampMan.DefineAll("lepEta_bin1", "abs(lep.Eta()) <= 1.4442")
    sampMan.DefineAll("lepEta_bin2", "abs(lep.Eta()) > 1.4442")

    if doMuon:
        etabins = ["lepEta_bin0"]
    else:
        etabins = ["lepEta_bin0", "lepEta_bin1", "lepEta_bin2"]
    etabins = ["lepEta_bin0"]
    
    mass_bins = mass_bins_test[0]
    mtbins = []
    for imt in range(len(mass_bins)-1):
        sampMan.DefineAll(
            f"mt{imt}", f"mtCorr >= {mass_bins[imt]} && mtCorr < {mass_bins[imt+1]}")
        mtbins.append(f"mt{imt}")
        
    def makeDraws(strname, ymax):
        ptbins = np.array([25.0, 30.0, 35.0, 45.0, 60.0, 80.0])
        metbins = np.array([0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 80.0])
        etabins = np.array([-2.4, -2.0, -1.6, -1.2, -0.8, -0.4, 0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4])
        phibins = np.array([0, 0.8, 1.6, 2.4, 3.2])
        
        sampMan.cacheDraw("RelIso", f"histo_wjets_{lepname}_RelIso_{strname}", 72, 0, 0.72, DrawConfig(xmin=0., xmax=0.75, xlabel="Relative Isolation", ylabel=f"Events / bin", dology=True, ymax=ymax, donormalizebin=False, addOverflow=True, addUnderflow=True, showratio=True, legendPos=[0.94, 0.88, 0.70, 0.68]), weightname=strname)
        
        sampMan.cacheDraw("Lep_pt", f"histo_wjets_{lepname}_Lep_pt_{strname}", ptbins, DrawConfig(xmin=20, xmax=80, xlabel="Lepton p_{T} [GeV]", ylabel=f"Events / bin", dology=True, ymax=ymax, donormalizebin=False, addOverflow=True, addUnderflow=True, showratio=True, legendPos=[0.94, 0.88, 0.70, 0.68]), weightname=strname)
        
        sampMan.cacheDraw("Lep_eta", f"histo_wjets_{lepname}_Lep_eta_{strname}", etabins, DrawConfig(xmin=-2.5, xmax=2.5, xlabel="Lepton #eta", ylabel=f"Events / bin", dology=True, ymax=ymax, donormalizebin=False, addOverflow=True, addUnderflow=True, showratio=True, legendPos=[0.94, 0.88, 0.70, 0.68]), weightname=strname)
        
        sampMan.cacheDraw("met_pt", f"histo_wjets_{lepname}_met_pt_{strname}", metbins, DrawConfig(xmin=0, xmax=80, xlabel="MET [GeV]", ylabel=f"Events / bin", dology=True, ymax=ymax, donormalizebin=False, addOverflow=True, addUnderflow=True, showratio=True, legendPos=[0.94, 0.88, 0.70, 0.68]), weightname=strname)
        
        sampMan.cacheDraw("deltaPhi", f"histo_wjets_{lepname}_deltaPhi_{strname}", phibins, DrawConfig(xmin=0, xmax=3.2, xlabel="#Delta#phi(#mu, MET)", ylabel=f"Events / bin", dology=True, ymax=ymax, donormalizebin=False, addOverflow=True, addUnderflow=True, showratio=True, legendPos=[0.94, 0.88, 0.70, 0.68]), weightname=strname)
        
        sampMan.cacheDraw2D("Lep_pt", "Lep_eta", f"histo_wjets_{lepname}_pt_vs_eta_{strname}", ptbins, etabins, DrawConfig(xmin=25, xmax=80, ymin=-2.4, ymax=2.4, xlabel="Lepton p_{T} [GeV]", ylabel="Lepton #eta", zlabel="Events / bin", dology=False, donormalizebin=False, addOverflow=True, addUnderflow=True, showratio=False, legendPos=[0.94, 0.88, 0.70, 0.68]), weightname=strname)
        sampMan.cacheDraw2D("Lep_pt", "met_pt", f"histo_wjets_{lepname}_pt_vs_met_{strname}", ptbins, metbins, DrawConfig(xmin=25, xmax=80, ymin=0, ymax=80, xlabel="Lepton p_{T} [GeV]", ylabel="MET [GeV]", zlabel="Events / bin", dology=False, donormalizebin=False, addOverflow=True, addUnderflow=True, showratio=False, legendPos=[0.94, 0.88, 0.70, 0.68]), weightname=strname)
        sampMan.cacheDraw2D("Lep_pt", "deltaPhi", f"histo_wjets_{lepname}_pt_vs_deltaPhi_{strname}", ptbins, phibins, DrawConfig(xmin=25, xmax=80, ymin=0, ymax=3.2, xlabel="Lepton p_{T} [GeV]", ylabel="#Delta#phi(#mu, MET)", zlabel="Events / bin", dology=False, donormalizebin=False, addOverflow=True, addUnderflow=True, showratio=False, legendPos=[0.94, 0.88, 0.70, 0.68]), weightname=strname)
        
        # make the pt vs eta and pt vs met plots in deltaPhi < pi/2 and deltaPhi > pi/2
        sampMan.DefineAll(f"{strname}_deltaPhiP", f"{strname} * (deltaPhi >= {ROOT.TMath.Pi()/2.})")
        sampMan.DefineAll(f"{strname}_deltaPhiM", f"{strname} * (deltaPhi <  {ROOT.TMath.Pi()/2.})")
        sampMan.cacheDraw2D("Lep_pt", "Lep_eta", f"histo_wjets_{lepname}_pt_vs_eta_{strname}_deltaPhiP", ptbins, etabins, DrawConfig(xmin=25, xmax=80, ymin=-2.4, ymax=2.4, xlabel="Lepton p_{T} [GeV]", ylabel="Lepton #eta", zlabel="Events / bin", dology=False, donormalizebin=False, addOverflow=True, addUnderflow=True, showratio=False, legendPos=[0.94, 0.88, 0.70, 0.68]), weightname=f"{strname}_deltaPhiP")
        sampMan.cacheDraw2D("Lep_pt", "met_pt", f"histo_wjets_{lepname}_pt_vs_met_{strname}_deltaPhiP", ptbins, metbins, DrawConfig(xmin=25, xmax=80, ymin=0, ymax=80, xlabel="Lepton p_{T} [GeV]", ylabel="MET [GeV]", zlabel="Events / bin", dology=False, donormalizebin=False, addOverflow=True, addUnderflow=True, showratio=False, legendPos=[0.94, 0.88, 0.70, 0.68]), weightname=f"{strname}_deltaPhiP")
        sampMan.cacheDraw2D("Lep_pt", "Lep_eta", f"histo_wjets_{lepname}_pt_vs_eta_{strname}_deltaPhiM", ptbins, etabins, DrawConfig(xmin=25, xmax=80, ymin=-2.4, ymax=2.4, xlabel="Lepton p_{T} [GeV]", ylabel="Lepton #eta", zlabel="Events / bin", dology=False, donormalizebin=False, addOverflow=True, addUnderflow=True, showratio=False, legendPos=[0.94, 0.88, 0.70, 0.68]), weightname=f"{strname}_deltaPhiM")
        sampMan.cacheDraw2D("Lep_pt", "met_pt", f"histo_wjets_{lepname}_pt_vs_met_{strname}_deltaPhiM", ptbins, metbins, DrawConfig(xmin=25, xmax=80, ymin=0, ymax=80, xlabel="Lepton p_{T} [GeV]", ylabel="MET [GeV]", zlabel="Events / bin", dology=False, donormalizebin=False, addOverflow=True, addUnderflow=True, showratio=False, legendPos=[0.94, 0.88, 0.70, 0.68]), weightname=f"{strname}_deltaPhiM")

    # draw the lepton isolation distribution
    ymax = 5e7 if doMuon else 1e6
    #for wpt in wptbins:
    #    for lepeta in etabins:
    #        strname = f"weight_{wpt}_{lepeta}"
    #        sampMan.DefineAll(strname, f"weight_WVpt * {wpt} * {lepeta}")
    #        makeDraws(strname, ymax)
    #        
    #        for mt in mtbins:
    #            strname = f"weight_{wpt}_{lepeta}_{mt}"
    #            sampMan.DefineAll(strname, f"weight_WVpt * {wpt} * {lepeta} * {mt}")
    #            makeDraws(strname, ymax)
                
    idx = 0
    # do the fine binning first; then rebin in the processHists
    #mass_bins = np.array([0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0,
    #                     70.0, 75.0, 80.0, 85.0, 90.0, 95.0, 100.0, 105.0, 110.0, 115.0, 120.0, 125.0, 130.0, 135.0, 140.])
    for iso in isobins:
        for wpt in wptbins:
            for lepeta in etabins:
                for chg in chgbins:
                    strname = "weight_{}_{}_{}_{}".format(
                        chg, iso, wpt, lepeta)

                    sampMan.DefineAll(
                        strname, "w_{} * weight_WVpt * {} * {} * {}".format(iso, wpt, lepeta, chg))
                   
                    sampMan.cacheDraw("mtCorr", f"histo_wjets_{lepname}_mtcorr_{strname}", mass_bins, DrawConfig(xmin=0., xmax=140.0, xlabel="m_{T} [GeV]", ylabel=f"Events / 5 GeV", dology=True, ymax=6e5, donormalizebin=False, addOverflow=True, addUnderflow=True, showratio=False, lheader=f"{isoCuts[idx]} < I < {isoCuts[idx + 1]}", legendPos=[0.94, 0.88, 0.70, 0.68]), weightname=strname)

                    makeDraws(strname, ymax)
                    
                    for mt in mtbins:
                        strname = "weight_{}_{}_{}_{}_{}".format(chg, iso, wpt, lepeta, mt)
                        sampMan.DefineAll(strname, "w_{} * weight_WVpt * {} * {} * {} * {}".format(iso, wpt, lepeta, chg, mt))
                        makeDraws(strname, ymax)
                        
        idx += 1

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
    
    hpt_vs_etas_deltaPhiP_data = OrderedDict()
    hpt_vs_etas_deltaPhiP_mc = OrderedDict()
    hpt_vs_etas_deltaPhiP_data_mt = OrderedDict()
    hpt_vs_etas_deltaPhiP_mc_mt = OrderedDict()
    
    hpt_vs_etas_deltaPhiM_data = OrderedDict()
    hpt_vs_etas_deltaPhiM_mc = OrderedDict()
    hpt_vs_etas_deltaPhiM_data_mt = OrderedDict()
    hpt_vs_etas_deltaPhiM_mc_mt = OrderedDict() 
    
    hpt_vs_mets_data = OrderedDict()
    hpt_vs_mets_mc = OrderedDict()
    hpt_vs_mets_data_mt = OrderedDict()
    hpt_vs_mets_mc_mt = OrderedDict()
    
    hpt_vs_mets_deltaPhiP_data = OrderedDict()
    hpt_vs_mets_deltaPhiP_mc = OrderedDict()
    hpt_vs_mets_deltaPhiP_data_mt = OrderedDict()
    hpt_vs_mets_deltaPhiP_mc_mt = OrderedDict() 
    
    hpt_vs_mets_deltaPhiM_data = OrderedDict()
    hpt_vs_mets_deltaPhiM_mc = OrderedDict()
    hpt_vs_mets_deltaPhiM_data_mt = OrderedDict()
    hpt_vs_mets_deltaPhiM_mc_mt = OrderedDict()  
    
    hpt_vs_deltaphis_data = OrderedDict()
    hpt_vs_deltaphis_mc = OrderedDict()
    hpt_vs_deltaphis_data_mt = OrderedDict()
    hpt_vs_deltaphis_mc_mt = OrderedDict()
    
    def getHisto(name):
        hdata = sampMan.hdatas[name]
        hstacked = THStack2TH1(sampMan.hsmcs[name])
        for ibin in range(hstacked.GetNbinsX()+1):
            hstacked.SetBinContent(ibin, max(hstacked.GetBinContent(ibin), 0))
        hdata.Add(hstacked, -1.0)
        hdata.SetName(name)
        return hdata
    
    def getHisto2D(name):
        """
        return both data and merged mc
        """
        hdata = sampMan.hdatas2D[name]
        hmc = sampMan.hmcs2D[name]
        hdata.SetName(name + "_Data")
        hmc.SetName(name + "_MC")
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
                   
                    outputname = f"histo_wjets_{lepname}_pt_vs_eta_{strname}_deltaPhiP"
                    hpt_vs_etas_deltaPhiP_data[strname], hpt_vs_etas_deltaPhiP_mc[strname] = getHisto2D(outputname)
                    
                    outputname = f"histo_wjets_{lepname}_pt_vs_met_{strname}_deltaPhiP"
                    hpt_vs_mets_deltaPhiP_data[strname], hpt_vs_mets_deltaPhiP_mc[strname] = getHisto2D(outputname) 
                    
                    outputname = f"histo_wjets_{lepname}_pt_vs_eta_{strname}_deltaPhiM"
                    hpt_vs_etas_deltaPhiM_data[strname], hpt_vs_etas_deltaPhiM_mc[strname] = getHisto2D(outputname)
                    
                    outputname = f"histo_wjets_{lepname}_pt_vs_met_{strname}_deltaPhiM"
                    hpt_vs_mets_deltaPhiM_data[strname], hpt_vs_mets_deltaPhiM_mc[strname] = getHisto2D(outputname)  
                    
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
                        
                        outputname = f"histo_wjets_{lepname}_pt_vs_eta_{strname}_deltaPhiP"
                        hpt_vs_etas_deltaPhiP_data_mt[strname], hpt_vs_etas_deltaPhiP_mc_mt[strname] = getHisto2D(outputname)
                        
                        outputname = f"histo_wjets_{lepname}_pt_vs_met_{strname}_deltaPhiP"
                        hpt_vs_mets_deltaPhiP_data_mt[strname], hpt_vs_mets_deltaPhiP_mc_mt[strname] = getHisto2D(outputname) 
                        
                        outputname = f"histo_wjets_{lepname}_pt_vs_eta_{strname}_deltaPhiM"
                        hpt_vs_etas_deltaPhiM_data_mt[strname], hpt_vs_etas_deltaPhiM_mc_mt[strname] = getHisto2D(outputname)
                        
                        outputname = f"histo_wjets_{lepname}_pt_vs_met_{strname}_deltaPhiM"
                        hpt_vs_mets_deltaPhiM_data_mt[strname], hpt_vs_mets_deltaPhiM_mc_mt[strname] = getHisto2D(outputname)  
                        
                        outputname = f"histo_wjets_{lepname}_pt_vs_deltaPhi_{strname}"
                        hpt_vs_deltaphis_data_mt[strname], hpt_vs_deltaphis_mc_mt[strname] = getHisto2D(outputname)
                        

    suffix += "_" + lepname + "nu"
    sqrtS = "5TeV" if is5TeV else "13TeV"
    suffix += f"_{sqrtS}.root"

    sampMan.dumpCounts()

    dumpHistos([hIsos, hIsos_mt], f"output_qcdIsoMean_{suffix}", doPrint=True)
    dumpHistos([hleppts, hleppts_mt], f"output_qcdLepPtMean_{suffix}", doPrint=False)
    dumpHistos([hlepetas, hlepetas_mt], f"output_qcdLepEtaMean_{suffix}", doPrint=False)
    dumpHistos([hmetpts, hmetpts_mt], f"output_qcdMetPtMean_{suffix}", doPrint=False)
    dumpHistos([hdeltaphis, hdeltaphis_mt], f"output_qcdDeltaPhiMean_{suffix}", doPrint=False)
    
    dumpHistos([hpt_vs_etas_data, hpt_vs_etas_mc, hpt_vs_etas_data_mt, hpt_vs_etas_mc_mt, hpt_vs_etas_deltaPhiP_data, hpt_vs_etas_deltaPhiP_mc, hpt_vs_etas_deltaPhiP_data_mt, hpt_vs_etas_deltaPhiP_mc_mt, hpt_vs_etas_deltaPhiM_data, hpt_vs_etas_deltaPhiM_mc, hpt_vs_etas_deltaPhiM_data_mt, hpt_vs_etas_deltaPhiM_mc_mt], f"output_qcdLepPtVsEtaMean_{suffix}", doPrint=False)
    
    dumpHistos([hpt_vs_mets_data, hpt_vs_mets_mc, hpt_vs_mets_data_mt, hpt_vs_mets_mc_mt, hpt_vs_mets_deltaPhiP_data, hpt_vs_mets_deltaPhiP_mc, hpt_vs_mets_deltaPhiP_data_mt, hpt_vs_mets_deltaPhiP_mc_mt, hpt_vs_mets_deltaPhiM_data, hpt_vs_mets_deltaPhiM_mc, hpt_vs_mets_deltaPhiM_data_mt, hpt_vs_mets_deltaPhiM_mc_mt], f"output_qcdLepPtVsMetMean_{suffix}", doPrint=False)
    
    dumpHistos([hpt_vs_deltaphis_data, hpt_vs_deltaphis_mc, hpt_vs_deltaphis_data_mt, hpt_vs_deltaphis_mc_mt], f"output_qcdLepPtVsDeltaPhiMean_{suffix}", doPrint=False)
    
    #
    # for the qcd background extrapolation
    # 
    hmts_comp = OrderedDict()
    # hetas_mtCut_comp = OrderedDict()
    for iso in isobins:
        for wpt in wptbins:
            for lepeta in etabins:
                for chg in chgbins:
                    strname = "weight_{}_{}_{}_{}".format(
                        chg, iso, wpt, lepeta)

                    # for mT
                    outputname = f"histo_wjets_{lepname}_mtcorr_" + strname
                    hstacked = THStack2TH1(sampMan.hsmcs[outputname])
                    hdata = sampMan.hdatas[outputname]
                    for ibin in range(hstacked.GetNbinsX()+1):
                        # hstacked should always be above 0
                        val_mc = hstacked.GetBinContent(ibin)
                        if val_mc < 0:
                            hstacked.SetBinContent(ibin, 0)
                            hstacked.SetBinError(ibin, 0)
                        # subtract mc from data
                        val_data = hdata.GetBinContent(ibin)
                        hdata.SetBinContent(ibin, max(val_data - val_mc, 0))
                        err_data = hdata.GetBinError(ibin)
                        err_mc = hstacked.GetBinError(ibin)
                        # apply 30% signal contamination unc
                        err_mcunc = val_mc * 0.3
                        hdata.SetBinError(ibin, np.sqrt(
                            err_data**2 + err_mc**2 + err_mcunc**2))

                    hmts_comp[strname] = hdata
                    hmts_comp[strname].SetName(outputname)

    outfile = ROOT.TFile.Open("root/output_qcdshape_backup_"+suffix, "recreate")

    for wpt in wptbins:
        # odir = outfile.mkdir(wpt)
        # outfile.cd(wpt)
        for iso in isobins:
            # skip the last iso bin as it is used for the uncertaintiy of the previous iso bin
            for lepeta in etabins:
                for chg in chgbins:
                    i = int(iso[3:])
                    if iso == isobins[-1]:
                        # a bit cheating,
                        # for the last bin, use the previous bin shape as the shape variation
                        iso_next = "iso" + str(i-1)
                    else:
                        iso_next = "iso" + str(i+1)
                    strname = "weight_{}_{}_{}_{}".format(
                        chg, iso, wpt, lepeta)
                    strname_next = "weight_{}_{}_{}_{}".format(
                        chg, iso_next, wpt, lepeta)

                    outputname = f"histo_wjets_{lepname}_mtcorr_" + strname

                    hcenter = hmts_comp[strname]
                    # shape uncertaintis are using the shape difference from the neighboring bins
                    # used in the original qcd background predictions
                    hup = hmts_comp[strname_next].Clone(outputname+"_shapeUp")
                    hdown = hmts_comp[strname_next].Clone(
                        outputname+"_shapeDown")

                    hup.Scale(hcenter.Integral() / (hup.Integral()+1e-6))

                    for ibin in range(1, hcenter.GetNbinsX()+1):
                        center = hcenter.GetBinContent(ibin)
                        up = hup.GetBinContent(ibin)
                        hdown.SetBinContent(ibin, max(2*center - up, 0))

                        hcenter.SetBinContent(ibin, max(center, 0))
                        hup.SetBinContent(ibin, max(up, 0))

                    # hcenter.SetDirectory(odir)
                    hcenter.SetDirectory(outfile)
                    hcenter.Write()
                    # hup.SetDirectory(odir)
                    hup.SetDirectory(outfile)
                    hup.Write()
                    # hdown.SetDirectory(odir)
                    hdown.SetDirectory(outfile)
                    hdown.Write()

    outfile.Close()

    sampMan.dumpCounts()

    print("Program end...")

    input()

    return


if __name__ == "__main__":
    main()

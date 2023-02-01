from xml import dom
import ROOT
import numpy as np
from collections import OrderedDict
import sys
from CMSPLOTS.myFunction import DrawHistos, THStack2TH1
import CMSPLOTS.CMS_lumi
import pickle
from modules.SampleManager import DrawConfig, Sample, SampleManager
from modules.Binnings import Vptbins
from modules.theoryUncIndices import theoryUnc_13TeV_zjets, theoryUnc_5TeV_zjets
import argparse

ROOT.gROOT.SetBatch(True)
ROOT.ROOT.EnableImplicitMT(10)


def main():
    print("Program start...")

    parser = argparse.ArgumentParser(
        description="Make plots for Wlnu analysis")
    parser.add_argument("--doTest", action="store_true", dest="doTest",
                        help="Run on a subset of samples for debugging; false runs on all dataset.")
    parser.add_argument("--is5TeV", action="store_true", dest="is5TeV",
                        help="Analyze the 5TeV data; false runs on 13TeV data")
    parser.add_argument("--doElectron", action="store_true", dest="doElectron",
                        help="Analyze the electron channel; false runs the muon channel")
    parser.add_argument("--doTheoryNorm", action="store_true", dest="doTheoryNorm",
                        help="Normalize the theory uncertainties to the nominal; false does not normalize")
    args = parser.parse_args()

    doTest = args.doTest
    is5TeV = args.is5TeV
    doMuon = not args.doElectron
    doTheoryNorm = args.doTheoryNorm

    print("doMuon: ", doMuon)
    print("doTest:", doTest)
    print("is5TeV:", is5TeV)
    print("doTheoryNorm:", doTheoryNorm)

    # ROOT.gROOT.ProcessLine('TFile* f_zpt = TFile::Open("data/zpt_ratio_amc2data.root")')
    # ROOT.gROOT.ProcessLine('TH1D* h_zpt_ratio  = (TH1D*)f_zpt->Get("hc")')
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

    if not is5TeV:
        if doMuon:
            if doTest:
                input_data = "inputs/zmumu/input_data.txt"
                input_ttbar = "inputs/zmumu/input_ttbar_dilepton.txt"
            else:
                input_data = "inputs/zmumu/input_data.txt"
                input_dy = "inputs/zmumu/input_zjets.txt"
                input_ttbar = "inputs/zmumu/input_ttbar_dilepton.txt"
                input_ttbar_1lep = "inputs/zmumu/input_ttbar_singlelepton.txt"
                input_ttbar_0lep = "inputs/zmumu/input_ttbar_hadronic.txt"
                input_ww = "inputs/zmumu/input_ww.txt"
                input_wz = "inputs/zmumu/input_wz.txt"
                input_zz = "inputs/zmumu/input_zz.txt"
                input_zxx = "inputs/zmumu/input_zxx.txt"
        else:
            # electron channel
            input_data = "inputs/zee/input_data.txt"
            input_dy = "inputs/zee/input_zjets.txt"
            input_ttbar = "inputs/zee/input_ttbar_dilepton.txt"
            input_ttbar_1lep = "inputs/zee/input_ttbar_singlelepton.txt"
            input_ttbar_0lep = "inputs/zee/input_ttbar_hadronic.txt"
            input_ww = "inputs/zee/input_ww.txt"
            input_wz = "inputs/zee/input_wz.txt"
            input_zz = "inputs/zee/input_zz.txt"
            input_zxx = "inputs/zee/input_zxx.txt"
    else:
        if doMuon:
            input_data = "inputs_5TeV/zmumu/input_data.txt"
            input_dy = "inputs_5TeV/zmumu/input_zjets.txt"
            input_ttbar = "inputs_5TeV/zmumu/input_ttbar.txt"
            input_ww = "inputs_5TeV/zmumu/input_ww.txt"
            input_wz = "inputs_5TeV/zmumu/input_wz.txt"
            input_zz2l = "inputs_5TeV/zmumu/input_zz2l.txt"
            input_zz4l = "inputs_5TeV/zmumu/input_zz4l.txt"
            input_zxx = "inputs_5TeV/zmumu/input_zxx.txt"
        else:
            input_data = "inputs_5TeV/zee/input_data.txt"
            input_dy = "inputs_5TeV/zee/input_zjets.txt"
            input_ttbar = "inputs_5TeV/zee/input_ttbar.txt"
            input_ww = "inputs_5TeV/zee/input_ww.txt"
            input_wz = "inputs_5TeV/zee/input_wz.txt"
            input_zz2l = "inputs_5TeV/zee/input_zz2l.txt"
            input_zz4l = "inputs_5TeV/zee/input_zz4l.txt"
            input_zxx = "inputs_5TeV/zee/input_zxx.txt"

    DataSamp = Sample(input_data, isMC=False, legend="Data",
                      name="Data", doTheoryVariation=False)
    lepchannel = "Z#rightarrow#mu^{+}#mu^{-}" if doMuon else "Z#rightarrow e^{+}e^{-}"
    if not is5TeV:
        TTbarSamp = Sample(input_ttbar, color=46,
                           legend="t#bar{t}", name="ttbar2lep")
        if not doTest:
            DYSamp = Sample(input_dy,    color=92,
                            legend=lepchannel, name="DY", reweightZpt=True)
            WWSamp = Sample(input_ww,  color=38,  legend="WW2L",     name="WW")
            WZSamp = Sample(input_wz,  color=39,  legend="WZ3L",     name="WZ")
            ZZSamp = Sample(input_zz,  color=37,  legend="ZZ2L",     name="ZZ")
            ZXXSamp = Sample(input_zxx, color=40,
                             legend="ZXX",      name="ZXX", reweightZpt=True)
            TT1LepSamp = Sample(input_ttbar_1lep,  color=47,
                                legend="t#bar{t}", name="ttbar1lep")
            TT0LepSamp = Sample(input_ttbar_0lep,  color=48,
                                legend="t#bar{t}", name="ttbar0lep")

        if not doTest:
            sampMan = SampleManager(DataSamp, [
                                    DYSamp, TTbarSamp, WWSamp, WZSamp, ZZSamp, ZXXSamp, TT1LepSamp, TT0LepSamp])
        else:
            sampMan = SampleManager(DataSamp, [TTbarSamp])
        sampMan.groupMCs(["WW", "WZ", "ZZ", "ZXX"], "EWK", 216, "EWK")
        sampMan.groupMCs(
            ["ttbar2lep", "ttbar1lep", "ttbar0lep"], "ttbar", 96, "t#bar{t}")
    else:
        #
        # 5TeV dataset
        #
        DYSamp = Sample(input_dy,    color=92,
                        legend=lepchannel, name="DY",    is5TeV=True, reweightZpt=True)
        TTbarSamp = Sample(input_ttbar, color=96,
                           legend="t#bar{t}", name="ttbar", is5TeV=True)
        WWSamp = Sample(input_ww,    color=38, legend="WW2L",
                        name="WW",    is5TeV=True)
        WZSamp = Sample(input_wz,    color=39, legend="WZ3L",
                        name="WZ",    is5TeV=True)
        ZZ2LSamp = Sample(input_zz2l,  color=37, legend="ZZ2L",
                          name="ZZ2L",  is5TeV=True)
        ZZ4LSamp = Sample(input_zz4l,  color=37, legend="ZZ4L",
                          name="ZZ4L",  is5TeV=True)
        ZXXSamp = Sample(input_zxx,   color=40, legend="ZXX",
                         name="ZXX",   is5TeV=True, reweightZpt=True)

        sampMan = SampleManager(DataSamp, [
                                DYSamp, TTbarSamp, WWSamp, WZSamp, ZZ2LSamp, ZZ4LSamp, ZXXSamp], is5TeV=True)
        sampMan.groupMCs(["WW", "WZ", "ZZ2L", "ZZ4L", "ZXX"],
                         "EWK", 216, "EWK")

    sampMan.DefineAll("zpt", "Z.Pt()")
    sampMan.DefineAll("zmass", "ZMass")
    sampMan.DefineAll("zy",  "Z.Rapidity()")
    sampMan.DefineAll("zeta",  "Z.Rapidity()")
    sampMan.DefineAll("zphi", "Z.Phi()")

    sampMan.DefineMC("u1_corr", "pU1")
    sampMan.DefineMC("u2_corr", "pU2")
    sampMan.DefineMC("met_corr", "metVars[1]")
    sampMan.DefineMC("met_phi_corr", "metVarsPhi[1]")
    # for data, no corrections applied
    DataSamp.Define("u1_corr", "u1")
    DataSamp.Define("u2_corr", "u2")
    DataSamp.Define("met_corr", "met")
    DataSamp.Define("met_phi_corr", "metPhi")

    sampMan.DefineAll("nPV", "npv")

    sampMan.DefineAll("m1pt",  "lep1.Pt()")
    sampMan.DefineAll("m1eta", "lep1.Eta()")
    sampMan.DefineAll("m2pt",  "lep2.Pt()")
    sampMan.DefineAll("m2eta", "lep2.Eta()")

    sampMan.DefineAll("LeadMuon_pt",     "(m1pt > m2pt ? m1pt : m2pt)")
    sampMan.DefineAll("LeadMuon_eta",    "(m1pt > m2pt ? m1eta : m2eta)")
    sampMan.DefineAll("SubleadLep_pt",  "(m1pt > m2pt ? m2pt : m1pt)")
    sampMan.DefineAll("SubleadLep_eta", "(m1pt > m2pt ? m2eta : m1eta)")

    # DYSamp.Define("u_pt_corr_central", "TMath::Sqrt(u1_corr_central*u1_corr_central + u2_corr_central*u2_corr_central)")
    # sampMan.DefineAll("u1_corr_central",     "u1"  , excludes=['DY'])

    # sample weights
    for i in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]:
        sampMan.DefineMC("weight_{}_tmp".format(str(i)),
                         f"evtWeight[{i}] * self.fnorm * ZptWeight")
        # fix a few very rare nan cases
        sampMan.DefineMC("weight_{}".format(
            str(i)), "TMath::IsNaN(weight_{}_tmp) ? 0.: weight_{}_tmp".format(str(i), str(i)))
        DataSamp.Define("weight_{}".format(str(i)),  "1.0")

    # scale the theory variations back
    # some theory variations are far from 1.0
    # use this to scale them back
    # this should be equavalent as using fnorm_varied
    for i in range(111):
        sampMan.DefineMC(
            f"lheweightS_{i}", f"lheweight[{i}] * self.nmcevt / self.nmcevts_varied[{i}]")

    vsamples = ["DY"]
    sampMan.DefineSpecificMCs("VpT", "genV.Pt()", sampnames=vsamples)
    theoryMaps = theoryUnc_13TeV_zjets if not is5TeV else theoryUnc_5TeV_zjets

    theoryVariations = OrderedDict()
    for vpt in range(len(Vptbins)-1):
        ptmin = Vptbins[vpt]
        ptmax = Vptbins[vpt+1]
        for var in ["MuRVPTUp", "MuRVPTDown", "MuFVPTUp", "MuFVPTDown", "MuFMuRVPTUp", "MuFMuRVPTDown"]:
            idx = theoryMaps.getIndex(var.replace("VPT", ""))
            sampMan.DefineSpecificMCs(f"{var}".replace("VPT", str(
                vpt)), f"(VpT >= {ptmin} && VpT < {ptmax}) ? lheweightS_{idx}: 1.0", sampnames=vsamples)
            sampMan.DefineAll(f"{var}".replace(
                "VPT", str(vpt)), "1.0", excludes=vsamples)

            theoryVariations[f"{var}".replace("VPT", str(
                vpt))] = f"{var}".replace("VPT", str(vpt))
    # pdf
    pdfStart = theoryMaps.getIndex("PDFStart")
    for i in range(1, 101):
        theoryVariations[f"PDF{i}Up"] = f"lheweightS_{i + pdfStart - 1}"
    # alphaS
    alphaSUp = theoryMaps.getIndex("alphaSUp")
    alphaSDown = theoryMaps.getIndex("alphaSDown")
    theoryVariations["alphaSUp"] = f"lheweightS_{alphaSUp}"
    theoryVariations["alphaSDown"] = f"lheweightS_{alphaSDown}"

    for var in theoryVariations.values():
        sampMan.DefineMC(
            f"weight_{var}", f"(TMath::IsNaN(evtWeight[0])) ? 0. : evtWeight[0] * ZptWeight * self.fnorm * {var}")
        DataSamp.Define(f"weight_{var}",  "1.0")

    sampMan.DefineAll(
        "weight_noEcal", "prefireEcal == 0 ? 1 : (weight_0 / prefireEcal)")
    sampMan.DefineAll(
        "weight_noMuon", "prefireMuon == 0 ? 1 : (weight_0 / prefireMuon)")

    if doMuon:
        lepname = "mumu"
        leplabel = "#mu"
    else:
        lepname = "ee"
        leplabel = "e"

    met_pt_bins = np.array([0., 2.0, 4., 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 33,
                           36, 39, 42, 45, 48, 51, 55, 60, 65, 70, 75, 80, 90, 100, 110, 120, 135, 150, 165, 180, 200])
    u1_bins = np.array([-20.0, -16, -12, -10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22,
                       24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 53, 56, 59, 64, 68, 72, 76, 80, 85, 90, 100])
    u2_bins = np.array([-40, -35, -30, -25., -22., -20, -18, -16, -14, -12, -
                       10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10, 12, 16, 18, 20, 22, 25, 30, 35, 40])
    # u_bins = np.array([0., 2., 4., 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 43, 46, 49, 52, 56, 60, 64, 68, 72, 76, 80, 85, 90, 95, 100, 105, 110, 115, 120, 130, 140, 150])
    phimin = -ROOT.TMath.Pi()
    phimax = ROOT.TMath.Pi()

    eta_bins = np.concatenate(
        (np.linspace(-2.4, -2.0, 3), np.linspace(-1.8, 1.8, 37), np.linspace(2.0, 2.4, 3)))
    pt_bins = np.concatenate((np.linspace(24, 35, 6), np.linspace(
        36, 55, 20), np.linspace(57, 63, 4), np.linspace(66, 70, 2)))
    mass_bins = np.array([60, 64, 68, 72, 74, 76, 78, 80, 82, 84.0, 86.0, 87.0, 88.0,
                         89.0, 90.0, 91., 92., 93., 94., 96., 98., 100., 102., 105., 108., 112., 116., 120.])
    u_bins = np.concatenate((np.linspace(0, 20, 11), np.linspace(
        24, 80, 15), np.linspace(85, 109, 4), np.linspace(120, 150, 3)))

    sampMan.cacheDraw("nPV", "histo_nPV_" + lepname, 10, 0, 10, DrawConfig(
        xmin=0, xmax=10, xlabel='# PV', ylabel='Events / 1', donormalizebin=False))

    # z pt befor and after pt reweighting
    sampMan.cacheDraw("zpt",   "histo_zjets_zpt_WoZptWeight_" + lepname, u_bins, DrawConfig(xmin=0, xmax=150,
                      xlabel='p_{{T}}^{{{leplabel}{leplabel}}} [GeV]'.format(leplabel=leplabel), yrmax=1.19, yrmin=0.81), weightname="weight_WoVpt")
    sampMan.cacheDraw("zpt",   "histo_zjets_zpt_WZpTWeight_" + lepname, u_bins, DrawConfig(xmin=0, xmax=150,
                      xlabel='p_{{T}}^{{{leplabel}{leplabel}}} [GeV]'.format(leplabel=leplabel), yrmax=1.19, yrmin=0.81), weightname="weight_WVpt")
    sampMan.cacheDraw("zmass", "histo_zjets_zmass_" + lepname, mass_bins, DrawConfig(
        xmin=70, xmax=110, xlabel='m_{{{leplabel}{leplabel}}} [GeV]'.format(leplabel=leplabel)))

    # z eta with prefiring variations
    sampMan.cacheDraw("zeta",  "histo_zjets_zeta_" + lepname, eta_bins, DrawConfig(xmin=-2.5, xmax=2.5,
                      xlabel='#eta_{{{leplabel}{leplabel}}}'.format(leplabel=leplabel), ymax=1e7, ylabel='Events / 1', yrmax=1.19, yrmin=0.81))
    # sampMan.cacheDraw("zeta",  "histo_zjets_zetaPrefUp_" + lepname, eta_bins, DrawConfig(xmin=-2.5, xmax=2.5, xlabel='#eta_{{{leplabel}{leplabel}}}'.format(leplabel=leplabel), ymax=1e7, ylabel='Events / 1', yrmax=1.19, yrmin=0.81), weightname = "weight_6")
    # sampMan.cacheDraw("zeta",  "histo_zjets_zetaPrefDown_" + lepname, eta_bins, DrawConfig(xmin=-2.5, xmax=2.5, xlabel='#eta_{{{leplabel}{leplabel}}}'.format(leplabel=leplabel), ymax=1e7, ylabel='Events / 1', yrmax=1.19, yrmin=0.81), weightname = "weight_7")
    # sampMan.cacheDraw("zeta",  "histo_zjets_zetaEcalPrefUp_" + lepname, eta_bins, DrawConfig(xmin=-2.5, xmax=2.5, xlabel='#eta_{{{leplabel}{leplabel}}}'.format(leplabel=leplabel), ymax=1e7, ylabel='Events / 1', yrmax=1.19, yrmin=0.81), weightname = "weight_8")
    # sampMan.cacheDraw("zeta",  "histo_zjets_zetaEcalPrefDown_" + lepname, eta_bins, DrawConfig(xmin=-2.5, xmax=2.5, xlabel='#eta_{{{leplabel}{leplabel}}}'.format(leplabel=leplabel), ymax=1e7, ylabel='Events / 1', yrmax=1.19, yrmin=0.81), weightname = "weight_9")
    # sampMan.cacheDraw("zeta",  "histo_zjets_zetaMuonPrefUp_" + lepname, eta_bins, DrawConfig(xmin=-2.5, xmax=2.5, xlabel='#eta_{{{leplabel}{leplabel}}}'.format(leplabel=leplabel), ymax=1e7, ylabel='Events / 1', yrmax=1.19, yrmin=0.81), weightname = "weight_10")
    # sampMan.cacheDraw("zeta",  "histo_zjets_zetaMuonPrefDown_" + lepname, eta_bins, DrawConfig(xmin=-2.5, xmax=2.5, xlabel='#eta_{{{leplabel}{leplabel}}}'.format(leplabel=leplabel), ymax=1e7, ylabel='Events / 1', yrmax=1.19, yrmin=0.81), weightname = "weight_11")

    # z y with prefiring variations
    sampMan.cacheDraw("zy",    "histo_zjets_zrapidity_" + lepname, eta_bins, DrawConfig(xmin=-2.5, xmax=2.5,
                      xlabel='y_{{{leplabel}{leplabel}}}'.format(leplabel=leplabel), ymax=1e7, ylabel='Events / 1', yrmin=0.81, yrmax=1.19))
    sampMan.cacheDraw("zy",    "histo_zjets_zrapidityPrefUp_" + lepname, eta_bins, DrawConfig(xmin=-2.5, xmax=2.5,
                      xlabel='y_{{{leplabel}{leplabel}}}'.format(leplabel=leplabel), ymax=1e7, ylabel='Events / 1', yrmin=0.81, yrmax=1.19), weightname="weight_6")
    sampMan.cacheDraw("zy",    "histo_zjets_zrapidityPrefDown_" + lepname, eta_bins, DrawConfig(xmin=-2.5, xmax=2.5,
                      xlabel='y_{{{leplabel}{leplabel}}}'.format(leplabel=leplabel), ymax=1e7, ylabel='Events / 1', yrmin=0.81, yrmax=1.19), weightname="weight_7")
    sampMan.cacheDraw("zy",    "histo_zjets_zrapidityEcalPrefUp_" + lepname, eta_bins, DrawConfig(xmin=-2.5, xmax=2.5,
                      xlabel='y_{{{leplabel}{leplabel}}}'.format(leplabel=leplabel), ymax=1e7, ylabel='Events / 1', yrmin=0.81, yrmax=1.19), weightname="weight_8")
    sampMan.cacheDraw("zy",    "histo_zjets_zrapidityEcalPrefDown_" + lepname, eta_bins, DrawConfig(xmin=-2.5, xmax=2.5,
                      xlabel='y_{{{leplabel}{leplabel}}}'.format(leplabel=leplabel), ymax=1e7, ylabel='Events / 1', yrmin=0.81, yrmax=1.19), weightname="weight_9")
    sampMan.cacheDraw("zy",    "histo_zjets_zrapidityMuonPrefUp_" + lepname, eta_bins, DrawConfig(xmin=-2.5, xmax=2.5,
                      xlabel='y_{{{leplabel}{leplabel}}}'.format(leplabel=leplabel), ymax=1e7, ylabel='Events / 1', yrmin=0.81, yrmax=1.19), weightname="weight_10")
    sampMan.cacheDraw("zy",    "histo_zjets_zrapidityMuonPrefDown_" + lepname, eta_bins, DrawConfig(xmin=-2.5, xmax=2.5,
                      xlabel='y_{{{leplabel}{leplabel}}}'.format(leplabel=leplabel), ymax=1e7, ylabel='Events / 1', yrmin=0.81, yrmax=1.19), weightname="weight_11")
    sampMan.cacheDraw("zy",    "histo_zjets_zrapidityWoEcalPrefire_" + lepname, eta_bins, DrawConfig(xmin=-2.5, xmax=2.5,
                      xlabel='y_{{{leplabel}{leplabel}}}'.format(leplabel=leplabel), ymax=1e7, ylabel='Events / 1', yrmin=0.81, yrmax=1.19), weightname="weight_noEcal")
    sampMan.cacheDraw("zy",    "histo_zjets_zrapidityWoMuonPrefire_" + lepname, eta_bins, DrawConfig(xmin=-2.5, xmax=2.5,
                      xlabel='y_{{{leplabel}{leplabel}}}'.format(leplabel=leplabel), ymax=1e7, ylabel='Events / 1', yrmin=0.81, yrmax=1.19), weightname="weight_noMuon")

    sampMan.cacheDraw("LeadMuon_pt", "histo_leadLep_pt_" + lepname, pt_bins, DrawConfig(
        xmin=20, xmax=70, xlabel='p_{{T}}(Leading {leplabel}) [GeV]'.format(leplabel=leplabel)))
    # leading lepton eta with prefiring variations
    sampMan.cacheDraw("LeadMuon_eta", "histo_leadLep_eta_" + lepname, eta_bins, DrawConfig(xmin=-2.6, xmax=2.6,
                      xlabel='#eta (Leading {leplabel})'.format(leplabel=leplabel), ymax=1e7, ylabel='Events / 1', yrmin=0.81, yrmax=1.19))
    sampMan.cacheDraw("LeadMuon_eta", "histo_leadLep_etaPrefUp_" + lepname, eta_bins, DrawConfig(xmin=-2.6, xmax=2.6,
                      xlabel='#eta (Leading {leplabel})'.format(leplabel=leplabel), ymax=1e7, ylabel='Events / 1', yrmin=0.81, yrmax=1.19), weightname="weight_6")
    sampMan.cacheDraw("LeadMuon_eta", "histo_leadLep_etaPrefDown_" + lepname, eta_bins, DrawConfig(xmin=-2.6, xmax=2.6,
                      xlabel='#eta (Leading {leplabel})'.format(leplabel=leplabel), ymax=1e7, ylabel='Events / 1', yrmin=0.81, yrmax=1.19), weightname="weight_7")

    sampMan.cacheDraw("SubleadLep_pt", "histo_subleadLep_pt_" + lepname, pt_bins, DrawConfig(
        xmin=20, xmax=70, xlabel='p_{{T}}(Subleading {leplabel}) [GeV]'.format(leplabel=leplabel)))
    # subleading lepton eta with prefiring variations
    sampMan.cacheDraw("SubleadLep_eta", "histo_subleadLep_eta_" + lepname, eta_bins, DrawConfig(xmin=-2.6, xmax=2.6,
                      xlabel='#eta (Subleading {leplabel})'.format(leplabel=leplabel), ymax=1e7, ylabel='Events / 1', yrmin=0.81, yrmax=1.19))
    sampMan.cacheDraw("SubleadLep_eta", "histo_subleadLep_etaPrefUp_" + lepname, eta_bins, DrawConfig(xmin=-2.6, xmax=2.6,
                      xlabel='#eta (Subleading {leplabel})'.format(leplabel=leplabel), ymax=1e7, ylabel='Events / 1', yrmin=0.81, yrmax=1.19), weightname="weight_6")
    sampMan.cacheDraw("SubleadLep_eta", "histo_subleadLep_etaPrefDown_" + lepname, eta_bins, DrawConfig(xmin=-2.6, xmax=2.6,
                      xlabel='#eta (Subleading {leplabel})'.format(leplabel=leplabel), ymax=1e7, ylabel='Events / 1', yrmin=0.81, yrmax=1.19), weightname="weight_7")

    sampMan.cacheDraw("met", "histo_zjets_pfmet_pt" + lepname,
                      met_pt_bins, DrawConfig(xmin=0, xmax=100, xlabel='PF MET [GeV]'))
    sampMan.cacheDraw("metPhi", "histo_zjets_pfmet_phi" + lepname, 30, phimin,
                      phimax, DrawConfig(xmin=phimin, xmax=phimax, xlabel='PF MET #phi'))

    sampMan.cacheDraw("met_corr", "histo_zjets_pfmet_pt_corr_" + lepname,
                      met_pt_bins, DrawConfig(xmin=0, xmax=100, xlabel='Corrected PF MET [GeV]'))
    sampMan.cacheDraw("met_phi_corr", "histo_zjets_pfmet_phi_corr_" + lepname, 30, phimin,
                      phimax, DrawConfig(xmin=phimin, xmax=phimax, xlabel='Corrected PF MET #phi'))

    sampMan.cacheDraw("u1_corr", "histo_zjets_pfmet_u1_corr_" + lepname, u1_bins,
                      DrawConfig(xmin=-20.0, xmax=50.0, xlabel='Corrected u_{#parallel} [GeV]'))
    sampMan.cacheDraw("u2_corr", "histo_zjets_pfmet_u2_corr_" + lepname, u2_bins,
                      DrawConfig(xmin=-40.0, xmax=40.0, xlabel='Corrected u_{#perp } [GeV]'))

    sampMan.cacheDraw("u1", "histo_zjets_pfmet_u1_" + lepname, u1_bins,
                      DrawConfig(xmin=-20.0, xmax=50.0, xlabel='u_{#parallel} [GeV]'))
    sampMan.cacheDraw("u2", "histo_zjets_pfmet_u2_" + lepname, u2_bins,
                      DrawConfig(xmin=-40.0, xmax=40.0, xlabel='u_{#perp } [GeV]'))

    # sampMan.cacheDraw("ucor1", "histo_zjets_mumu_pfmet_ucor1_pt", u1_bins, DrawConfig(xmin=-20.0, xmax=100.0, xlabel="u_{#parallel} [GeV]", extraText='#mu#mu channel'))
    # sampMan.cacheDraw("ucor2", "histo_zjets_mumu_pfmet_ucor2_pt", u2_bins, DrawConfig(xmin=-30.0, xmax=30.0, xlabel= "u_{#perp  } [GeV]", extraText='#mu#mu channel'))

    sampMan.cacheDraw("zmass", "histo_zjets_zmass_" + lepname, mass_bins, DrawConfig(
        xmin=60, xmax=120, xlabel='m_{{{leplabel}{leplabel}}} [GeV]'.format(leplabel=leplabel)))
    for i in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]:
        sampMan.cacheDraw("zmass", "histo_zjets_zmass_{}_weight_{}".format(lepname, str(i)),  15, 60, 120, DrawConfig(xmin=60, xmax=120, xlabel="m_{{{leplabel}{leplabel}}} [GeV]".format(
            leplabel=leplabel), dology=True, donormalizebin=False, addOverflow=False, addUnderflow=False), weightname="weight_{}".format(str(i)))
    for var in theoryVariations.values():
        sampMan.cacheDraw("zmass", f"histo_zjets_zmass_{lepname}_{var}", 15, 60, 120, DrawConfig(
            xmin=60, xmax=120, xlabel=f"m_{{{leplabel}{leplabel}}} [GeV]", dology=True, donormalizebin=False, addOverflow=False, addUnderflow=False), weightname=f"weight_{var}")

    sampMan.launchDraw()

    # Draw the prefire correction histogram with uncertainty band
    hname = "histo_zjets_zrapidity_" + lepname
    hnameUp = "histo_zjets_zrapidityPrefUp_" + lepname
    hnameDown = "histo_zjets_zrapidityPrefDown_" + lepname
    hratiopanel = sampMan.hratios[hname][0]
    hratiopanelUp = sampMan.hratios[hnameUp][0]
    hratiopanelDown = sampMan.hratios[hnameDown][0]
    print((hratiopanel.GetNbinsX()))
    for ibin in range(hratiopanel.GetNbinsX() + 1):
        # print(("ibin {} bin content {} up {} down {}".format(ibin, hratiopanel.GetBinContent(ibin), hratiopanelUp.GetBinContent(ibin), hratiopanelDown.GetBinContent(ibin))))
        hratiopanel.SetBinContent(ibin, 1.0)
        hratiopanel.SetBinError(
            ibin, 0.5*abs(hratiopanelUp.GetBinContent(ibin) - hratiopanelDown.GetBinContent(ibin)))
    drawconfigs = DrawConfig(xmin=-2.5, xmax=2.5, xlabel='y_{{{leplabel}{leplabel}}}'.format(
        leplabel=leplabel), ymax=1e7, ylabel='Events / 1', yrmin=0.81, yrmax=1.19, outputname=hname + "_withUnc")
    DrawHistos([sampMan.hdatas[hname], sampMan.hsmcs[hname]], ["Data", lepchannel, "t#bar{t}", "EWK"], drawconfigs.xmin, drawconfigs.xmax, drawconfigs.xlabel, drawconfigs.ymin, drawconfigs.ymax, drawconfigs.ylabel, drawconfigs.outputname, is5TeV=is5TeV, hratiopanel=hratiopanel, dology=drawconfigs.dology, dologx=drawconfigs.dologx, showratio=drawconfigs.showratio, yrmax=drawconfigs.yrmax,
               yrmin=drawconfigs.yrmin, yrlabel=drawconfigs.yrlabel, donormalize=drawconfigs.donormalize, ratiobase=drawconfigs.ratiobase, legendPos=drawconfigs.legendPos, redrawihist=drawconfigs.redrawihist, extraText=drawconfigs.extraText, noCMS=drawconfigs.noCMS, addOverflow=drawconfigs.addOverflow, addUnderflow=drawconfigs.addUnderflow, nMaxDigits=drawconfigs.nMaxDigits)

    # htests = sampMan.hsmcs["histo_zjets_zmass_{}_weight_8".format(lepname)]
    # htest = list(htests.GetHists())
    # for h in htest:
    #    print("integral: ", h.Integral())
    #    print("bin contents: ", h.Print("all"))

    sqrtS = "5TeV" if is5TeV else "13TeV"

    # zpt ratio
    if False:
        output_suffix = f"{lepname}_{sqrtS}.root"
        outfile = ROOT.TFile.Open(
            f"data/output_Zpt_{output_suffix}", "recreate")
        hdata_zpt = sampMan.hdatas[f"histo_zjets_zpt_WoZptWeight_{lepname}"]
        hsmcs_zpt = sampMan.hsmcs[f"histo_zjets_zpt_WoZptWeight_{lepname}"]
        hratio_zpt = list(
            sampMan.hratios[f"histo_zjets_zpt_WoZptWeight_{lepname}"].values())[0]

        hmc_zpt = THStack2TH1(hsmcs_zpt)

        hdata_zpt.SetName(f"histo_zjets_zpt_{lepname}_Data")
        hmc_zpt.SetName(f"histo_zjets_zpt_{lepname}_MC")
        hratio_zpt.SetName(f"histo_zjets_zpt_{lepname}_Ratio")

        hdata_zpt.SetDirectory(outfile)
        hdata_zpt.Write()
        hmc_zpt.SetDirectory(outfile)
        hmc_zpt.Write()
        hratio_zpt.SetDirectory(outfile)
        hratio_zpt.Write()

        outfile.Close()

    output_suffix = f"{lepname}_{sqrtS}_noTheoryNorm.root"
    if doTheoryNorm:
        output_suffix = f"{lepname}_{sqrtS}_TheoryNormed.root"
    outfile = ROOT.TFile.Open("root/output_shapes_"+output_suffix, "recreate")

    # Data
    hdata = sampMan.hdatas["histo_zjets_zmass_{}_weight_0".format(lepname)]
    hdata.SetDirectory(outfile)
    hdata.Write()

    hsmcs_central = sampMan.hsmcs["histo_zjets_zmass_{}_weight_0".format(
        lepname)]
    hlists_central = list(hsmcs_central.GetHists())

    # central values for MC
    for ih in range(len(hlists_central)):
        hcen = hlists_central[ih]
        hcen.SetDirectory(outfile)
        hcen.Write()

    if lepname == "mumu":
        channel = "mu"
    else:
        channel = "e"
    # weights/corrections
    for i in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]:
        hsmcs_up = sampMan.hsmcs["histo_zjets_zmass_{}_weight_{}".format(
            lepname, str(i))]
        hlists_up = list(hsmcs_up.GetHists())
        for ih in range(len(hlists_up)):
            # loop over different processes
            hcen = hlists_central[ih]
            hup = hlists_up[ih]

            prefix = ""
            if i == 2 or i == 3 or i == 4:
                # the systematics that should be decorrelated between lepton flavors
                prefix = channel + "_"

            if i <= 6 or i == 8 or i == 10:
                suffix = prefix + "SysWeight{}_{}Up".format(str(i), sqrtS)
            else:
                suffix = prefix + "SysWeight{}_{}Down".format(str(i-1), sqrtS)

            # if i==5:
            #    # need to uncorrelated statistical uncertainties (weight5)
            #    hup.SetName("{}_{}_SysWeight{}Up".format(hcen.GetName(), lepname, str(i)))
            # elif i<=6 or i==8 or i==10:
            #    hup.SetName("{}_{}_SysWeight{}Up".format(hcen.GetName(), lepname, str(i)))
            # else:
            #    hup.SetName("{}_{}_SysWeight{}Down".format(hcen.GetName(), lepname, str(i-1)))
            hup.SetName("{}_{}".format(hcen.GetName(), suffix))

            # if i==5:
            #    # need to uncorrelated statistical uncertainties (weight5)
            #    hdn  = hcen.Clone("{}_{}_SysWeight{}Down".format(hcen.GetName(), lepname, str(i)))
            # elif i<6:
            #    # for prefire, the up and down weights are calculated seperately
            #    hdn  = hcen.Clone("{}_{}_SysWeight{}Down".format(hcen.GetName(), lepname, str(i)))
            if i < 6:
                suffix = prefix + "SysWeight{}_{}Down".format(str(i), sqrtS)
                hdn = hcen.Clone("{}_{}".format(hcen.GetName(), suffix))
                for ibin in range(1, hup.GetNbinsX()+1):
                    hdn.SetBinContent(
                        ibin, 2*hcen.GetBinContent(ibin) - hup.GetBinContent(ibin))

            hup.SetDirectory(outfile)
            hup.Write()
            if i < 6:
                hdn.SetDirectory(outfile)
                hdn.Write()

    # theory uncertainties
    for sysname, var in theoryVariations.items():
        hsmcs_up = sampMan.hsmcs[f"histo_zjets_zmass_{lepname}_{var}"]
        hlists_up = list(hsmcs_up.GetHists())
        for ih in range(len(hlists_up)):
            # loop over different processes/samps
            hcen = hlists_central[ih]
            hup = hlists_up[ih]

            suffix = sqrtS + "_" + sysname
            hup.SetName("{}_{}".format(hcen.GetName(), suffix))
            for ibin in range(1, hup.GetNbinsX()+1):
                hup.SetBinError(ibin, 0.)

            if doTheoryNorm:
                # normalize the number of events to the central values
                # if the theory variations are normalized
                # only apply to signal process
                if "DY" in hcen.GetName():
                    hup.Scale(hcen.Integral()/hup.Integral())

            if "PDF" in sysname:
                suffix = sqrtS + "_" + sysname.replace("Up", "Down")
                hdn.SetName("{}_{}".format(hcen.GetName(), suffix))
                for ibin in range(1, hup.GetNbinsX()+1):
                    hdn.SetBinContent(
                        ibin, 2*hcen.GetBinContent(ibin) - hup.GetBinContent(ibin))
                    hdn.SetBinError(ibin, 0.)
                hdn.SetDirectory(outfile)
                hdn.Write()

            hup.SetDirectory(outfile)
            hup.Write()

    outfile.Close()

    sampMan.dumpCounts()
    print("Program end...")

    return


if __name__ == "__main__":
    main()

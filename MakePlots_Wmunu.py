import ROOT
import numpy as np
from collections import OrderedDict
import sys
sys.path.append("/uscms_data/d3/yfeng/CMSPLOTS")
from tdrstyle import setTDRStyle
from myFunction import DrawHistos, THStack2TH1
import CMS_lumi
import pickle

from SampleManager import DrawConfig, Sample, SampleManager

ROOT.gROOT.SetBatch(True)
ROOT.ROOT.EnableImplicitMT(15)

dotest = 0
VetoB = False
doSys = False

def main():
    print "Program start..."

    #ROOT.gROOT.ProcessLine('TFile* f_zpt = TFile::Open("results/zpt_weight.root")')
    #ROOT.gROOT.ProcessLine('TH1D* h_zpt_ratio  = (TH1D*)f_zpt->Get("h_zpt_ratio")')

    input_data    = "inputs/wmunu/input_data.txt"
    input_wm0     = "inputs/wmunu/input_wm0.txt"
    input_wm1     = "inputs/wmunu/input_wm1.txt"
    input_wm2     = "inputs/wmunu/input_wm2.txt"
    input_dy      = "inputs/wmunu/input_zjets.txt"
    input_ttbar   = "inputs/wmunu/input_ttbar_dilepton.txt"
    input_ttbar_1lep = "inputs/wmunu/input_ttbar_singlelepton.txt"
    input_ttbar_0lep = "inputs/wmunu/input_ttbar_hadronic.txt"
    input_ww      = "inputs/wmunu/input_ww.txt"
    input_wz      = "inputs/wmunu/input_wz.txt"
    input_zz      = "inputs/wmunu/input_zz.txt"
    input_zxx     = "inputs/wmunu/input_zxx.txt"
    input_wx0     = "inputs/wmunu/input_wx0.txt"
    input_wx1     = "inputs/wmunu/input_wx1.txt"
    input_wx2     = "inputs/wmunu/input_wx2.txt"

    input_antiiso_data    = "inputs/awmunu/input_data.txt"
    input_antiiso_wm0     = "inputs/awmunu/input_wm0.txt"
    input_antiiso_wm1     = "inputs/awmunu/input_wm1.txt"
    input_antiiso_wm2     = "inputs/awmunu/input_wm2.txt"
    input_antiiso_ttbar   = "inputs/awmunu/input_ttbar_dilepton.txt"
    input_antiiso_ttbar_1lep = "inputs/awmunu/input_ttbar_singlelepton.txt"
    input_antiiso_ttbar_0lep = "inputs/awmunu/input_ttbar_hadronic.txt"
    input_antiiso_ww      = "inputs/awmunu/input_ww.txt"
    input_antiiso_wz      = "inputs/awmunu/input_wz.txt"
    input_antiiso_zz      = "inputs/awmunu/input_zz.txt"
    input_antiiso_zxx     = "inputs/awmunu/input_zxx.txt"
    input_antiiso_wx0     = "inputs/awmunu/input_wx0.txt"
    input_antiiso_wx1     = "inputs/awmunu/input_wx1.txt"
    input_antiiso_wx2     = "inputs/awmunu/input_wx2.txt"

    DataSamp  = Sample(input_data, isMC=False, legend="Data", name="Data", isWSR=True)
    # W -> munu
    wjetsnorm = 1.0
    Wm0Samp   = Sample(input_wm0, xsec = 49155*1e3, isMC=True, name = "wm0", isWSR=True, additionalnorm = wjetsnorm, isNLO=True)
    Wm1Samp   = Sample(input_wm1, xsec = 8047*1e3,  isMC=True, name = "wm1", isWSR=True, additionalnorm = wjetsnorm, isNLO=True)
    Wm2Samp   = Sample(input_wm2, xsec = 3160*1e3,  isMC=True, name = "wm2", isWSR=True, additionalnorm = wjetsnorm, isNLO=True)
    # DY
    DYSamp    = Sample(input_dy,    xsec = 6077.22*1e3,      isMC=True, name="DY", isWSR=True, isNLO=True)
    # ttbar
    TTbarSamp  = Sample(input_ttbar,      xsec = 88.29*1e3,  isMC=True, name = "ttbar_dilepton", isWSR=True)
    TT1LepSamp = Sample(input_ttbar_1lep, xsec = 365.35*1e3, isMC=True, name = "ttbar_1lepton",  isWSR=True)
    TT0LepSamp = Sample(input_ttbar_0lep, xsec = 377.96*1e3, isMC=True, name = "ttbar_0lepton",  isWSR=True)
    ## dibosons
    WWSamp = Sample(input_ww, xsec = 12.6*1e3,  isMC=True, name = "WW", isWSR=True)
    WZSamp = Sample(input_wz, xsec = 4.912*1e3, isMC=True, name = "WZ", isWSR=True)
    ZZSamp = Sample(input_zz, xsec = 16.52*1e3, isMC=True, name = "ZZ", isWSR=True)
    ## tau
    #ZXXSamp = Sample(input_zxx, isMC=True, name = "ZXX", isWSR=True)
    #Wx0Samp = Sample(input_wx0, isMC=True, name = "wx0", isWSR=True)
    #Wx1Samp = Sample(input_wx1, isMC=True, name = "wx1", isWSR=True)
    #Wx2Samp = Sample(input_wx2, isMC=True, name = "wx2", isWSR=True)

    ### for the QCD background estimation (data-driven)
    #qcdnorm = 1.0243
    #DataAisoSamp  = Sample(input_antiiso_data, isMC=False, name="Data_aiso", isWSR=True, additionalnorm = qcdnorm, legend = 'QCD', color='226')
    ## W -> munu
    #Wm0AisoSamp   = Sample(input_antiiso_wm0, isMC=True, name = "wm0_aiso", isWSR=True, additionalnorm=-1.0 * qcdnorm)
    #Wm1AisoSamp   = Sample(input_antiiso_wm1, isMC=True, name = "wm1_aiso", isWSR=True, additionalnorm=-1.0 * qcdnorm)
    #Wm2AisoSamp   = Sample(input_antiiso_wm2, isMC=True, name = "wm2_aiso", isWSR=True, additionalnorm=-1.0 * qcdnorm)
    ## ttbar
    #TTbarAisoSamp  = Sample(input_antiiso_ttbar, isMC=True, name = "ttbar_dilepton_aiso", isWSR=True, additionalnorm=-1.0 * qcdnorm)
    #TT1LepAisoSamp = Sample(input_antiiso_ttbar_1lep, isMC=True, name = "ttbar_1lepton_aiso", isWSR=True, additionalnorm=-1.0 * qcdnorm)
    #TT0LepAisoSamp = Sample(input_antiiso_ttbar_0lep, isMC=True, name = "ttbar_0lepton_aiso", isWSR=True, additionalnorm=-1.0 * qcdnorm)
    ### dibosons
    #WWAisoSamp = Sample(input_antiiso_ww, isMC=True, name = "WW_aiso", isWSR=True, additionalnorm=-1.0 * qcdnorm)
    #WZAisoSamp = Sample(input_antiiso_wz, isMC=True, name = "WZ_aiso", isWSR=True, additionalnorm=-1.0 * qcdnorm)
    #ZZAisoSamp = Sample(input_antiiso_zz, isMC=True, name = "ZZ_aiso", isWSR=True, additionalnorm=-1.0 * qcdnorm)
    ## tau
    #ZXXAisoSamp = Sample(input_antiiso_zxx, isMC=True, name = "ZXX_aiso", isWSR=True, additionalnorm=-1.0 * qcdnorm)
    #Wx0AisoSamp = Sample(input_antiiso_wx0, isMC=True, name = "wx0_aiso", isWSR=True, additionalnorm=-1.0 * qcdnorm)
    #Wx1AisoSamp = Sample(input_antiiso_wx1, isMC=True, name = "wx1_aiso", isWSR=True, additionalnorm=-1.0 * qcdnorm)
    #Wx2AisoSamp = Sample(input_antiiso_wx2, isMC=True, name = "wx2_aiso", isWSR=True, additionalnorm=-1.0 * qcdnorm)

    #sampMan = SampleManager(DataSamp, [Wm0Samp, Wm1Samp, Wm2Samp, TTbarSamp, TT1LepSamp, TT0LepSamp, WWSamp, WZSamp, ZZSamp, ZXXSamp, Wx0Samp, Wx1Samp, Wx2Samp, 
    #        DataAisoSamp, Wm0AisoSamp, Wm1AisoSamp, Wm2AisoSamp, TTbarAisoSamp, TT1LepAisoSamp, TT0LepAisoSamp, WWAisoSamp, WZAisoSamp, ZZAisoSamp, ZXXAisoSamp, Wx0AisoSamp, Wx1AisoSamp, Wx2AisoSamp])
    sampMan = SampleManager(DataSamp, [Wm0Samp, Wm1Samp, Wm2Samp, DYSamp, TTbarSamp, TT1LepSamp, TT0LepSamp, WWSamp, WZSamp, ZZSamp])
    #sampMan.groupMCs(["WW", "WZ", "ZZ", "ZXX", "wx0", "wx1", "wx2"], "EWK", 216, "EWK")
    sampMan.groupMCs(["DY", "WW", "WZ", "ZZ"], "EWK", 216, "EWK")
    sampMan.groupMCs(["ttbar_dilepton", "ttbar_1lepton", "ttbar_0lepton"], "ttbar", 96, "t#bar{t}")
    sampMan.groupMCs(['wm0', 'wm1', 'wm2'], "wmunu", 92,"W#rightarrow#mu#nu")
    #sampMan.groupMCs(['Data_aiso', 'wm0_aiso', 'wm1_aiso', 'wm2_aiso', 'ttbar_dilepton_aiso', 'ttbar_1lepton_aiso', 'ttbar_0lepton_aiso', 'WW_aiso', 'WZ_aiso', 'ZZ_aiso', 'ZXX_aiso', 'wx0_aiso', 'wx1_aiso' ,'wx2_aiso'], 'QCD', 226, 'QCD')

    sampMan.ApplyCutAll("passSelection > 0")
    sampMan.DefineAll("muon_pt",  "Muon_corrected_pt[0]")
    sampMan.DefineAll("muon_eta", "Muon_eta[0]")
    sampMan.DefineAll("met_pt",   "RawMET_corrected_pt")
    sampMan.DefineAll("met_phi",  "RawMET_corrected_phi")

    sampMan.DefineAll("weight_wplus", "weight_WoVpt * (Muon_charge[0]>0)")
    sampMan.DefineAll("weight_wminus", "weight_WoVpt * (Muon_charge[0]<0)")

    # deltaPhi
    sampMan.DefineAll("deltaPhi", "TVector2::Phi_0_2pi(Muon_phi[0] - met_phi)")
    sampMan.DefineAll("mT", "TMath::Sqrt(2 *muon_pt * met_pt * (1 - TMath::Cos(deltaPhi)))")

    # muon isolation selection
    #sampMan.ApplyCutAll("Muon_pfRelIso04_all[0] < 0.15")
    sampMan.DefineAll("passIso",  "Muon_pfRelIso04_all[0] < 0.15")
    sampMan.DefineAll("IsoCR1",   "Muon_pfRelIso04_all[0] > 0.20 && Muon_pfRelIso04_all[0] < 0.25")
    sampMan.DefineAll("IsoCR2",   "Muon_pfRelIso04_all[0] > 0.25 && Muon_pfRelIso04_all[0] < 0.30")
    sampMan.DefineAll("IsoCR3",   "Muon_pfRelIso04_all[0] > 0.30 && Muon_pfRelIso04_all[0] < 0.35")
    sampMan.DefineAll("IsoCR4",   "Muon_pfRelIso04_all[0] > 0.35 && Muon_pfRelIso04_all[0] < 0.40")
    # Aram's isolation region
    sampMan.DefineAll("IsoCR5",   "Muon_pfRelIso04_all[0] > 0.30 && Muon_pfRelIso04_all[0] < 0.45")
    sampMan.DefineAll("IsoCR6",   "Muon_pfRelIso04_all[0] > 0.45 && Muon_pfRelIso04_all[0] < 0.55")
    for iso in ["passIso", "IsoCR1", "IsoCR2", "IsoCR3", "IsoCR4", "IsoCR5", "IsoCR6"]:
        sampMan.DefineAll("weight_" + iso,        "weight_WoVpt * {}".format(iso))
        sampMan.DefineAll("weight_wplus_" + iso,  "weight_wplus * {}".format(iso))
        sampMan.DefineAll("weight_wminus_" + iso, "weight_wminus * {}".format(iso))
    sampMan.SetDefaultWeightStr("weight_passIso")
    # mT selection
    #sampMan.ApplyCutAll("mT>40.0")

    # recoil variables
    sampMan.DefineAll("V2W", "UVec(muon_pt, Muon_phi[0], met_pt, met_phi)")
    sampMan.DefineAll("WpT", "V2W.Mod()")
    sampMan.DefineAll("Wphi", "TVector2::Phi_mpi_pi(V2W.Phi())")

    # pTOvermT
    sampMan.DefineAll("pTOvermT", "muon_pt / (mT+1e-6)")

    met_pt_bins = np.array([0., 2.0, 4., 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 33, 36, 39, 42, 45, 48, 51, 55, 60, 65, 70, 75, 80])
    u1_bins = np.array([-40.,-36.,-32., -28., -25., -22.0, -20.0, -18, -16, -14, -12, -10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 53, 56, 59, 64, 68, 72, 76, 80, 85, 90, 100])
    u2_bins = np.array([-80., -70., -65., -60., -56., -52, -48, -44, -40, -37, -34, -31, -28, -25., -22., -20, -18, -16, -14, -12, -10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10, 12, 16, 18, 20, 22, 25, 28, 31, 34, 37, 40, 44, 48, 52, 56, 60, 65, 70, 80])
    #u_bins = np.array([0., 2., 4., 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 43, 46, 49, 52, 56, 60, 64, 68, 72, 76, 80, 85, 90, 95, 100, 105, 110, 115, 120, 130, 140, 150])
    mt_bins = np.concatenate((np.linspace(0, 15, 4), np.linspace(18, 36, 10), np.linspace(37,39, 3), np.linspace(40, 90, 51), np.linspace(92, 100, 5), np.linspace(104, 120, 5)))
    phimin = -ROOT.TMath.Pi()
    phimax = ROOT.TMath.Pi()
    dphi_bins = np.concatenate((np.linspace(0,1, 3), np.linspace(1.2, 2.4, 7),np.linspace(2.5, 3.5, 11), np.linspace(3.6, 4.8, 7), np.linspace(5.0, 6.0, 3))) * ROOT.TMath.Pi() / 3

    eta_bins = np.concatenate((np.linspace(-2.4, -2.0, 3),np.linspace(-1.8,1.8,37), np.linspace(2.0, 2.4, 3)))
    pt_bins = np.concatenate((np.linspace(25, 35, 6),np.linspace(36,55,20), np.linspace(57, 63, 4), np.linspace(66,70,2)))
    mass_bins = np.concatenate((np.linspace(70, 76, 3),np.linspace(78,86,5), np.linspace(87, 96, 10), np.linspace(98,104,5), np.linspace(106, 110, 2)))
    u_bins = np.concatenate((np.linspace(0, 20, 11),np.linspace(24,80,15), np.linspace(85, 109, 4), np.linspace(120,150,3)))
    u_bins = np.concatenate((np.linspace(0, 20, 11),np.linspace(24,60,10)))
    ptOvermT_bins = np.concatenate((np.linspace(0,0.2,3), np.linspace(0.3,1.5,25), np.linspace(1.6,2.0, 5), np.linspace(2.2,3.0, 5)))

    # nPV
    #sampMan.cacheDraw("npv", "histo_nPV_munu", 10, 0, 10, DrawConfig(xmin=0, xmax=10, xlabel='# PV', ylabel='Events / 1', dology=False, ymax = 9e5))
    #sampMan.cacheDraw("npv", "histo_nPV_muplusnu", 10, 0, 10, DrawConfig(xmin=0, xmax=10, xlabel='# PV', ylabel='Events / 1', dology=False, ymax = 5e5), weightname = "weight_wplus")
    #sampMan.cacheDraw("npv", "histo_nPV_muminusnu", 10, 0, 10, DrawConfig(xmin=0, xmax=10, xlabel='# PV', ylabel='Events / 1', dology=False, ymax = 4.5e5), weightname = "weight_wminus")
    # lepton pt
    sampMan.cacheDraw("muon_pt", "histo_muon_pt_munu", pt_bins, DrawConfig(xmin=20, xmax=70, xlabel='p_{T}(#mu) [GeV]', dology=False, ymax=1.05e5))
    sampMan.cacheDraw("muon_eta", "histo_muon_eta_munu", eta_bins, DrawConfig(xmin=-2.6, xmax=2.6, xlabel='#eta (#mu) [GeV]', ylabel='Events / 1', dology=False, ymax=6e5))
    sampMan.cacheDraw("muon_pt", "histo_muon_pt_muplusnu", pt_bins, DrawConfig(xmin=20, xmax=70, xlabel='p_{T}(#mu) [GeV]', dology=False, ymax=6.05e4), weightname = "weight_wplus_passIso")
    sampMan.cacheDraw("muon_eta", "histo_muon_eta_muplusnu", eta_bins, DrawConfig(xmin=-2.6, xmax=2.6, xlabel='#eta (#mu) [GeV]', ylabel='Events / 1', dology=False, ymax=4e5), weightname = "weight_wplus_passIso")
    sampMan.cacheDraw("muon_pt", "histo_muon_pt_muminusnu", pt_bins, DrawConfig(xmin=20, xmax=70, xlabel='p_{T}(#mu) [GeV]', dology=False, ymax=5.05e4), weightname = "weight_wminus_passIso")
    sampMan.cacheDraw("muon_eta", "histo_muon_eta_muminusnu", eta_bins, DrawConfig(xmin=-2.6, xmax=2.6, xlabel='#eta (#mu) [GeV]', ylabel='Events / 1', dology=False, ymax=3.5e5), weightname = "weight_wminus_passIso")

    # met and raw_met
    sampMan.cacheDraw("met_pt", "histo_wjets_munu_pfmet_pt", met_pt_bins, DrawConfig(xmin=0, xmax=80, xlabel='PF MET [GeV]', dology=False, ymax=8e4))
    sampMan.cacheDraw("met_phi", "histo_wjets_munu_pfmet_phi", 30, phimin, phimax, DrawConfig(xmin=phimin, xmax=phimax, xlabel='PF MET #phi', dology=False, ymax=6e5))
    sampMan.cacheDraw("RawMET_pt", "histo_wjets_munu_pfMET_pt", met_pt_bins, DrawConfig(xmin=0, xmax=80, xlabel='Raw PF MET [GeV]', dology=False, ymax=8e4))
    sampMan.cacheDraw("RawMET_phi", "histo_wjets_munu_pfMET_phi", 30, phimin, phimax, DrawConfig(xmin=phimin, xmax=phimax, xlabel='Raw PF MET #phi', dology=False, ymax=6e5))

    #sampMan.cacheDraw("met_pt", "histo_wjets_muplusnu_pfmet_pt", met_pt_bins, DrawConfig(xmin=0, xmax=80, xlabel='PF MET [GeV]', dology=False, ymax=5e4), weightname = "weight_wplus_passIso")
    #sampMan.cacheDraw("met_phi", "histo_wjets_muplusnu_pfmet_phi", 30, phimin, phimax, DrawConfig(xmin=phimin, xmax=phimax, xlabel='PF MET #phi', dology=False, ymax=4e5), weightname = "weight_wplus_passIso")
    #sampMan.cacheDraw("RawMET_pt", "histo_wjets_muplusnu_pfMET_pt", met_pt_bins, DrawConfig(xmin=0, xmax=80, xlabel='Raw PF MET [GeV]', dology=False, ymax=5e4), weightname = "weight_wplus_passIso")
    #sampMan.cacheDraw("RawMET_phi", "histo_wjets_muplusnu_pfMET_phi", 30, phimin, phimax, DrawConfig(xmin=phimin, xmax=phimax, xlabel='Raw PF MET #phi', dology=False, ymax=3.5e5), weightname = "weight_wplus_passIso")

    #sampMan.cacheDraw("met_pt", "histo_wjets_muminusnu_pfmet_pt", met_pt_bins, DrawConfig(xmin=0, xmax=80, xlabel='PF MET [GeV]', dology=False, ymax=4e4), weightname = "weight_wminus_passIso")
    #sampMan.cacheDraw("met_phi", "histo_wjets_muminusnu_pfmet_phi", 30, phimin, phimax, DrawConfig(xmin=phimin, xmax=phimax, xlabel='PF MET #phi', dology=False, ymax=3e5), weightname = "weight_wminus_passIso")
    #sampMan.cacheDraw("RawMET_pt", "histo_wjets_muminusnu_pfMET_pt", met_pt_bins, DrawConfig(xmin=0, xmax=80, xlabel='Raw PF MET [GeV]', dology=False, ymax=4e4), weightname = "weight_wminus_passIso")
    #sampMan.cacheDraw("RawMET_phi", "histo_wjets_muminusnu_pfMET_phi", 30, phimin, phimax, DrawConfig(xmin=phimin, xmax=phimax, xlabel='Raw PF MET #phi', dology=False, ymax=3e5), weightname = "weight_wminus_passIso")

    # deltaPhi
    sampMan.cacheDraw("deltaPhi", "histo_wjets_munu_deltaphi",      dphi_bins,  DrawConfig(xmin=0., xmax=2*phimax, xlabel='#Delta #phi(#mu, MET)', dology=False, ymax=2e6))
    #sampMan.cacheDraw("deltaPhi", "histo_wjets_muplusnu_deltaphi",  dphi_bins,  DrawConfig(xmin=0., xmax=2*phimax, xlabel='#Delta #phi(#mu, MET)', dology=False, ymax=1.4e6), weightname = "weight_wplus_passIso")
    #sampMan.cacheDraw("deltaPhi", "histo_wjets_muminusnu_deltaphi", dphi_bins,  DrawConfig(xmin=0., xmax=2*phimax, xlabel='#Delta #phi(#mu, MET)', dology=False, ymax=1.4e6), weightname = "weight_wminus_passIso")

    # mt
    sampMan.cacheDraw("mT", "histo_wjets_munu_mtcorr_passIso",      mt_bins, DrawConfig(xmin=0, xmax=120, xlabel="m_{T} [GeV]", dology=False, ymax=7e4, donormalizebin=False))
    sampMan.cacheDraw("mT", "histo_wjets_muplusnu_mtcorr_passIso",  mt_bins, DrawConfig(xmin=0, xmax=120, xlabel="m_{T} [GeV]", dology=False, ymax=4e4, donormalizebin=False), weightname = "weight_wplus_passIso")
    sampMan.cacheDraw("mT", "histo_wjets_muminusnu_mtcorr_passIso", mt_bins, DrawConfig(xmin=0, xmax=120, xlabel="m_{T} [GeV]", dology=False, ymax=3.5e4, donormalizebin=False), weightname = "weight_wminus_passIso")
    for iso in ["IsoCR1", "IsoCR2", "IsoCR3", "IsoCR4", "IsoCR5", "IsoCR6"]:
        sampMan.cacheDraw("mT", "histo_wjets_munu_mtcorr_" + iso, mt_bins, DrawConfig(xmin=0, xmax=120, xlabel="m_{T} [GeV]", dology=False, ymax=7e4, donormalizebin=False), weightname = "weight_"+ iso)
        sampMan.cacheDraw("mT", "histo_wjets_muplusnu_mtcorr_" + iso, mt_bins, DrawConfig(xmin=0, xmax=120, xlabel="m_{T} [GeV]", dology=False, ymax=7e4, donormalizebin=False), weightname = "weight_wplus_"+ iso)
        sampMan.cacheDraw("mT", "histo_wjets_muminusnu_mtcorr_" + iso, mt_bins, DrawConfig(xmin=0, xmax=120, xlabel="m_{T} [GeV]", dology=False, ymax=7e4, donormalizebin=False), weightname = "weight_wminus_"+ iso)

    # W pt
    sampMan.cacheDraw("WpT", "histo_wjets_munu_wpt", u_bins, DrawConfig(xmin=0., xmax=60, xlabel="p^{W}_{T} [GeV]", dology=False, ymax=1.1e5))
    #sampMan.cacheDraw("WpT", "histo_wjets_muplusnu_wpt", u_bins, DrawConfig(xmin=0., xmax=60, xlabel="p^{W}_{T} [GeV]", dology=False, ymax=7e4), weightname = "weight_wplus_passIso")
    #sampMan.cacheDraw("WpT", "histo_wjets_muminusnu_wpt", u_bins, DrawConfig(xmin=0., xmax=60, xlabel="p^{W}_{T} [GeV]", dology=False, ymax=5.0e4), weightname = "weight_wminus_passIso")

    sampMan.cacheDraw("Wphi", "histo_wjets_munu_wphi", 30, phimin, phimax, DrawConfig(xmin=phimin, xmax=phimax,      xlabel="#phi^{W}", dology=False, ymax=6e5))
    #sampMan.cacheDraw("Wphi", "histo_wjets_muplusnu_wphi", 30, phimin, phimax, DrawConfig(xmin=phimin, xmax=phimax,  xlabel="#phi^{W}", dology=False, ymax=4e5), weightname = "weight_wplus_passIso")
    #sampMan.cacheDraw("Wphi", "histo_wjets_muminusnu_wphi", 30, phimin, phimax, DrawConfig(xmin=phimin, xmax=phimax, xlabel="#phi^{W}", dology=False, ymax=3e5), weightname = "weight_wminus_passIso")

    # lep pt / mT
    sampMan.cacheDraw("pTOvermT", "histo_wjets_munu_pTOvermT", ptOvermT_bins, DrawConfig(ximin=0, xmax=3, xlabel="p^{#mu}_{T}/m_{T}", dology=True, ymax=1e8))

    sampMan.launchDraw()

    sampMan.dumpCounts()

    # save the histograms for fitting
    output = ROOT.TFile("result.root", "recreate")
    for channel in ["munu", "muplusnu", "muminusnu"]:
        for iso in ["passIso", "IsoCR1", "IsoCR2", "IsoCR3", "IsoCR4", "IsoCR5", "IsoCR6"]:
            hdata = sampMan.hdatas["histo_wjets_" + channel + "_mtcorr_" + iso]
            hmc   = THStack2TH1(sampMan.hsmcs[  "histo_wjets_" + channel + "_mtcorr_" + iso])
            hdata.SetDirectory(output)
            hdata.SetName("histo_wjets_" + channel + "_mtCorr_Data_" + iso)
            hdata.Write()
            hmc.SetDirectory(output)
            hmc.SetName("histo_wjets_" + channel + "_mtCorr_MC_" + iso)
            hmc.Write()
    output.Close()

    print "Program end..."

    raw_input()
    
    return 

if __name__ == "__main__":
   main()

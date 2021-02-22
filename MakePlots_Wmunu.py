import ROOT
import numpy as np
from collections import OrderedDict
import sys
sys.path.append("/afs/cern.ch/work/y/yofeng/public/CMSPLOTS")
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
    wjetsnorm = 1.06
    Wm0Samp   = Sample(input_wm0, isMC=True, name = "wm0", isWSR=True, additionalnorm = wjetsnorm)
    Wm1Samp   = Sample(input_wm1, isMC=True, name = "wm1", isWSR=True, additionalnorm = wjetsnorm)
    Wm2Samp   = Sample(input_wm2, isMC=True, name = "wm2", isWSR=True, additionalnorm = wjetsnorm)
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

    ## for the QCD background estimation (data-driven)
    qcdnorm = 1.0243
    DataAisoSamp  = Sample(input_antiiso_data, isMC=False, name="Data_aiso", isWSR=True, additionalnorm = qcdnorm, legend = 'QCD', color='226')
    # W -> munu
    Wm0AisoSamp   = Sample(input_antiiso_wm0, isMC=True, name = "wm0_aiso", isWSR=True, additionalnorm=-1.0 * qcdnorm)
    Wm1AisoSamp   = Sample(input_antiiso_wm1, isMC=True, name = "wm1_aiso", isWSR=True, additionalnorm=-1.0 * qcdnorm)
    Wm2AisoSamp   = Sample(input_antiiso_wm2, isMC=True, name = "wm2_aiso", isWSR=True, additionalnorm=-1.0 * qcdnorm)
    # ttbar
    TTbarAisoSamp  = Sample(input_antiiso_ttbar, isMC=True, name = "ttbar_dilepton_aiso", isWSR=True, additionalnorm=-1.0 * qcdnorm)
    TT1LepAisoSamp = Sample(input_antiiso_ttbar_1lep, isMC=True, name = "ttbar_1lepton_aiso", isWSR=True, additionalnorm=-1.0 * qcdnorm)
    TT0LepAisoSamp = Sample(input_antiiso_ttbar_0lep, isMC=True, name = "ttbar_0lepton_aiso", isWSR=True, additionalnorm=-1.0 * qcdnorm)
    ## dibosons
    WWAisoSamp = Sample(input_antiiso_ww, isMC=True, name = "WW_aiso", isWSR=True, additionalnorm=-1.0 * qcdnorm)
    WZAisoSamp = Sample(input_antiiso_wz, isMC=True, name = "WZ_aiso", isWSR=True, additionalnorm=-1.0 * qcdnorm)
    ZZAisoSamp = Sample(input_antiiso_zz, isMC=True, name = "ZZ_aiso", isWSR=True, additionalnorm=-1.0 * qcdnorm)
    # tau
    ZXXAisoSamp = Sample(input_antiiso_zxx, isMC=True, name = "ZXX_aiso", isWSR=True, additionalnorm=-1.0 * qcdnorm)
    Wx0AisoSamp = Sample(input_antiiso_wx0, isMC=True, name = "wx0_aiso", isWSR=True, additionalnorm=-1.0 * qcdnorm)
    Wx1AisoSamp = Sample(input_antiiso_wx1, isMC=True, name = "wx1_aiso", isWSR=True, additionalnorm=-1.0 * qcdnorm)
    Wx2AisoSamp = Sample(input_antiiso_wx2, isMC=True, name = "wx2_aiso", isWSR=True, additionalnorm=-1.0 * qcdnorm)

    sampMan = SampleManager(DataSamp, [Wm0Samp, Wm1Samp, Wm2Samp, TTbarSamp, TT1LepSamp, TT0LepSamp, WWSamp, WZSamp, ZZSamp, ZXXSamp, Wx0Samp, Wx1Samp, Wx2Samp, 
            DataAisoSamp, Wm0AisoSamp, Wm1AisoSamp, Wm2AisoSamp, TTbarAisoSamp, TT1LepAisoSamp, TT0LepAisoSamp, WWAisoSamp, WZAisoSamp, ZZAisoSamp, ZXXAisoSamp, Wx0AisoSamp, Wx1AisoSamp, Wx2AisoSamp])
    #sampMan = SampleManager(DataSamp, [Wm0Samp, Wm1Samp, Wm2Samp, TTbarSamp, TT1LepSamp, TT0LepSamp, WWSamp, WZSamp, ZZSamp, ZXXSamp, Wx0Samp, Wx1Samp, Wx2Samp,
    #        DataAisoSamp])
    sampMan.groupMCs(["WW", "WZ", "ZZ", "ZXX", "wx0", "wx1", "wx2"], "EWK", 216, "EWK")
    sampMan.groupMCs(["ttbar_dilepton", "ttbar_1lepton", "ttbar_0lepton"], "ttbar", 96, "t#bar{t}")
    sampMan.groupMCs(['wm0', 'wm1', 'wm2'], "wmunu", 92,"W#rightarrow#mu#nu")
    sampMan.groupMCs(['Data_aiso', 'wm0_aiso', 'wm1_aiso', 'wm2_aiso', 'ttbar_dilepton_aiso', 'ttbar_1lepton_aiso', 'ttbar_0lepton_aiso', 'WW_aiso', 'WZ_aiso', 'ZZ_aiso', 'ZXX_aiso', 'wx0_aiso', 'wx1_aiso' ,'wx2_aiso'], 'QCD', 226, 'QCD')

    sampMan.ApplyCutSpecificMCs('relIso > 0.35 && relIso < 0.45', sampgroupnames=['QCD'])

    #sampMan.DefineAll("zpt", "Z.Pt()")
    #sampMan.DefineAll("zmass", "ZMass")
    #sampMan.DefineAll("zy",  "Z.Rapidity()")
    #sampMan.DefineAll("zeta",  "Z.Rapidity()")
    #sampMan.DefineAll("nPV", "npv")

    sampMan.DefineAll("Muon_pt", "lep.Pt()")
    sampMan.DefineAll("Muon_eta", "lep.Eta()")
    sampMan.DefineMC("met_pt", "metVars[1]", excludes=['Data_aiso'])
    sampMan.DefineMC("met_phi", "TVector2::Phi_mpi_pi(metVarsPhi[1])", excludes=['Data_aiso'])
    # Data does not run any recoil corrections
    DataSamp.Define("met_pt", "metVars[0]")
    DataSamp.Define("met_phi", "TVector2::Phi_mpi_pi(metVarsPhi[0])")
    DataAisoSamp.Define("met_pt", "metVars[0]")
    DataAisoSamp.Define("met_phi", "TVector2::Phi_mpi_pi(metVarsPhi[0])")
    sampMan.DefineAll("met_raw_pt", "metVars[0]")
    sampMan.DefineAll("met_raw_phi", "TVector2::Phi_mpi_pi(metVarsPhi[0])")

    sampMan.DefineAll("weight_wplus", "weight_WoVpt * (q>0)", excludeGroups=['QCD'])
    sampMan.DefineAll("weight_wminus", "weight_WoVpt * (q<0)", excludeGroups=['QCD'])
    sampMan.DefineSpecificMCs("weight_wplus", "weight_WoVpt * (q>0)", sampgroupnames=['QCD'])
    sampMan.DefineSpecificMCs("weight_wminus", "weight_WoVpt * (q<0)", sampgroupnames=['QCD'])

    # deltaPhi
    sampMan.DefineAll("deltaPhi", "TVector2::Phi_0_2pi(lep.Phi() - met_phi)")

    # recoil variables
    sampMan.DefineAll("V2W", "UVec(lep.Pt(), lep.Phi(), met_pt, met_phi)")
    sampMan.DefineAll("WpT", "V2W.Mod()")
    sampMan.DefineAll("Wphi", "TVector2::Phi_mpi_pi(V2W.Phi())")


    met_pt_bins = np.array([0., 2.0, 4., 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 33, 36, 39, 42, 45, 48, 51, 55, 60, 65, 70, 75, 80])
    u1_bins = np.array([-40.,-36.,-32., -28., -25., -22.0, -20.0, -18, -16, -14, -12, -10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 53, 56, 59, 64, 68, 72, 76, 80, 85, 90, 100])
    u2_bins = np.array([-80., -70., -65., -60., -56., -52, -48, -44, -40, -37, -34, -31, -28, -25., -22., -20, -18, -16, -14, -12, -10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10, 12, 16, 18, 20, 22, 25, 28, 31, 34, 37, 40, 44, 48, 52, 56, 60, 65, 70, 80])
    #u_bins = np.array([0., 2., 4., 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 43, 46, 49, 52, 56, 60, 64, 68, 72, 76, 80, 85, 90, 95, 100, 105, 110, 115, 120, 130, 140, 150])
    mt_bins = np.concatenate((np.linspace(40, 90, 51), np.linspace(92, 100, 5), np.linspace(104, 120, 5)))
    phimin = -ROOT.TMath.Pi()
    phimax = ROOT.TMath.Pi()
    dphi_bins = np.concatenate((np.linspace(0,1, 3), np.linspace(1.2, 2.4, 7),np.linspace(2.5, 3.5, 11), np.linspace(3.6, 4.8, 7), np.linspace(5.0, 6.0, 3))) * ROOT.TMath.Pi() / 3

    eta_bins = np.concatenate((np.linspace(-2.4, -2.0, 3),np.linspace(-1.8,1.8,37), np.linspace(2.0, 2.4, 3)))
    pt_bins = np.concatenate((np.linspace(25, 35, 6),np.linspace(36,55,20), np.linspace(57, 63, 4), np.linspace(66,70,2)))
    mass_bins = np.concatenate((np.linspace(70, 76, 3),np.linspace(78,86,5), np.linspace(87, 96, 10), np.linspace(98,104,5), np.linspace(106, 110, 2)))
    u_bins = np.concatenate((np.linspace(0, 20, 11),np.linspace(24,80,15), np.linspace(85, 109, 4), np.linspace(120,150,3)))
    u_bins = np.concatenate((np.linspace(0, 20, 11),np.linspace(24,60,10)))

    #sampMan.cacheDraw("zpt",   "histo_wjets_zpt_WoZptWeight_munu", u_bins, DrawConfig(xmin=0, xmax=150, xlabel='p_{T}^{#mu#mu} [GeV]'), weightname="weight_WoVpt")
    #sampMan.cacheDraw("zmass", "histo_wjets_zmass_munu", mass_bins, DrawConfig(xmin=70, xmax=110, xlabel='m_{#mu#mu} [GeV]'))
    #sampMan.cacheDraw("zeta",  "histo_wjets_zeta_munu", eta_bins, DrawConfig(xmin=-2.5, xmax=2.5, xlabel='#eta_{#mu#mu} [GeV]', ymax=1e7, ylabel='Events / 1'))
    #sampMan.cacheDraw("zy",    "histo_wjets_zrapidity_munu", eta_bins, DrawConfig(xmin=-2.5, xmax=2.5, xlabel='y_{#mu#mu} [GeV]', ymax=1e7, ylabel='Events / 1'))

    # nPV
    sampMan.cacheDraw("npv", "histo_nPV_munu", 10, 0, 10, DrawConfig(xmin=0, xmax=10, xlabel='# PV', ylabel='Events / 1', dology=False, ymax = 9e5))
    sampMan.cacheDraw("npv", "histo_nPV_muplusnu", 10, 0, 10, DrawConfig(xmin=0, xmax=10, xlabel='# PV', ylabel='Events / 1', dology=False, ymax = 5e5), weightname = "weight_wplus")
    sampMan.cacheDraw("npv", "histo_nPV_muminusnu", 10, 0, 10, DrawConfig(xmin=0, xmax=10, xlabel='# PV', ylabel='Events / 1', dology=False, ymax = 4.5e5), weightname = "weight_wminus")
    # lepton pt
    sampMan.cacheDraw("Muon_pt", "histo_Muon_pt_munu", pt_bins, DrawConfig(xmin=20, xmax=70, xlabel='p_{T}(#mu) [GeV]', dology=False, ymax=1.05e5))
    sampMan.cacheDraw("Muon_eta", "histo_Muon_eta_munu", eta_bins, DrawConfig(xmin=-2.6, xmax=2.6, xlabel='#eta (#mu) [GeV]', ylabel='Events / 1', dology=False, ymax=6e5))
    sampMan.cacheDraw("Muon_pt", "histo_Muon_pt_muplusnu", pt_bins, DrawConfig(xmin=20, xmax=70, xlabel='p_{T}(#mu) [GeV]', dology=False, ymax=6.05e4), weightname = "weight_wplus")
    sampMan.cacheDraw("Muon_eta", "histo_Muon_eta_muplusnu", eta_bins, DrawConfig(xmin=-2.6, xmax=2.6, xlabel='#eta (#mu) [GeV]', ylabel='Events / 1', dology=False, ymax=4e5), weightname = "weight_wminus")
    sampMan.cacheDraw("Muon_pt", "histo_Muon_pt_muminusnu", pt_bins, DrawConfig(xmin=20, xmax=70, xlabel='p_{T}(#mu) [GeV]', dology=False, ymax=5.05e4), weightname = "weight_wminus")
    sampMan.cacheDraw("Muon_eta", "histo_Muon_eta_muminusnu", eta_bins, DrawConfig(xmin=-2.6, xmax=2.6, xlabel='#eta (#mu) [GeV]', ylabel='Events / 1', dology=False, ymax=3.5e5), weightname = "weight_wminus")

    # met and raw_met
    sampMan.cacheDraw("met_pt", "histo_wjets_munu_pfmet_pt", met_pt_bins, DrawConfig(xmin=0, xmax=80, xlabel='PF MET [GeV]', dology=False, ymax=8e4))
    sampMan.cacheDraw("met_phi", "histo_wjets_munu_pfmet_phi", 30, phimin, phimax, DrawConfig(xmin=phimin, xmax=phimax, xlabel='PF MET #phi', dology=False, ymax=6e5))
    sampMan.cacheDraw("met_raw_pt", "histo_wjets_munu_pfmet_raw_pt", met_pt_bins, DrawConfig(xmin=0, xmax=80, xlabel='Raw PF MET [GeV]', dology=False, ymax=8e4))
    sampMan.cacheDraw("met_raw_phi", "histo_wjets_munu_pfmet_raw_phi", 30, phimin, phimax, DrawConfig(xmin=phimin, xmax=phimax, xlabel='Raw PF MET #phi', dology=False, ymax=6e5))

    sampMan.cacheDraw("met_pt", "histo_wjets_muplusnu_pfmet_pt", met_pt_bins, DrawConfig(xmin=0, xmax=80, xlabel='PF MET [GeV]', dology=False, ymax=5e4), weightname = "weight_wplus")
    sampMan.cacheDraw("met_phi", "histo_wjets_muplusnu_pfmet_phi", 30, phimin, phimax, DrawConfig(xmin=phimin, xmax=phimax, xlabel='PF MET #phi', dology=False, ymax=4e5), weightname = "weight_wplus")
    sampMan.cacheDraw("met_raw_pt", "histo_wjets_muplusnu_pfmet_raw_pt", met_pt_bins, DrawConfig(xmin=0, xmax=80, xlabel='Raw PF MET [GeV]', dology=False, ymax=5e4), weightname = "weight_wplus")
    sampMan.cacheDraw("met_raw_phi", "histo_wjets_muplusnu_pfmet_raw_phi", 30, phimin, phimax, DrawConfig(xmin=phimin, xmax=phimax, xlabel='Raw PF MET #phi', dology=False, ymax=3.5e5), weightname = "weight_wplus")

    sampMan.cacheDraw("met_pt", "histo_wjets_muminusnu_pfmet_pt", met_pt_bins, DrawConfig(xmin=0, xmax=80, xlabel='PF MET [GeV]', dology=False, ymax=4e4), weightname = "weight_wminus")
    sampMan.cacheDraw("met_phi", "histo_wjets_muminusnu_pfmet_phi", 30, phimin, phimax, DrawConfig(xmin=phimin, xmax=phimax, xlabel='PF MET #phi', dology=False, ymax=3e5), weightname = "weight_wminus")
    sampMan.cacheDraw("met_raw_pt", "histo_wjets_muminusnu_pfmet_raw_pt", met_pt_bins, DrawConfig(xmin=0, xmax=80, xlabel='Raw PF MET [GeV]', dology=False, ymax=4e4), weightname = "weight_wminus")
    sampMan.cacheDraw("met_raw_phi", "histo_wjets_muminusnu_pfmet_raw_phi", 30, phimin, phimax, DrawConfig(xmin=phimin, xmax=phimax, xlabel='Raw PF MET #phi', dology=False, ymax=3e5), weightname = "weight_wminus")

    # deltaPhi
    sampMan.cacheDraw("deltaPhi", "histo_wjets_munu_deltaphi",      dphi_bins,  DrawConfig(xmin=0., xmax=2*phimax, xlabel='#Delta #phi(#mu, MET)', dology=False, ymax=2e6))
    sampMan.cacheDraw("deltaPhi", "histo_wjets_muplusnu_deltaphi",  dphi_bins,  DrawConfig(xmin=0., xmax=2*phimax, xlabel='#Delta #phi(#mu, MET)', dology=False, ymax=1.4e6), weightname = "weight_wplus")
    sampMan.cacheDraw("deltaPhi", "histo_wjets_muminusnu_deltaphi", dphi_bins,  DrawConfig(xmin=0., xmax=2*phimax, xlabel='#Delta #phi(#mu, MET)', dology=False, ymax=1.4e6), weightname = "weight_wminus")

    # mt
    sampMan.cacheDraw("mtCorr", "histo_wjets_munu_mtcorr",      mt_bins, DrawConfig(xmin=40, xmax=120, xlabel="m_{T} [GeV]", dology=False, ymax=7e4))
    sampMan.cacheDraw("mtCorr", "histo_wjets_muplusnu_mtcorr",  mt_bins, DrawConfig(xmin=40, xmax=120, xlabel="m_{T} [GeV]", dology=False, ymax=4e4), weightname = "weight_wplus")
    sampMan.cacheDraw("mtCorr", "histo_wjets_muminusnu_mtcorr", mt_bins, DrawConfig(xmin=40, xmax=120, xlabel="m_{T} [GeV]", dology=False, ymax=3.5e4), weightname = "weight_wminus")

    # W pt
    sampMan.cacheDraw("WpT", "histo_wjets_munu_wpt", u_bins, DrawConfig(xmin=0., xmax=60, xlabel="p^{W}_{T} [GeV]", dology=False, ymax=1.1e5))
    sampMan.cacheDraw("WpT", "histo_wjets_muplusnu_wpt", u_bins, DrawConfig(xmin=0., xmax=60, xlabel="p^{W}_{T} [GeV]", dology=False, ymax=7e4), weightname = "weight_wplus")
    sampMan.cacheDraw("WpT", "histo_wjets_muminusnu_wpt", u_bins, DrawConfig(xmin=0., xmax=60, xlabel="p^{W}_{T} [GeV]", dology=False, ymax=5.0e4), weightname = "weight_wminus")

    sampMan.cacheDraw("Wphi", "histo_wjets_munu_wphi", 30, phimin, phimax, DrawConfig(xmin=phimin, xmax=phimax,      xlabel="#phi^{W}", dology=False, ymax=6e5))
    sampMan.cacheDraw("Wphi", "histo_wjets_muplusnu_wphi", 30, phimin, phimax, DrawConfig(xmin=phimin, xmax=phimax,  xlabel="#phi^{W}", dology=False, ymax=4e5), weightname = "weight_wplus")
    sampMan.cacheDraw("Wphi", "histo_wjets_muminusnu_wphi", 30, phimin, phimax, DrawConfig(xmin=phimin, xmax=phimax, xlabel="#phi^{W}", dology=False, ymax=3e5), weightname = "weight_wminus")

    sampMan.launchDraw()

    sampMan.dumpCounts()

    print "Program end..."

    raw_input()
    
    return 

if __name__ == "__main__":
   main()

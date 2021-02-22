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

    input_data    = "inputs/wenu/input_data.txt"
    input_we0     = "inputs/wenu/input_we0.txt"
    input_we1     = "inputs/wenu/input_we1.txt"
    input_we2     = "inputs/wenu/input_we2.txt"
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

    input_antiiso_data    = "inputs/awenu/input_data.txt"
    input_antiiso_we0     = "inputs/awenu/input_we0.txt"
    input_antiiso_we1     = "inputs/awenu/input_we1.txt"
    input_antiiso_we2     = "inputs/awenu/input_we2.txt"
    input_antiiso_ttbar   = "inputs/awenu/input_ttbar_dilepton.txt"
    input_antiiso_ttbar_1lep = "inputs/awenu/input_ttbar_singlelepton.txt"
    input_antiiso_ttbar_0lep = "inputs/awenu/input_ttbar_hadronic.txt"
    input_antiiso_ww      = "inputs/awenu/input_ww.txt"
    input_antiiso_wz      = "inputs/awenu/input_wz.txt"
    input_antiiso_zz      = "inputs/awenu/input_zz.txt"
    input_antiiso_zxx     = "inputs/awenu/input_zxx.txt"
    input_antiiso_wx0     = "inputs/awenu/input_wx0.txt"
    input_antiiso_wx1     = "inputs/awenu/input_wx1.txt"
    input_antiiso_wx2     = "inputs/awenu/input_wx2.txt"

    DataSamp  = Sample(input_data, isMC=False, legend="Data", name="Data", isWSR=True)
    # W -> enu
    wjetsnorm = 1.075
    Wm0Samp   = Sample(input_we0, isMC=True, name = "we0", isWSR=True, additionalnorm = wjetsnorm)
    Wm1Samp   = Sample(input_we1, isMC=True, name = "we1", isWSR=True, additionalnorm = wjetsnorm)
    Wm2Samp   = Sample(input_we2, isMC=True, name = "we2", isWSR=True, additionalnorm = wjetsnorm)
    # ttbar
    ttbarnorm = 1.2054
    TTbarSamp  = Sample(input_ttbar, isMC=True, name = "ttbar_dilepton", isWSR=True, additionalnorm = ttbarnorm)
    TT1LepSamp = Sample(input_ttbar_1lep, isMC=True, name = "ttbar_1lepton", isWSR=True, additionalnorm = ttbarnorm)
    TT0LepSamp = Sample(input_ttbar_0lep, isMC=True, name = "ttbar_0lepton", isWSR=True, additionalnorm = ttbarnorm)
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
    qcdnorm = 16.7254
    DataAisoSamp  = Sample(input_antiiso_data, isMC=False, name="Data_aiso", isWSR=True, additionalnorm = qcdnorm, legend = 'QCD', color='226')
    # W -> enu
    Wm0AisoSamp   = Sample(input_antiiso_we0, isMC=True, name = "we0_aiso", isWSR=True, additionalnorm=-1.0 * qcdnorm)
    Wm1AisoSamp   = Sample(input_antiiso_we1, isMC=True, name = "we1_aiso", isWSR=True, additionalnorm=-1.0 * qcdnorm)
    Wm2AisoSamp   = Sample(input_antiiso_we2, isMC=True, name = "we2_aiso", isWSR=True, additionalnorm=-1.0 * qcdnorm)
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
    sampMan.groupMCs(['we0', 'we1', 'we2'], "wenu", 92,"W#rightarrow e#nu")
    sampMan.groupMCs(['Data_aiso', 'we0_aiso', 'we1_aiso', 'we2_aiso', 'ttbar_dilepton_aiso', 'ttbar_1lepton_aiso', 'ttbar_0lepton_aiso', 'WW_aiso', 'WZ_aiso', 'ZZ_aiso', 'ZXX_aiso', 'wx0_aiso', 'wx1_aiso' ,'wx2_aiso'], 'QCD', 226, 'QCD')

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

    #sampMan.cacheDraw("zpt",   "histo_wjets_zpt_WoZptWeight_enu", u_bins, DrawConfig(xmin=0, xmax=150, xlabel='p_{T}^{#mu#mu} [GeV]'), weightname="weight_WoVpt")
    #sampMan.cacheDraw("zmass", "histo_wjets_zmass_enu", mass_bins, DrawConfig(xmin=70, xmax=110, xlabel='m_{#mu#mu} [GeV]'))
    #sampMan.cacheDraw("zeta",  "histo_wjets_zeta_enu", eta_bins, DrawConfig(xmin=-2.5, xmax=2.5, xlabel='#eta_{#mu#mu} [GeV]', ymax=1e7, ylabel='Events / 1'))
    #sampMan.cacheDraw("zy",    "histo_wjets_zrapidity_enu", eta_bins, DrawConfig(xmin=-2.5, xmax=2.5, xlabel='y_{#mu#mu} [GeV]', ymax=1e7, ylabel='Events / 1'))

    sf = 0.7

    # nPV
    sampMan.cacheDraw("npv", "histo_nPV_enu", 10, 0, 10, DrawConfig(xmin=0, xmax=10, xlabel='# PV', ylabel='Events / 1', dology=False, ymax = 9e5 * sf))
    sampMan.cacheDraw("npv", "histo_nPV_eplusnu", 10, 0, 10, DrawConfig(xmin=0, xmax=10, xlabel='# PV', ylabel='Events / 1', dology=False, ymax = 5e5 * sf), weightname = "weight_wplus")
    sampMan.cacheDraw("npv", "histo_nPV_eminusnu", 10, 0, 10, DrawConfig(xmin=0, xmax=10, xlabel='# PV', ylabel='Events / 1', dology=False, ymax = 4.5e5 * sf), weightname = "weight_wminus")
    # lepton pt
    sampMan.cacheDraw("Muon_pt", "histo_Muon_pt_enu", pt_bins, DrawConfig(xmin=20, xmax=70, xlabel='p_{T}(e) [GeV]', dology=False, ymax=1.05e5 * sf))
    sampMan.cacheDraw("Muon_eta", "histo_Muon_eta_enu", eta_bins, DrawConfig(xmin=-2.6, xmax=2.6, xlabel='#eta (e) [GeV]', ylabel='Events / 1', dology=False, ymax=6e5 * sf))
    sampMan.cacheDraw("Muon_pt", "histo_Muon_pt_eplusnu", pt_bins, DrawConfig(xmin=20, xmax=70, xlabel='p_{T}(e) [GeV]', dology=False, ymax=6.05e4 * sf), weightname = "weight_wplus")
    sampMan.cacheDraw("Muon_eta", "histo_Muon_eta_eplusnu", eta_bins, DrawConfig(xmin=-2.6, xmax=2.6, xlabel='#eta (e) [GeV]', ylabel='Events / 1', dology=False, ymax=4e5 * sf), weightname = "weight_wminus")
    sampMan.cacheDraw("Muon_pt", "histo_Muon_pt_eminusnu", pt_bins, DrawConfig(xmin=20, xmax=70, xlabel='p_{T}(e) [GeV]', dology=False, ymax=5.05e4 * sf), weightname = "weight_wminus")
    sampMan.cacheDraw("Muon_eta", "histo_Muon_eta_eminusnu", eta_bins, DrawConfig(xmin=-2.6, xmax=2.6, xlabel='#eta (e) [GeV]', ylabel='Events / 1', dology=False, ymax=3.5e5 * sf), weightname = "weight_wminus")

    # met and raw_met
    sampMan.cacheDraw("met_pt", "histo_wjets_enu_pfmet_pt", met_pt_bins, DrawConfig(xmin=0, xmax=80, xlabel='PF MET [GeV]', dology=False, ymax=8e4 * sf ))
    sampMan.cacheDraw("met_phi", "histo_wjets_enu_pfmet_phi", 30, phimin, phimax, DrawConfig(xmin=phimin, xmax=phimax, xlabel='PF MET #phi', dology=False, ymax=6e5 * sf))
    sampMan.cacheDraw("met_raw_pt", "histo_wjets_enu_pfmet_raw_pt", met_pt_bins, DrawConfig(xmin=0, xmax=80, xlabel='Raw PF MET [GeV]', dology=False, ymax=8e4 * sf))
    sampMan.cacheDraw("met_raw_phi", "histo_wjets_enu_pfmet_raw_phi", 30, phimin, phimax, DrawConfig(xmin=phimin, xmax=phimax, xlabel='Raw PF MET #phi', dology=False, ymax=6e5 * sf))

    sampMan.cacheDraw("met_pt", "histo_wjets_eplusnu_pfmet_pt", met_pt_bins, DrawConfig(xmin=0, xmax=80, xlabel='PF MET [GeV]', dology=False, ymax=5e4 * sf), weightname = "weight_wplus")
    sampMan.cacheDraw("met_phi", "histo_wjets_eplusnu_pfmet_phi", 30, phimin, phimax, DrawConfig(xmin=phimin, xmax=phimax, xlabel='PF MET #phi', dology=False, ymax=4e5 * sf), weightname = "weight_wplus")
    sampMan.cacheDraw("met_raw_pt", "histo_wjets_eplusnu_pfmet_raw_pt", met_pt_bins, DrawConfig(xmin=0, xmax=80, xlabel='Raw PF MET [GeV]', dology=False, ymax=5e4 * sf), weightname = "weight_wplus")
    sampMan.cacheDraw("met_raw_phi", "histo_wjets_eplusnu_pfmet_raw_phi", 30, phimin, phimax, DrawConfig(xmin=phimin, xmax=phimax, xlabel='Raw PF MET #phi', dology=False, ymax=3.5e5 * sf), weightname = "weight_wplus")

    sampMan.cacheDraw("met_pt", "histo_wjets_eminusnu_pfmet_pt", met_pt_bins, DrawConfig(xmin=0, xmax=80, xlabel='PF MET [GeV]', dology=False, ymax=4e4 * sf), weightname = "weight_wminus")
    sampMan.cacheDraw("met_phi", "histo_wjets_eminusnu_pfmet_phi", 30, phimin, phimax, DrawConfig(xmin=phimin, xmax=phimax, xlabel='PF MET #phi', dology=False, ymax=3e5 * sf), weightname = "weight_wminus")
    sampMan.cacheDraw("met_raw_pt", "histo_wjets_eminusnu_pfmet_raw_pt", met_pt_bins, DrawConfig(xmin=0, xmax=80, xlabel='Raw PF MET [GeV]', dology=False, ymax=4e4 * sf), weightname = "weight_wminus")
    sampMan.cacheDraw("met_raw_phi", "histo_wjets_eminusnu_pfmet_raw_phi", 30, phimin, phimax, DrawConfig(xmin=phimin, xmax=phimax, xlabel='Raw PF MET #phi', dology=False, ymax=3e5 * sf), weightname = "weight_wminus")

    # deltaPhi
    sampMan.cacheDraw("deltaPhi", "histo_wjets_enu_deltaphi",      dphi_bins,  DrawConfig(xmin=0., xmax=2*phimax, xlabel='#Delta #phi(e, MET)', dology=False, ymax=2e6 * sf))
    sampMan.cacheDraw("deltaPhi", "histo_wjets_eplusnu_deltaphi",  dphi_bins,  DrawConfig(xmin=0., xmax=2*phimax, xlabel='#Delta #phi(e, MET)', dology=False, ymax=1.4e6 * sf), weightname = "weight_wplus")
    sampMan.cacheDraw("deltaPhi", "histo_wjets_eminusnu_deltaphi", dphi_bins,  DrawConfig(xmin=0., xmax=2*phimax, xlabel='#Delta #phi(e, MET)', dology=False, ymax=1.4e6 *sf), weightname = "weight_wminus")

    # mt
    sampMan.cacheDraw("mtCorr", "histo_wjets_enu_mtcorr",      mt_bins, DrawConfig(xmin=40, xmax=120, xlabel="m_{T} [GeV]", dology=False, ymax=7e4 *sf))
    sampMan.cacheDraw("mtCorr", "histo_wjets_eplusnu_mtcorr",  mt_bins, DrawConfig(xmin=40, xmax=120, xlabel="m_{T} [GeV]", dology=False, ymax=4e4 *sf), weightname = "weight_wplus")
    sampMan.cacheDraw("mtCorr", "histo_wjets_eminusnu_mtcorr", mt_bins, DrawConfig(xmin=40, xmax=120, xlabel="m_{T} [GeV]", dology=False, ymax=3.5e4 *sf ), weightname = "weight_wminus")

    # W pt
    sampMan.cacheDraw("WpT", "histo_wjets_enu_wpt", u_bins, DrawConfig(xmin=0., xmax=60, xlabel="p^{W}_{T} [GeV]", dology=False, ymax=1.1e5*sf))
    sampMan.cacheDraw("WpT", "histo_wjets_eplusnu_wpt", u_bins, DrawConfig(xmin=0., xmax=60, xlabel="p^{W}_{T} [GeV]", dology=False, ymax=7e4*sf), weightname = "weight_wplus")
    sampMan.cacheDraw("WpT", "histo_wjets_eminusnu_wpt", u_bins, DrawConfig(xmin=0., xmax=60, xlabel="p^{W}_{T} [GeV]", dology=False, ymax=5.0e4*sf), weightname = "weight_wminus")

    sampMan.cacheDraw("Wphi", "histo_wjets_enu_wphi", 30, phimin, phimax, DrawConfig(xmin=phimin, xmax=phimax,      xlabel="#phi^{W}", dology=False, ymax=6e5*sf))
    sampMan.cacheDraw("Wphi", "histo_wjets_eplusnu_wphi", 30, phimin, phimax, DrawConfig(xmin=phimin, xmax=phimax,  xlabel="#phi^{W}", dology=False, ymax=4e5*sf), weightname = "weight_wplus")
    sampMan.cacheDraw("Wphi", "histo_wjets_eminusnu_wphi", 30, phimin, phimax, DrawConfig(xmin=phimin, xmax=phimax, xlabel="#phi^{W}", dology=False, ymax=3e5*sf), weightname = "weight_wminus")

    sampMan.launchDraw()

    sampMan.dumpCounts()

    print "Program end..."

    raw_input()
    
    return 

if __name__ == "__main__":
   main()

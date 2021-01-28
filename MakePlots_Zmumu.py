'''
plot the data-MC comparisons pre/post DeepMET corrections.
and the systematic uncs.
'''
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
ROOT.ROOT.EnableImplicitMT(10)

dotest = 0
VetoB = False
doSys = False

def main():
    print "Program start..."

    #ROOT.gROOT.ProcessLine('TFile* f_zpt = TFile::Open("results/zpt_weight.root")')
    #ROOT.gROOT.ProcessLine('TH1D* h_zpt_ratio  = (TH1D*)f_zpt->Get("h_zpt_ratio")')

    if dotest:
        input_data    = "inputs/zmumu/input_data.txt"
        input_dy      = "inputs/zmumu/input_zjets.txt"
        input_ttbar   = "inputs/zmumu/input_ttbar_dilepton.txt"
    else:
        input_data    = "inputs/zmumu/input_data.txt"
        input_dy      = "inputs/zmumu/input_zjets.txt"
        input_ttbar   = "inputs/zmumu/input_ttbar_dilepton.txt"
        input_ttbar_1lep = "inputs/zmumu/input_ttbar_singlelepton.txt"
        input_ttbar_0lep = "inputs/zmumu/input_ttbar_hadronic.txt"
        input_ww      = "inputs/zmumu/input_ww.txt"
        input_wz      = "inputs/zmumu/input_wz.txt"
        input_zz      = "inputs/zmumu/input_zz.txt"
        input_zxx     = "inputs/zmumu/input_zxx.txt"

    DataSamp  = Sample(input_data, isMC=False, legend="Data", name="Data", bjetVeto=VetoB)
    DYSamp    = Sample(input_dy,    xsec = 6077.22*1e3,      color=92,  reweightzpt = True, legend="DY", name="DY", bjetVeto=VetoB, nmcevt=311461009200.830078)
    TTbarSamp = Sample(input_ttbar, xsec = 831.76*0.105*1e3, color=46, reweightzpt = False, legend="t#bar{t}", name="ttbar2lep", bjetVeto=VetoB, nmcevt=360411094.548295)
    if not dotest:
        WWSamp  = Sample(input_ww,  xsec = 12.178*1e3,       color=38, reweightzpt = False, legend="WW2L",     name="WW",  bjetVeto=VetoB, nmcevt=1000000.0)
        WZSamp  = Sample(input_wz,  xsec = 5.26*1e3,         color=39, reweightzpt = False, legend="WZ3L",     name="WZ",  bjetVeto=VetoB, nmcevt=885000.00)
        ZZSamp  = Sample(input_zz,  xsec = 0.564*1e3,        color=37, reweightzpt = False, legend="ZZ2L",     name="ZZ",  bjetVeto=VetoB, nmcevt=1000000.0)
        ZXXSamp = Sample(input_zxx, xsec = 1.0,              color=40, reweightzpt = False, legend="ZXX",      name="ZXX",   bjetVeto=VetoB, nmcevt=3.1146101e+11)
        TT1LepSamp = Sample(input_ttbar_1lep,  xsec = 831.76*0.219*1e3, color=47, reweightzpt = False, legend="t#bar{t}", name="ttbar1lep", bjetVeto=VetoB, nmcevt=1.2561491e+09)
        TT0LepSamp = Sample(input_ttbar_0lep,  xsec = 831.76*0.219*1e3, color=48, reweightzpt = False, legend="t#bar{t}", name="ttbar0lep", bjetVeto=VetoB, nmcevt=6.2414786e+09)

    if not dotest:
        sampMan = SampleManager(DataSamp, [DYSamp, TTbarSamp, WWSamp, WZSamp, ZZSamp, ZXXSamp, TT1LepSamp, TT0LepSamp])
    else:
        sampMan = SampleManager(DataSamp, [DYSamp, TTbarSamp])
    sampMan.groupMCs(["WW", "WZ", "ZZ", "ZXX"], "EWK", 216, "EWK")
    sampMan.groupMCs(["ttbar2lep", "ttbar1lep", "ttbar0lep"], "ttbar", 96, "t#bar{t}")

    sampMan.DefineAll("zpt", "Z.Pt()")
    sampMan.DefineAll("zmass", "ZMass")
    sampMan.DefineAll("zy",  "Z.Rapidity()")
    sampMan.DefineAll("zeta",  "Z.Rapidity()")
    sampMan.DefineAll("nPV", "npv")

    # define u1 and u2 for PUPPI and PF
    # u1 and u2 for DeepMET has already been defined in SampleManager
    sampMan.DefineAll("LeadMuon_pt", "lep1_corr.Pt()")
    sampMan.DefineAll("LeadMuon_eta", "lep1_corr.Eta()")
    sampMan.DefineAll("SubleadMuon_pt", "lep2_corr.Pt()")
    sampMan.DefineAll("SubleadMuon_eta", "lep2_corr.Eta()")

    #DYSamp.Define("u_pt_corr_central", "TMath::Sqrt(u1_corr_central*u1_corr_central + u2_corr_central*u2_corr_central)")
    #sampMan.DefineAll("u1_corr_central",     "u1"  , excludes=['DY'])


    met_pt_bins = np.array([0., 2.0, 4., 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 33, 36, 39, 42, 45, 48, 51, 55, 60, 65, 70, 75, 80, 90, 100, 110, 120, 135, 150, 165, 180, 200])
    u1_bins = np.array([-40.,-36.,-32., -28., -25., -22.0, -20.0, -18, -16, -14, -12, -10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 53, 56, 59, 64, 68, 72, 76, 80, 85, 90, 100])
    u2_bins = np.array([-80., -70., -65., -60., -56., -52, -48, -44, -40, -37, -34, -31, -28, -25., -22., -20, -18, -16, -14, -12, -10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10, 12, 16, 18, 20, 22, 25, 28, 31, 34, 37, 40, 44, 48, 52, 56, 60, 65, 70, 80])
    #u_bins = np.array([0., 2., 4., 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 43, 46, 49, 52, 56, 60, 64, 68, 72, 76, 80, 85, 90, 95, 100, 105, 110, 115, 120, 130, 140, 150])
    phimin = -ROOT.TMath.Pi()
    phimax = ROOT.TMath.Pi()

    eta_bins = np.concatenate((np.linspace(-2.4, -2.0, 3),np.linspace(-1.8,1.8,37), np.linspace(2.0, 2.4, 3)))
    pt_bins = np.concatenate((np.linspace(25, 35, 6),np.linspace(36,55,20), np.linspace(57, 63, 4), np.linspace(66,70,2)))
    mass_bins = np.concatenate((np.linspace(70, 76, 3),np.linspace(78,86,5), np.linspace(87, 96, 10), np.linspace(98,104,5), np.linspace(106, 110, 2)))
    u_bins = np.concatenate((np.linspace(0, 20, 11),np.linspace(24,80,15), np.linspace(85, 109, 4), np.linspace(120,150,3)))

    # z pt befor and after pt reweighting
    sampMan.cacheDraw("zpt",   "histo_zjets_zpt_WoZptWeight", u_bins, DrawConfig(xmin=0, xmax=150, xlabel='p_{T}^{#mu#mu} [GeV]'), weightname="weight_WoVpt")
    sampMan.cacheDraw("zmass", "histo_zjets_zmass", mass_bins, DrawConfig(xmin=70, xmax=110, xlabel='m_{#mu#mu} [GeV]'))
    sampMan.cacheDraw("zeta",  "histo_zjets_zeta", eta_bins, DrawConfig(xmin=-2.5, xmax=2.5, xlabel='#eta_{#mu#mu} [GeV]', ymax=1e7, ylabel='Events / 1'))
    sampMan.cacheDraw("zy",    "histo_zjets_zrapidity", eta_bins, DrawConfig(xmin=-2.5, xmax=2.5, xlabel='y_{#mu#mu} [GeV]', ymax=1e7, ylabel='Events / 1'))

    sampMan.cacheDraw("nPV", "histo_nPV", 10, 0, 10, DrawConfig(xmin=0, xmax=10, xlabel='# PV', ylabel='Events / 1'))

    sampMan.cacheDraw("LeadMuon_pt", "histo_leadMuon_pt", pt_bins, DrawConfig(xmin=20, xmax=70, xlabel='p_{T}(Leading #mu) [GeV]'))
    sampMan.cacheDraw("LeadMuon_eta", "histo_leadMuon_eta", eta_bins, DrawConfig(xmin=-2.6, xmax=2.6, xlabel='#eta (Leading #mu) [GeV]', ymax=1e7, ylabel='Events / 1'))
    sampMan.cacheDraw("SubleadMuon_pt", "histo_subleadMuon_pt", pt_bins, DrawConfig(xmin=20, xmax=70, xlabel='p_{T}(Subleading #mu) [GeV]'))
    sampMan.cacheDraw("SubleadMuon_eta", "histo_subleadMuon_eta", eta_bins, DrawConfig(xmin=-2.6, xmax=2.6, xlabel='#eta (Subleading #mu) [GeV]', ymax=1e7, ylabel='Events / 1'))

    sampMan.cacheDraw("met", "histo_zjets_pfmet_pt", met_pt_bins, DrawConfig(xmin=0, xmax=100, xlabel='PF MET [GeV]'))
    sampMan.cacheDraw("metPhi", "histo_zjets_pfmet_phi", 30, phimin, phimax, DrawConfig(xmin=phimin, xmax=phimax, xlabel='PF MET #phi'))

    sampMan.launchDraw()

    sampMan.dumpCounts()

    print "Program end..."

    raw_input()
    
    return 

if __name__ == "__main__":
   main()

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

# boolean flag to set either muon or electron channel
# doMuon = False means the electron channel
doMuon = True


def main():
    print "Program start..."

    #ROOT.gROOT.ProcessLine('TFile* f_zpt = TFile::Open("results/zpt_weight.root")')
    #ROOT.gROOT.ProcessLine('TH1D* h_zpt_ratio  = (TH1D*)f_zpt->Get("h_zpt_ratio")')

    if doMuon:
        input_data    = "inputs/wmunu/input_data.txt"
        input_wl0     = "inputs/wmunu/input_wm0.txt"
        input_wl1     = "inputs/wmunu/input_wm1.txt"
        input_wl2     = "inputs/wmunu/input_wm2.txt"
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

    #input_antiiso_data    = "inputs/awmunu/input_data.txt"
    #input_antiiso_wm0     = "inputs/awmunu/input_wm0.txt"
    #input_antiiso_wm1     = "inputs/awmunu/input_wm1.txt"
    #input_antiiso_wm2     = "inputs/awmunu/input_wm2.txt"
    #input_antiiso_ttbar   = "inputs/awmunu/input_ttbar_dilepton.txt"
    #input_antiiso_ttbar_1lep = "inputs/awmunu/input_ttbar_singlelepton.txt"
    #input_antiiso_ttbar_0lep = "inputs/awmunu/input_ttbar_hadronic.txt"
    #input_antiiso_ww      = "inputs/awmunu/input_ww.txt"
    #input_antiiso_wz      = "inputs/awmunu/input_wz.txt"
    #input_antiiso_zz      = "inputs/awmunu/input_zz.txt"
    #input_antiiso_zxx     = "inputs/awmunu/input_zxx.txt"
    #input_antiiso_wx0     = "inputs/awmunu/input_wx0.txt"
    #input_antiiso_wx1     = "inputs/awmunu/input_wx1.txt"
    #input_antiiso_wx2     = "inputs/awmunu/input_wx2.txt"

    DataSamp  = Sample(input_data, isMC=False, legend="Data", name="Data", isWSR=True)
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
    global sampMan
    sampMan = SampleManager(DataSamp, [Wl0Samp, Wl1Samp, Wl2Samp, TTbarSamp, TT1LepSamp, TT0LepSamp, WWSamp, WZSamp, ZZSamp, ZXXSamp, Wx0Samp, Wx1Samp, Wx2Samp])
    sampMan.groupMCs(["wx0", "wx1", "wx2"], 'wx', 216, 'wx')
    sampMan.groupMCs(["WW", "WZ", "ZZ"], 'vv', 216, 'vv')
    sampMan.groupMCs(["ttbar_dilepton", "ttbar_1lepton", "ttbar_0lepton"], "ttbar", 96, "t#bar{t}")
    label = "W#rightarrow#mu#nu" if doMuon else "W#rightarrow e#nu"
    sampMan.groupMCs(['wl0', 'wl1', 'wl2'], "wlnu", 92,"W#rightarrow#mu#nu")
    #sampMan.groupMCs(['Data_aiso', 'wm0_aiso', 'wm1_aiso', 'wm2_aiso', 'ttbar_dilepton_aiso', 'ttbar_1lepton_aiso', 'ttbar_0lepton_aiso', 'WW_aiso', 'WZ_aiso', 'ZZ_aiso', 'ZXX_aiso', 'wx0_aiso', 'wx1_aiso' ,'wx2_aiso'], 'QCD', 226, 'QCD')

    sampMan.DefineAll("Lep_pt",  "lep.Pt()")
    sampMan.DefineAll("Lep_eta", "lep.Eta()")
    sampMan.DefineAll("Lep_phi", "lep.Phi()")
    for i in [1, 2, 3, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]:
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

    # WpT bins
    # bin0 is the inclusive bin
    sampMan.DefineAll("WpT_bin0", "WpT>=0.")
    sampMan.DefineAll("WpT_bin1", "WpT<8.0")
    sampMan.DefineAll("WpT_bin2", "WpT>8.0 &&  WpT<16.0")
    sampMan.DefineAll("WpT_bin3", "WpT>16.0 && WpT<24.0")
    sampMan.DefineAll("WpT_bin4", "WpT>24.0 && WpT<32.0")
    sampMan.DefineAll("WpT_bin5", "WpT>32.0 && WpT<40.0")
    sampMan.DefineAll("WpT_bin6", "WpT>40.0 && WpT<50.0")
    sampMan.DefineAll("WpT_bin7", "WpT>50.0 && WpT<70.0")
    sampMan.DefineAll("WpT_bin8", "WpT>70.0 && WpT<100.0")

    # eta bins for electrons: barral and endcap
    sampMan.DefineAll("lepEta_bin0", "1.0")
    sampMan.DefineAll("lepEta_bin1", "abs(lep.Eta()) <= 1.4442") 
    sampMan.DefineAll("lepEta_bin2", "abs(lep.Eta()) > 1.4442")

    if doMuon:
        etabins = ["lepEta_bin0"]
    else:
        etabins = ["lepEta_bin0", "lepEta_bin1", "lepEta_bin2"]
    
    #sampMan.DefineAll("weight_wplus", "weight_WoVpt * (q>0)")
    #sampMan.DefineAll("weight_wminus", "weight_WoVpt * (q<0)")
    #for wpt in ["WpT_bin0", "WpT_bin1", "WpT_bin2", "WpT_bin3", "WpT_bin4", "WpT_bin5", "WpT_bin6", "WpT_bin7", "WpT_bin8"]:
    #    sampMan.DefineAll("weight_wplus_{}".format(wpt),  "weight_wplus * {}".format(wpt))
    #    sampMan.DefineAll("weight_wminus_{}".format(wpt), "weight_wminus * {}".format(wpt))

    # sample weight
    # 0 is the central one with all corrections
    for i in [0, 1, 2, 3, 4, 6, 7]:
        for wpt in ["WpT_bin0", "WpT_bin1", "WpT_bin2", "WpT_bin3", "WpT_bin4", "WpT_bin5", "WpT_bin6", "WpT_bin7", "WpT_bin8"]:
            for lepeta in etabins:
                sampMan.DefineMC("weight_{}_{}_{}".format(str(i), wpt, lepeta), "evtWeight[{}] * mcnorm * {} * {}".format(str(i), wpt, lepeta))
                DataSamp.Define("weight_{}_{}_{}".format(str(i), wpt, lepeta),  "1.0 * {} * {}".format(wpt, lepeta))

                sampMan.DefineAll("weight_wplus_{}_{}_{}".format(str(i),  wpt, lepeta), "weight_{}_{}_{} * (q>0)".format(str(i), wpt, lepeta))
                sampMan.DefineAll("weight_wminus_{}_{}_{}".format(str(i), wpt, lepeta), "weight_{}_{}_{} * (q<0)".format(str(i), wpt, lepeta))

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

    lepname = "mu" if doMuon else "e"

    # nPV
    sampMan.cacheDraw("npv", "histo_nPV_"+lepname+"nu", 10, 0, 10, DrawConfig(xmin=0, xmax=10, xlabel='# PV', ylabel='Events / 1', dology=False, ymax = 9e5), weightname = "weight_0_WpT_bin0_lepEta_bin0")
    # lepton pt
    sampMan.cacheDraw("Lep_pt", "histo_Lep_pt_"+lepname+"nu", pt_bins, DrawConfig(xmin=20, xmax=70, xlabel='p_{T}(#mu) [GeV]', dology=False, ymax=1.05e5), weightname="weight_0_WpT_bin0_lepEta_bin0")
    sampMan.cacheDraw("Lep_eta", "histo_Lep_eta_"+lepname+"nu", eta_bins, DrawConfig(xmin=-2.6, xmax=2.6, xlabel='#eta (#mu) [GeV]', ylabel='Events / 1', dology=False, ymax=6e5), weightname = "weight_0_WpT_bin0_lepEta_bin0")
    sampMan.cacheDraw("Lep_pt", "histo_Lep_pt_"+lepname+"plusnu", pt_bins, DrawConfig(xmin=20, xmax=70, xlabel='p_{T}(#mu) [GeV]', dology=False, ymax=6.05e4), weightname = "weight_wplus_0_WpT_bin0_lepEta_bin0")
    sampMan.cacheDraw("Lep_pt", "histo_Lep_pt_"+lepname+"minusnu", pt_bins, DrawConfig(xmin=20, xmax=70, xlabel='p_{T}(#mu) [GeV]', dology=False, ymax=5.05e4), weightname = "weight_wminus_0_WpT_bin0_lepEta_bin0")

    # met and raw_met
    sampMan.cacheDraw("met_1_pt", "histo_wjets_"+lepname+"nu_pfmet_pt", met_pt_bins, DrawConfig(xmin=0, xmax=80, xlabel='PF MET [GeV]', dology=False, ymax=8e4), weightname="weight_0_WpT_bin0_lepEta_bin0")
    sampMan.cacheDraw("met_1_phi", "histo_wjets_"+lepname+"nu_pfmet_phi", 30, phimin, phimax, DrawConfig(xmin=phimin, xmax=phimax, xlabel='PF MET #phi', dology=False, ymax=6e5), weightname = "weight_0_WpT_bin0_lepEta_bin0")

    # deltaPhi
    sampMan.cacheDraw("deltaPhi", "histo_wjets_"+lepname+"nu_deltaphi",      dphi_bins,  DrawConfig(xmin=0., xmax=2*phimax, xlabel='#Delta #phi(#mu, MET)', dology=False, ymax=2e6), weightname="weight_0_WpT_bin0_lepEta_bin0")

    # mt
    sampMan.cacheDraw("mtCorr", "histo_wjets_"+lepname+"nu_mtcorr",      mt_bins, DrawConfig(xmin=40, xmax=120, xlabel="m_{T} [GeV]", dology=False, ymax=7e4), weightname="weight_0_WpT_bin0_lepEta_bin0")
    sampMan.cacheDraw("mtCorr", "histo_wjets_"+lepname+"plusnu_mtcorr",  mt_bins, DrawConfig(xmin=40, xmax=120, xlabel="m_{T} [GeV]", dology=False, ymax=4e4), weightname = "weight_wplus_0_WpT_bin0_lepEta_bin0")
    sampMan.cacheDraw("mtCorr", "histo_wjets_"+lepname+"minusnu_mtcorr", mt_bins, DrawConfig(xmin=40, xmax=120, xlabel="m_{T} [GeV]", dology=False, ymax=3.5e4), weightname = "weight_wminus_0_WpT_bin0_lepEta_bin0")

    # variation on recoil corections
    for i in [1, 2, 3, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]:
        for wpt in ["WpT_bin0", "WpT_bin1", "WpT_bin2", "WpT_bin3", "WpT_bin4", "WpT_bin5", "WpT_bin6", "WpT_bin7", "WpT_bin8"]:
            for lepeta in etabins:
                sampMan.cacheDraw("mT_{}".format(str(i)), "histo_wjets_{}plusnu_mT_{}_{}_{}".format(lepname, str(i), wpt, lepeta),  70, 0, 140, DrawConfig(xmin=0, xmax=140, xlabel="m_{T} [GeV]", dology=False, ymax=4e4, donormalizebin=False, addOverflow=False, addUnderflow=False), weightname = "weight_wplus_0_{}_{}".format(wpt, lepeta))
                sampMan.cacheDraw("mT_{}".format(str(i)), "histo_wjets_{}minusnu_mT_{}_{}_{}".format(lepname, str(i), wpt, lepeta), 70, 0, 140, DrawConfig(xmin=0, xmax=140, xlabel="m_{T} [GeV]", dology=False, ymax=3.5e4, donormalizebin=False, addOverflow=False, addUnderflow=False), weightname = "weight_wminus_0_{}_{}".format(wpt, lepeta))

    # variation on the scale factors
    for i in [0, 1, 2, 3, 4, 6, 7]:
        for wpt in ["WpT_bin0", "WpT_bin1", "WpT_bin2", "WpT_bin3", "WpT_bin4", "WpT_bin5", "WpT_bin6", "WpT_bin7", "WpT_bin8"]:
            for lepeta in etabins:
                sampMan.cacheDraw("mT_1", "histo_wjets_{}plusnu_mtcorr_weight_{}_{}_{}".format(lepname, str(i), wpt, lepeta),  70, 0, 140, DrawConfig(xmin=0, xmax=140, xlabel="m_{T} [GeV]", dology=False, ymax=4e4, donormalizebin=False, addOverflow=False, addUnderflow=False), weightname = "weight_wplus_{}_{}_{}".format(str(i), wpt, lepeta))
                sampMan.cacheDraw("mT_1", "histo_wjets_{}minusnu_mtcorr_weight_{}_{}_{}".format(lepname, str(i), wpt, lepeta), 70, 0, 140, DrawConfig(xmin=0, xmax=140, xlabel="m_{T} [GeV]", dology=False, ymax=3.5e4, donormalizebin=False, addOverflow=False, addUnderflow=False), weightname = "weight_wminus_{}_{}_{}".format(str(i), wpt, lepeta))


    sampMan.launchDraw()
    sampMan.dumpCounts()


    outfile = ROOT.TFile("output_shapes_"+lepname+"nu.root", "recreate")

    for chg in [lepname+"plusnu", lepname+"minusnu"]:
        for wpt in ["WpT_bin0", "WpT_bin1", "WpT_bin2", "WpT_bin3", "WpT_bin4", "WpT_bin5", "WpT_bin6", "WpT_bin7", "WpT_bin8"]:
            for lepeta in etabins:
                hdata = sampMan.hdatas["histo_wjets_{}_mT_1_{}_{}".format(chg, wpt, lepeta)]
                hdata.SetDirectory(outfile)
                hdata.Write()

                hsmcs_central = sampMan.hsmcs["histo_wjets_{}_mT_1_{}_{}".format(chg, wpt, lepeta)]
                hlists_central = list(hsmcs_central.GetHists())

                # recoil systematics
                for i in [2, 3, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]:
                    hsmcs_up = sampMan.hsmcs["histo_wjets_{}_mT_{}_{}_{}".format(chg, str(i), wpt, lepeta)]

                    hlists_up = list(hsmcs_up.GetHists())
                    for ih in xrange(len(hlists_up)):
                        # loop over different simulated processes
                        hcen = hlists_central[ih]
                        hup  = hlists_up[ih]
                        hup.SetName("{}_SysRecoil{}Up".format(hcen.GetName(), str(i)))

                        hdn  = hcen.Clone("{}_SysRecoil{}Down".format(hcen.GetName(), str(i)))
                        for ibin in xrange(1, hup.GetNbinsX()+1):
                            hdn.SetBinContent(ibin, 2*hcen.GetBinContent(ibin) - hup.GetBinContent(ibin))

                        hcen.SetDirectory(outfile)
                        hcen.Write()
                        hup.SetDirectory(outfile)
                        hup.Write()
                        hdn.SetDirectory(outfile)
                        hdn.Write()

                # weights/corrections
                for i in [0, 1, 2, 3, 4, 6, 7]:
                    hsmcs_up = sampMan.hsmcs["histo_wjets_{}_mtcorr_weight_{}_{}_{}".format(chg, str(i), wpt, lepeta)]
                    hlists_up = list(hsmcs_up.GetHists())
                    for ih in xrange(len(hlists_up)):
                        # loop over different processes
                        hcen = hlists_central[ih]
                        hup  = hlists_up[ih]
                        if i<=6:
                            hup.SetName("{}_SysWeight{}Up".format(hcen.GetName(), str(i)))
                        else:
                            hup.SetName("{}_SysWeight6Down".format(hcen.GetName()))

                        if i<6:
                            # for prefire, the up and down weights are calculated seperately
                            hdn  = hcen.Clone("{}_SysWeight{}Down".format(hcen.GetName(), str(i)))
                            for ibin in xrange(1, hup.GetNbinsX()+1):
                                hdn.SetBinContent(ibin, 2*hcen.GetBinContent(ibin) - hup.GetBinContent(ibin))
                        
                        hup.SetDirectory(outfile)
                        hup.Write()
                        if i<6:
                            hdn.SetDirectory(outfile)
                            hdn.Write()

    outfile.Close()

    print "Program end..."

    raw_input()
    
    return 

if __name__ == "__main__":
   main()

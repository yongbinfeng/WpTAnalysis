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

# boolean flag to set either muon or electron channel
# doMuon = False means the electron channel
doMuon = True

def main():
    print "Program start..."

    #ROOT.gROOT.ProcessLine('TFile* f_zpt = TFile::Open("results/zpt_weight.root")')
    #ROOT.gROOT.ProcessLine('TH1D* h_zpt_ratio  = (TH1D*)f_zpt->Get("h_zpt_ratio")')

    if doMuon:
        input_antiiso_data    = "inputs/awmunu/input_data.txt"
        input_antiiso_wl0     = "inputs/awmunu/input_wm0.txt"
        input_antiiso_wl1     = "inputs/awmunu/input_wm1.txt"
        input_antiiso_wl2     = "inputs/awmunu/input_wm2.txt"
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
    else:
        input_antiiso_data    = "inputs/awenu/input_data.txt"
        input_antiiso_wl0     = "inputs/awenu/input_we0.txt"
        input_antiiso_wl1     = "inputs/awenu/input_we1.txt"
        input_antiiso_wl2     = "inputs/awenu/input_we2.txt"
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

    ## for the QCD background estimation (data-driven)
    qcdnorm  = -1.0
    DataAisoSamp  = Sample(input_antiiso_data, isMC=False, name="Data_aiso", isWSR=True, additionalnorm = -1.0 * qcdnorm, legend = 'QCD', color='226')
    # W -> lnu
    Wl0AisoSamp   = Sample(input_antiiso_wl0, isMC=True, name = "wl0_aiso", isWSR=True, additionalnorm=-1.0 * qcdnorm)
    Wl1AisoSamp   = Sample(input_antiiso_wl1, isMC=True, name = "wl1_aiso", isWSR=True, additionalnorm=-1.0 * qcdnorm)
    Wl2AisoSamp   = Sample(input_antiiso_wl2, isMC=True, name = "wl2_aiso", isWSR=True, additionalnorm=-1.0 * qcdnorm)
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

    sampMan = SampleManager(DataAisoSamp, [Wl0AisoSamp, Wl1AisoSamp, Wl2AisoSamp, TTbarAisoSamp, TT1LepAisoSamp, TT0LepAisoSamp, WWAisoSamp, WZAisoSamp, ZZAisoSamp, ZXXAisoSamp, Wx0AisoSamp, Wx1AisoSamp, Wx2AisoSamp])

    sampMan.groupMCs(["WW_aiso", "WZ_aiso", "ZZ_aiso", "ZXX_aiso", "wx0_aiso", "wx1_aiso", "wx2_aiso"], "EWK", 216, "EWK")
    sampMan.groupMCs(["ttbar_dilepton_aiso", "ttbar_1lepton_aiso", "ttbar_0lepton_aiso"], "ttbar", 96, "t#bar{t}")
    label = "W#rightarrow#mu#nu" if doMuon else "W#rightarrow e#nu"
    sampMan.groupMCs(['wl0_aiso', 'wl1_aiso', 'wl2_aiso'], "wlnu", 92, label)

    sampMan.DefineAll("Lep_pt", "lep.Pt()")
    sampMan.ApplyCutAll("Lep_pt > 25.0")

    sampMan.DefineAll("w_iso1", "(relIso > 0.15 && relIso < 0.30)")
    sampMan.DefineAll("w_iso2", "(relIso > 0.30 && relIso < 0.45)")
    sampMan.DefineAll("w_iso3", "(relIso > 0.45 && relIso < 0.55)")
    sampMan.DefineAll("w_iso4", "(relIso > 0.55 && relIso < 0.65)")

    # muon and electron isolation distributions are different
    # more coarse binning for electrons to make sure enough statistics
    if doMuon:
        sampMan.DefineAll("w_iso5", "(relIso > 0.20 && relIso < 0.25)")
        sampMan.DefineAll("w_iso6", "(relIso > 0.25 && relIso < 0.30)")
        sampMan.DefineAll("w_iso7", "(relIso > 0.30 && relIso < 0.35)")
        sampMan.DefineAll("w_iso8", "(relIso > 0.35 && relIso < 0.40)")
        sampMan.DefineAll("w_iso9", "(relIso > 0.40 && relIso < 0.45)")
        sampMan.DefineAll("w_iso10", "(relIso > 0.45 && relIso < 0.50)")
        sampMan.DefineAll("w_iso11", "(relIso > 0.50 && relIso < 0.55)")
        sampMan.DefineAll("w_iso12", "(relIso > 0.55 && relIso < 0.60)")
        sampMan.DefineAll("w_iso13", "(relIso > 0.60 && relIso < 0.65)")
        isobins = ["iso1", "iso2", "iso3", "iso4", "iso5", "iso6", "iso7", "iso8", "iso9", "iso10", "iso11", "iso12", "iso13"]
    else:
        sampMan.DefineAll("w_iso5", "(relIso > 0.20 && relIso < 0.25)")
        sampMan.DefineAll("w_iso6", "(relIso > 0.25 && relIso < 0.30)")
        sampMan.DefineAll("w_iso7", "(relIso > 0.30 && relIso < 0.35)")
        sampMan.DefineAll("w_iso8", "(relIso > 0.35 && relIso < 0.40)")
        sampMan.DefineAll("w_iso9", "(relIso > 0.40 && relIso < 0.50)")
        sampMan.DefineAll("w_iso10", "(relIso > 0.50 && relIso < 0.70)")
        isobins = ["iso1", "iso2", "iso3", "iso4", "iso5", "iso6", "iso7", "iso8", "iso9", "iso10"]

    sampMan.DefineAll("Lep_eta", "lep.Eta()")
    sampMan.DefineMC("met_pt", "metVars[1]", excludes=['Data_aiso'])
    sampMan.DefineMC("met_phi", "TVector2::Phi_mpi_pi(metVarsPhi[1])", excludes=['Data_aiso'])
    # Data does not run any recoil corrections
    DataAisoSamp.Define("met_pt", "metVars[0]")
    DataAisoSamp.Define("met_phi", "TVector2::Phi_mpi_pi(metVarsPhi[0])")
    sampMan.DefineAll("met_raw_pt", "metVars[0]")
    sampMan.DefineAll("met_raw_phi", "TVector2::Phi_mpi_pi(metVarsPhi[0])")

    #sampMan.DefineAll("weight_wplus", "weight_WoVpt * (q>0)", excludeGroups=['QCD'])
    #sampMan.DefineAll("weight_wminus", "weight_WoVpt * (q<0)", excludeGroups=['QCD'])
    #sampMan.DefineSpecificMCs("weight_wplus", "weight_WoVpt * (q>0)", sampgroupnames=['QCD'])
    #sampMan.DefineSpecificMCs("weight_wminus", "weight_WoVpt * (q<0)", sampgroupnames=['QCD'])
    sampMan.DefineAll("weight_wplus", "weight_WoVpt * (q>0)")
    sampMan.DefineAll("weight_wminus", "weight_WoVpt * (q<0)")

    #sampMan.DefineAll("weight_eta0to0p8",   "fabs(Lep_eta)<0.8")
    #sampMan.DefineAll("weight_eta0p8to1p4", "fabs(Lep_eta)>0.8 && fabs(Lep_eta)<1.4")
    #sampMan.DefineAll("weight_eta1p4to2",   "fabs(Lep_eta)<1.4 && fabs(Lep_eta)<2.0")
    #sampMan.DefineAll("weight_eta2to2p4",   "fabs(Lep_eta)>2.0 && fabs(Lep_eta)<2.4")

    # recoil variables
    sampMan.DefineAll("V2W", "UVec(lep.Pt(), lep.Phi(), met_pt, met_phi)")
    sampMan.DefineAll("WpT", "V2W.Mod()")
    sampMan.DefineAll("Wphi", "TVector2::Phi_mpi_pi(V2W.Phi())")

    # WpT bins
    sampMan.DefineAll("WpT_bin0", "WpT>0.0")
    sampMan.DefineAll("WpT_bin1", "WpT<8.0")
    sampMan.DefineAll("WpT_bin2", "WpT>8.0 &&  WpT<16.0")
    sampMan.DefineAll("WpT_bin3", "WpT>16.0 && WpT<24.0")
    sampMan.DefineAll("WpT_bin4", "WpT>24.0 && WpT<32.0")
    sampMan.DefineAll("WpT_bin5", "WpT>32.0 && WpT<40.0")
    sampMan.DefineAll("WpT_bin6", "WpT>40.0 && WpT<50.0")
    sampMan.DefineAll("WpT_bin7", "WpT>50.0 && WpT<70.0")
    sampMan.DefineAll("WpT_bin8", "WpT>70.0 && WpT<100.0")
    wptbins = ["WpT_bin0", "WpT_bin1", "WpT_bin2", "WpT_bin3", "WpT_bin4", "WpT_bin5", "WpT_bin6", "WpT_bin7", "WpT_bin8"]

    # eta bins for electrons: barral and endcap
    sampMan.DefineAll("lepEta_bin0", "1.0")
    sampMan.DefineAll("lepEta_bin1", "abs(lep.Eta()) <= 1.4442")
    sampMan.DefineAll("lepEta_bin2", "abs(lep.Eta()) > 1.4442")

    if doMuon:
        etabins = ["lepEta_bin0"]
    else:
        etabins = ["lepEta_bin0", "lepEta_bin1", "lepEta_bin2"]

    # pT_ell / mT
    sampMan.DefineAll("ptOvermT", "Lep_pt / mtCorr")

    #for iso in ['iso1', 'iso2', 'iso3', 'iso4', 'iso5', 'iso6', 'iso7', 'iso8', 'iso9', 'iso10', 'iso11', 'iso12', 'iso13']:
    for iso in isobins:
        for chg in ['WoVpt', 'wplus', 'wminus']:
            #for eta in ['0to0p8', '0p8to1p4', '1p4to2', '2to2p4']:
            #for wpt in ["WpT_bin0", "WpT_bin1", "WpT_bin2", "WpT_bin3", "WpT_bin4", "WpT_bin5", "WpT_bin6", "WpT_bin7", "WpT_bin8"]:
            for wpt in wptbins:
                for lepeta in etabins:
                    sampMan.DefineAll("weight_{}_{}_{}_{}".format(chg, iso, wpt, lepeta), "w_{} * weight_{} * {} * {}".format(iso, chg, wpt, lepeta))

    met_pt_bins = np.array([0., 2.0, 4., 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 33, 36, 39, 42, 45, 48, 51, 55, 60, 65, 70, 75, 80, 90, 100, 110, 120, 135, 150, 165, 180, 200])
    u1_bins = np.array([-40.,-36.,-32., -28., -25., -22.0, -20.0, -18, -16, -14, -12, -10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 53, 56, 59, 64, 68, 72, 76, 80, 85, 90, 100])
    u2_bins = np.array([-80., -70., -65., -60., -56., -52, -48, -44, -40, -37, -34, -31, -28, -25., -22., -20, -18, -16, -14, -12, -10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10, 12, 16, 18, 20, 22, 25, 28, 31, 34, 37, 40, 44, 48, 52, 56, 60, 65, 70, 80])
    #u_bins = np.array([0., 2., 4., 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 43, 46, 49, 52, 56, 60, 64, 68, 72, 76, 80, 85, 90, 95, 100, 105, 110, 115, 120, 130, 140, 150])
    u_bins = np.concatenate((np.linspace(0, 20, 11),np.linspace(24,60,10)))
    phimin = -ROOT.TMath.Pi()
    phimax = ROOT.TMath.Pi()

    eta_bins = np.concatenate((np.linspace(-2.4, -2.0, 3),np.linspace(-1.8,1.8,37), np.linspace(2.0, 2.4, 3)))
    pt_bins = np.concatenate((np.linspace(0,15, 6), np.linspace(17,23,4), np.linspace(25, 35, 6),np.linspace(36,55,20), np.linspace(57, 63, 4), np.linspace(66,70,2),np.linspace(75,100, 6)))
    mass_bins = np.concatenate((np.linspace(0,52,14), np.linspace(55,67, 5), np.linspace(70, 76, 3),np.linspace(78,86,5), np.linspace(87, 96, 10), np.linspace(98,104,5), np.linspace(106, 110, 2)))
    u_bins = np.concatenate((np.linspace(0, 20, 11),np.linspace(24,80,15), np.linspace(85, 109, 4), np.linspace(120,150,3)))
    mt_bins = np.concatenate((np.linspace(0,20,6),np.linspace(22,38,9),np.linspace(40, 60, 21), np.linspace(62, 80, 10), np.linspace(84, 120, 10)))
    met_pt_bins = np.array([0., 2.0, 4., 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 33, 36, 39, 42, 45, 48, 51, 55, 60, 65, 70, 75, 80])

    #for iso in ['iso1', 'iso2', 'iso3', 'iso4', 'iso5', 'iso6', 'iso7', 'iso8', 'iso9']:
    #for iso in ['iso1', 'iso2', 'iso3', 'iso4', 'iso5', 'iso6', 'iso7', 'iso8', 'iso9', 'iso10', 'iso11', 'iso12', 'iso13']:
    for iso in isobins:
        for chg in ['WoVpt', 'wplus', 'wminus']:
            #for wpt in ["WpT_bin0", "WpT_bin1", "WpT_bin2", "WpT_bin3", "WpT_bin4", "WpT_bin5", "WpT_bin6", "WpT_bin7", "WpT_bin8"]:
            for wpt in wptbins:
                for lepeta in etabins:
                    if chg == 'WoVpt':
                        continue
                    strname = "weight_{}_{}_{}_{}".format(chg, iso, wpt, lepeta)
                    #outputname = "histo_wjetsAntiIso_mtcorr_" + strname
                    ##sampMan.cacheDraw("mtCorr", outputname, mt_bins, DrawConfig(xmin=0, xmax=120, xlabel="m_{T} [GeV]", dology=False, ymax=1.0e4), weightname = strname)
                    #sampMan.cacheDraw("mtCorr", outputname, 25, 40, 140, DrawConfig(xmin=40, xmax=140, xlabel="m_{T} [GeV]", dology=False, ymax=1e4, donormalizebin=False, addOverflow=False, addUnderflow=False, showratio=True), weightname = strname)    

                    outputname = "histo_wjetsAntiIso_fullrange_mtcorr_" + strname 
                    sampMan.cacheDraw("mtCorr", outputname, 70,  0, 140, DrawConfig(xmin=0, xmax=140, xlabel="m_{T} [GeV]", dology=False, ymax=2e4, donormalizebin=False, addOverflow=False, addUnderflow=False, showratio=True), weightname = strname)
                    #outputname = "histo_wjetsAntiIso_WpTcorr_" + strname
                    #sampMan.cacheDraw("WpT", outputname, u_bins, DrawConfig(xmin=0, xmax=60, xlabel="p^{W}_{T} [GeV]", dology=False, ymax=1.0e4), weightname = strname)
                    #outputname = "histo_wjetsAntiIso_leppt_" + strname
                    #sampMan.cacheDraw("Lep_pt", outputname, pt_bins, DrawConfig(xmin=25, xmax=100, xlabel="p^{#mu}_{T} [GeV]", dology=False, ymax=1.0e4), weightname = strname)
                    #outputname = "histo_wjetsAntiIso_metpt_" + strname
                    #sampMan.cacheDraw("met_pt", outputname, met_pt_bins, DrawConfig(xmin=0, xmax=80, xlabel="PF MET [GeV]", dology=False, ymax=1e4), weightname = strname)

    sampMan.launchDraw()


    #hmts_comp = OrderedDict()
    ##for iso in ['iso1', 'iso2', 'iso3', 'iso4', 'iso5', 'iso6', 'iso7', 'iso8', 'iso9']:
    #for iso in ['iso1', 'iso2', 'iso3', 'iso4', 'iso5', 'iso6', 'iso7', 'iso8', 'iso9', 'iso10', 'iso11', 'iso12', 'iso13']:
    #    for chg in ['WoVpt', 'wplus', 'wminus']:
    #        for wpt in ["WpT_bin0", "WpT_bin1", "WpT_bin2", "WpT_bin3", "WpT_bin4", "WpT_bin5", "WpT_bin6", "WpT_bin7", "WpT_bin8"]:
    #            for lepeta in etabins:
    #                if chg == "WoVpt":
    #                    continue
    #                strname = "weight_{}_{}_{}_{}".format(chg, iso, wpt, lepeta)

    #                outputname = "histo_wjetsAntiIso_mtcorr_" + strname
    #                hstacked = THStack2TH1(sampMan.hsmcs[outputname])
    #                for ibin in xrange(hstacked.GetNbinsX()+1):
    #                    # hstacked should always be above 0
    #                    hstacked.SetBinContent(ibin, max(hstacked.GetBinContent(ibin), 0))
    #                sampMan.hdatas[outputname].Add( hstacked, -1.0 )
    #                hmts_comp[strname] = sampMan.hdatas[outputname]
    #                hmts_comp[strname].SetName(outputname)

    #outfile = ROOT.TFile("output_qcdshape.root", "recreate")

    #for chg in ["wplus", "wminus"]:
    #    for i in xrange(1,13):
    #        for wpt in ["WpT_bin0", "WpT_bin1", "WpT_bin2", "WpT_bin3", "WpT_bin4", "WpT_bin5", "WpT_bin6", "WpT_bin7", "WpT_bin8"]:
    #            iso = "iso" + str(i)
    #            iso_next = "iso" + str(i+1)
    #            strname = "weight_{}_{}_{}".format(chg, iso, wpt)
    #            strname_next = "weight_{}_{}_{}".format(chg, iso_next, wpt)

    #            outputname = "histo_wjetsAntiIso_mtcorr_" + strname

    #            hcenter = hmts_comp[strname]
    #            hup = hmts_comp[strname_next].Clone(outputname+"_shapeUp")
    #            hdown = hmts_comp[strname_next].Clone(outputname+"_shapeDown")

    #            hup.Scale(hcenter.Integral() / hup.Integral())

    #            for ibin in xrange(1, hcenter.GetNbinsX()+1):
    #                center = hcenter.GetBinContent(ibin)
    #                up = hup.GetBinContent(ibin)
    #                hdown.SetBinContent(ibin, max(2*center - up, 0))

    #                hcenter.SetBinContent(ibin, max(center, 0))
    #                hup.SetBinContent(ibin, max(up, 0))

    #            hcenter.SetDirectory(outfile)
    #            hcenter.Write()
    #            hup.SetDirectory(outfile)
    #            hup.Write()
    #            hdown.SetDirectory(outfile)
    #            hdown.Write()

    #outfile.Close()


    # repeat the same for the full mT one
    hmts_comp = OrderedDict()
    #for iso in ['iso1', 'iso2', 'iso3', 'iso4', 'iso5', 'iso6', 'iso7', 'iso8', 'iso9']:
    #for iso in ['iso1', 'iso2', 'iso3', 'iso4', 'iso5', 'iso6', 'iso7', 'iso8', 'iso9', 'iso10', 'iso11', 'iso12', 'iso13']:
    for iso in isobins:
        for chg in ['WoVpt', 'wplus', 'wminus']:
            #for wpt in ["WpT_bin0", "WpT_bin1", "WpT_bin2", "WpT_bin3", "WpT_bin4", "WpT_bin5", "WpT_bin6", "WpT_bin7", "WpT_bin8"]:
            for wpt in wptbins:
                for lepeta in etabins:
                    if chg == "WoVpt":
                        continue
                    strname = "weight_{}_{}_{}_{}".format(chg, iso, wpt, lepeta)

                    outputname = "histo_wjetsAntiIso_fullrange_mtcorr_" + strname
                    hstacked = THStack2TH1(sampMan.hsmcs[outputname])
                    for ibin in xrange(hstacked.GetNbinsX()+1):
                        # hstacked should always be above 0
                        hstacked.SetBinContent(ibin, max(hstacked.GetBinContent(ibin), 0))
                    sampMan.hdatas[outputname].Add( hstacked, -1.0 )
                    hmts_comp[strname] = sampMan.hdatas[outputname]
                    hmts_comp[strname].SetName(outputname)

    lepname = "mu" if doMuon else "e"
    outfile = ROOT.TFile("output_qcdshape_fullrange_"+lepname+"nu.root", "recreate")

    for chg in ["wplus", "wminus"]:
        #for i in xrange(1,13):
        #for i in xrange(1,10):
        for i in xrange(1, len(isobins)):
            #for wpt in ["WpT_bin0", "WpT_bin1", "WpT_bin2", "WpT_bin3", "WpT_bin4", "WpT_bin5", "WpT_bin6", "WpT_bin7", "WpT_bin8"]:
            for wpt in wptbins:
                for lepeta in etabins:
                    # skip iso10 as it is used for the uncertainty cal for iso9
                    iso = "iso" + str(i)
                    iso_next = "iso" + str(i+1)
                    strname = "weight_{}_{}_{}_{}".format(chg, iso, wpt, lepeta)
                    strname_next = "weight_{}_{}_{}_{}".format(chg, iso_next, wpt, lepeta)

                    outputname = "histo_wjetsAntiIso_fullrange_mtcorr_" + strname

                    hcenter = hmts_comp[strname]
                    hup = hmts_comp[strname_next].Clone(outputname+"_shapeUp")
                    hdown = hmts_comp[strname_next].Clone(outputname+"_shapeDown")

                    hup.Scale(hcenter.Integral() / hup.Integral())

                    for ibin in xrange(1, hcenter.GetNbinsX()+1):
                        center = hcenter.GetBinContent(ibin)
                        up = hup.GetBinContent(ibin)
                        hdown.SetBinContent(ibin, max(2*center - up, 0))

                        hcenter.SetBinContent(ibin, max(center, 0))
                        hup.SetBinContent(ibin, max(up, 0))

                    hcenter.SetDirectory(outfile)
                    hcenter.Write()
                    hup.SetDirectory(outfile)
                    hup.Write()
                    hdown.SetDirectory(outfile)
                    hdown.Write()

    outfile.Close()



    # calculate the shape uncertainty and save them to root file
    #hmts_comp = OrderedDict()
    ##hwpts_comp = OrderedDict()
    ##hmupts_comp = OrderedDict()
    ##hmetpts_comp = OrderedDict()
    #for iso in ['iso1', 'iso2', 'iso3', 'iso4', 'iso5', 'iso6', 'iso7', 'iso8', 'iso9']:
    #    for chg in ['WoVpt', 'wplus', 'wminus']:
    #        if chg == "WoVpt":
    #            continue
    #        strname = "weight_{}_{}".format(chg, iso)

    #        outputname = "histo_wjetsAntiIso_mtcorr_" + strname
    #        sampMan.hdatas[outputname].Add( THStack2TH1(sampMan.hsmcs[outputname]), -1 )
    #        hmts_comp[strname] = sampMan.hdatas[outputname]

            #outputname = "histo_wjetsAntiIso_WpTcorr_" + strname
            #sampMan.hdatas[outputname].Add( THStack2TH1(sampMan.hsmcs[outputname]), -1 )
            #hwpts_comp[strname] = sampMan.hdatas[outputname]

            #outputname = "histo_wjetsAntiIso_leppt_" + strname
            #sampMan.hdatas[outputname].Add( THStack2TH1(sampMan.hsmcs[outputname]), -1 )
            #hmupts_comp[strname] = sampMan.hdatas[outputname]

            #outputname = "histo_wjetsAntiIso_metpt_" + strname
            #sampMan.hdatas[outputname].Add( THStack2TH1(sampMan.hsmcs[outputname]), -1 )
            #hmetpts_comp[strname] = sampMan.hdatas[outputname]

            #labels = ["|#eta|<0.8", "0.8<|#eta|<1.4", "1.4<|#eta|<2.0", "2.0<|#eta|<2.4"]
            #DrawHistos(hmts_comp.values(), labels, 0, 120, 'm_{T} [GeV]', 2e-4, 5e-1, 'A.U.', 'histo_wjetsAntiIso_comp_mT', dology=True, showratio=True, donormalize=True, mycolors=[1,2,3, 4], ratiobase=1, drawoptions=["e", "e", "e", "e"], legendoptions=['LPE', 'LPE', 'LPE', 'LPE'], yrmin=0.71, yrmax=1.29)

            #DrawHistos(hwpts_comp.values(), labels, 0, 60, 'p^{W}_{T} [GeV]', 2e-3, 5e-1, 'A.U.', 'histo_wjetsAntiIso_comp_WpT', dology=True, showratio=True, donormalize=True, mycolors=[1,2,3, 4], ratiobase=1, drawoptions=["e", "e", "e", "e"], legendoptions=['LPE', 'LPE', 'LPE', 'LPE'], yrmin=0.71, yrmax=1.29)

            #DrawHistos(hmupts_comp.values(), labels, 25, 100, 'p^{#mu}_{T} [GeV]', 2e-3, 5e-1, 'A.U.', 'histo_wjetsAntiIso_comp_mupt', dology=True, showratio=True, donormalize=True, mycolors=[1,2,3, 4], ratiobase=1, drawoptions=["e", "e", "e", "e"], legendoptions=['LPE', 'LPE', 'LPE', 'LPE'], yrmin=0.71, yrmax=1.29)

            #DrawHistos(hmetpts_comp.values(), labels, 0, 80, 'p^{miss}_{T} [GeV]', 2e-3, 5e-1, 'A.U.', 'histo_wjetsAntiIso_comp_metpt', dology=True, showratio=True, donormalize=True, mycolors=[1,2,3, 4], ratiobase=1, drawoptions=["e", "e", "e", "e"], legendoptions=['LPE', 'LPE', 'LPE', 'LPE'], yrmin=0.71, yrmax=1.29)



    #DrawHistos(hmts_comp.values(), ['0.15<I<0.30', '0.30<I<0.45', '0.45<I<0.55'], 40, 120, 'm_{T} [GeV]', 2e-4, 5e-1, 'A.U.', 'histo_wjetsAntiIso_comp_mT', dology=True, showratio=True, donormalize=True, mycolors=[1,2,3], ratiobase=1, drawoptions=["e", "e", "e"], legendoptions=['LPE', 'LPE', 'LPE'], yrmin=0.71, yrmax=1.29)

    #DrawHistos(hwpts_comp.values(), ['0.15<I<0.30', '0.30<I<0.45', '0.45<I<0.55'], 0, 60, 'p^{W}_{T} [GeV]', 2e-3, 5e-1, 'A.U.', 'histo_wjetsAntiIso_comp_WpT', dology=True, showratio=True, donormalize=True, mycolors=[1,2,3], ratiobase=1, drawoptions=["e", "e", "e"], legendoptions=['LPE', 'LPE', 'LPE'], yrmin=0.71, yrmax=1.29)

    #DrawHistos(hmupts_comp.values(), ['0.15<I<0.30', '0.30<I<0.45', '0.45<I<0.55'], 25, 70, 'p^{#mu}_{T} [GeV]', 2e-3, 5e-1, 'A.U.', 'histo_wjetsAntiIso_comp_mupt', dology=True, showratio=True, donormalize=True, mycolors=[1,2,3], ratiobase=1, drawoptions=["e", "e", "e"], legendoptions=['LPE', 'LPE', 'LPE'], yrmin=0.71, yrmax=1.29)

    #DrawHistos(hmetpts_comp.values(), ['0.15<I<0.30', '0.30<I<0.45', '0.45<I<0.55'], 0, 80, 'p^{W}_{T} [GeV]', 2e-3, 5e-1, 'A.U.', 'histo_wjetsAntiIso_comp_metpt', dology=True, showratio=True, donormalize=True, mycolors=[1,2,3], ratiobase=1, drawoptions=["e", "e", "e"], legendoptions=['LPE', 'LPE', 'LPE'], yrmin=0.71, yrmax=1.29)


    #sampMan.launchDraw()

    sampMan.dumpCounts()

    print "Program end..."

    raw_input()
    
    return 

if __name__ == "__main__":
   main()

"""
test on the anti-Isolated region mixed with signal region
"""
import ROOT
import numpy as np
from collections import OrderedDict
import sys
import argparse
from CMSPLOTS.myFunction import THStack2TH1

from modules.SampleManager import DrawConfig, Sample, SampleManager

ROOT.gROOT.SetBatch(True)
ROOT.ROOT.EnableImplicitMT(10)

def main():
    print("Program start...")

    parser = argparse.ArgumentParser(description="Make plots for Wlnu analysis in the lepton anti-isolated region")
    parser.add_argument("--doTest", action="store_true", dest="doTest", help="Run on a subset of samples for debugging; false runs on all dataset.")
    parser.add_argument("--doWpT", action="store_true", dest="doWpT", help="Bin in different W pt bins; false runs inclusively")
    parser.add_argument("--do5TeV", action="store_true", dest="do5TeV", help="Analyze the 5TeV data; false runs on 13TeV data")
    parser.add_argument("--doElectron", action="store_true", dest="doElectron", help="Analyze the electron channel; false runs the muon channel")
    parser.add_argument("--applyScaling", action="store_true", dest="applyScaling", help="Scale the MC cross section by 30 percents and subtract")
    args = parser.parse_args()

    doMuon = not args.doElectron
    applyScaling = args.applyScaling
    doTest = args.doTest
    doWpT  = args.doWpT
    do5TeV = args.do5TeV

    print("doMuon: ", doMuon)
    print("applyScaling: ", applyScaling)
    print("doTest:", doTest)
    print("doWpT:", doWpT)
    print("do5TeV:", do5TeV)

    if not do5TeV:
        if doTest:
            input_antiiso_data    = "inputs/awmunu_test/input_data.txt"
            input_antiiso_wl0     = "inputs/awmunu_test/input_wm0.txt"
            input_antiiso_wl1     = "inputs/awmunu_test/input_wm1.txt"
            input_antiiso_wl2     = "inputs/awmunu_test/input_wm2.txt"
            input_antiiso_ttbar   = "inputs/awmunu_test/input_ttbar_dilepton.txt"
            input_antiiso_ttbar_1lep = "inputs/awmunu_test/input_ttbar_singlelepton.txt"
            input_antiiso_ttbar_0lep = "inputs/awmunu_test/input_ttbar_hadronic.txt"
            input_antiiso_ww      = "inputs/awmunu_test/input_ww.txt"
            input_antiiso_wz      = "inputs/awmunu_test/input_wz.txt"
            input_antiiso_zz      = "inputs/awmunu_test/input_zz.txt"
            input_antiiso_zxx     = "inputs/awmunu_test/input_zxx.txt"
            input_antiiso_wx0     = "inputs/awmunu_test/input_wx0.txt"
            input_antiiso_wx1     = "inputs/awmunu_test/input_wx1.txt"
            input_antiiso_wx2     = "inputs/awmunu_test/input_wx2.txt"
            
        elif doMuon:
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
    else:
        # analyze 5TeV data
        if doMuon:
            input_antiiso_data    = "inputs_5TeV/awmunu/input_data.txt"
            input_antiiso_wl      = "inputs_5TeV/awmunu/input_wm.txt"
            input_antiiso_ttbar   = "inputs_5TeV/awmunu/input_ttbar.txt"
            input_antiiso_ww      = "inputs_5TeV/awmunu/input_ww.txt"
            input_antiiso_wz      = "inputs_5TeV/awmunu/input_wz.txt"
            input_antiiso_zz2l    = "inputs_5TeV/awmunu/input_zz2l.txt"
            input_antiiso_zz4l    = "inputs_5TeV/awmunu/input_zz4l.txt"
            input_antiiso_zxx     = "inputs_5TeV/awmunu/input_zxx.txt"
            input_antiiso_wx      = "inputs_5TeV/awmunu/input_wx.txt"
        else:   
            input_antiiso_data    = "inputs_5TeV/awenu/input_data.txt"
            input_antiiso_wl      = "inputs_5TeV/awenu/input_we.txt"
            input_antiiso_ttbar   = "inputs_5TeV/awenu/input_ttbar.txt"
            input_antiiso_ww      = "inputs_5TeV/awenu/input_ww.txt"
            input_antiiso_wz      = "inputs_5TeV/awenu/input_wz.txt"
            input_antiiso_zz2l    = "inputs_5TeV/awenu/input_zz2l.txt"
            input_antiiso_zz4l    = "inputs_5TeV/awenu/input_zz4l.txt"
            input_antiiso_zxx     = "inputs_5TeV/awenu/input_zxx.txt"
            input_antiiso_wx      = "inputs_5TeV/awenu/input_wx.txt"

    ## for the QCD background estimation (data-driven)
    qcdnorm = 1.0
    mcscale = 1.0
    if applyScaling:
        # scale up the MC for 30%
        mcscale = 1.3
    DataAisoSamp  = Sample(input_antiiso_data, isMC=False, name="Data_aiso", isWSR=True, additionalnorm = qcdnorm, legend = 'QCD', color='226')

    if not do5TeV:
        # W -> lnu
        Wl0AisoSamp   = Sample(input_antiiso_wl0, isMC=True, name = "wl0_aiso", isWSR=True, additionalnorm= qcdnorm * mcscale)
        Wl1AisoSamp   = Sample(input_antiiso_wl1, isMC=True, name = "wl1_aiso", isWSR=True, additionalnorm= qcdnorm * mcscale)
        Wl2AisoSamp   = Sample(input_antiiso_wl2, isMC=True, name = "wl2_aiso", isWSR=True, additionalnorm= qcdnorm * mcscale)
        # ttbar
        TTbarAisoSamp  = Sample(input_antiiso_ttbar, isMC=True, name = "ttbar_dilepton_aiso",     isWSR=True, additionalnorm= qcdnorm * mcscale)
        TT1LepAisoSamp = Sample(input_antiiso_ttbar_1lep, isMC=True, name = "ttbar_1lepton_aiso", isWSR=True, additionalnorm= qcdnorm * mcscale)
        TT0LepAisoSamp = Sample(input_antiiso_ttbar_0lep, isMC=True, name = "ttbar_0lepton_aiso", isWSR=True, additionalnorm= qcdnorm * mcscale)
        ## dibosons
        WWAisoSamp = Sample(input_antiiso_ww, isMC=True, name = "WW_aiso", isWSR=True, additionalnorm= qcdnorm * mcscale)
        WZAisoSamp = Sample(input_antiiso_wz, isMC=True, name = "WZ_aiso", isWSR=True, additionalnorm= qcdnorm * mcscale)
        ZZAisoSamp = Sample(input_antiiso_zz, isMC=True, name = "ZZ_aiso", isWSR=True, additionalnorm= qcdnorm * mcscale)
        # tau
        ZXXAisoSamp = Sample(input_antiiso_zxx, isMC=True, name = "ZXX_aiso", isWSR=True, additionalnorm= qcdnorm * mcscale)
        Wx0AisoSamp = Sample(input_antiiso_wx0, isMC=True, name = "wx0_aiso", isWSR=True, additionalnorm= qcdnorm * mcscale)
        Wx1AisoSamp = Sample(input_antiiso_wx1, isMC=True, name = "wx1_aiso", isWSR=True, additionalnorm= qcdnorm * mcscale)
        Wx2AisoSamp = Sample(input_antiiso_wx2, isMC=True, name = "wx2_aiso", isWSR=True, additionalnorm= qcdnorm * mcscale)

        sampMan = SampleManager(DataAisoSamp, [Wl0AisoSamp, Wl1AisoSamp, Wl2AisoSamp, TTbarAisoSamp, TT1LepAisoSamp, TT0LepAisoSamp, WWAisoSamp, WZAisoSamp, ZZAisoSamp, ZXXAisoSamp, Wx0AisoSamp, Wx1AisoSamp, Wx2AisoSamp])

        sampMan.groupMCs(["WW_aiso", "WZ_aiso", "ZZ_aiso", "ZXX_aiso", "wx0_aiso", "wx1_aiso", "wx2_aiso"], "EWK", 216, "EWK")
        sampMan.groupMCs(["ttbar_dilepton_aiso", "ttbar_1lepton_aiso", "ttbar_0lepton_aiso"], "ttbar", 96, "t#bar{t}")
        label = "W#rightarrow#mu#nu" if doMuon else "W#rightarrow e#nu"
        sampMan.groupMCs(['wl0_aiso', 'wl1_aiso', 'wl2_aiso'], "wlnu", 92, label)

    else:
        # W -> lnu
        label = "W#rightarrow#mu#nu" if doMuon else "W#rightarrow e#nu"
        WlAisoSamp   = Sample(input_antiiso_wl, isMC=True, name = "wlnu", isWSR=True, additionalnorm= qcdnorm * mcscale, is5TeV = True, color=92, legend=label)
        # ttbar
        TTbarAisoSamp  = Sample(input_antiiso_ttbar, isMC=True, name = "ttbar",     isWSR=True, additionalnorm= qcdnorm * mcscale, is5TeV = True, color=96, legend="t#bar{t}")
        ## dibosons
        WWAisoSamp = Sample(input_antiiso_ww, isMC=True, name = "WW_aiso", isWSR=True, additionalnorm= qcdnorm * mcscale, is5TeV = True)
        WZAisoSamp = Sample(input_antiiso_wz, isMC=True, name = "WZ_aiso", isWSR=True, additionalnorm= qcdnorm * mcscale, is5TeV = True)
        ZZ2LAisoSamp = Sample(input_antiiso_zz2l, isMC=True, name = "ZZ2L_aiso", isWSR=True, additionalnorm= qcdnorm * mcscale, is5TeV = True)
        ZZ4LAisoSamp = Sample(input_antiiso_zz4l, isMC=True, name = "ZZ4L_aiso", isWSR=True, additionalnorm= qcdnorm * mcscale, is5TeV = True)
        # tau
        ZXXAisoSamp = Sample(input_antiiso_zxx, isMC=True, name = "ZXX_aiso", isWSR=True, additionalnorm= qcdnorm * mcscale, is5TeV = True)
        WxAisoSamp  = Sample(input_antiiso_wx,  isMC=True, name = "wx_aiso",  isWSR=True, additionalnorm= qcdnorm * mcscale, is5TeV = True)

        sampMan = SampleManager(DataAisoSamp, [WlAisoSamp, TTbarAisoSamp, WWAisoSamp, WZAisoSamp, ZZ2LAisoSamp, ZZ4LAisoSamp, ZXXAisoSamp, WxAisoSamp], is5TeV = True)

        sampMan.groupMCs(["WW_aiso", "WZ_aiso", "ZZ2L_aiso", "ZZ4L_aiso", "ZXX_aiso", "wx_aiso", ], "EWK", 216, "EWK")

    # customize legends
    label_plus = "W^{+}#rightarrow#mu^{+}#nu" if doMuon else "W^{+}#rightarrow e^{+}#nu"
    label_minus = "W^{-}#rightarrow#mu^{-}#bar{#nu}" if doMuon else "W^{-}#rightarrow e^{-}#bar{#nu}"
    legends = ["Data", label_plus, "t#bar{t}", "EWK"]

    sampMan.DefineAll("Lep_pt", "lep.Pt()")
    sampMan.ApplyCutAll("Lep_pt > 25.0")
    sampMan.DefineAll("Lep_eta", "lep.Eta()")

    # muon and electron isolation distributions are different
    # more coarse binning for electrons to make sure enough statistics
    if doMuon:
        isoCuts = [0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95]
        isobins = []
        sampMan.DefineAll("RelIso", "relIso")
        for isobin in range(4, 20):
            sampMan.DefineAll(f"w_iso{isobin}", f"(relIso > {isoCuts[isobin-4]} && relIso < {isoCuts[isobin-3]})")
            isobins.append(f"iso{isobin}")

    else:
        isoCuts = [0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.50, 0.65, 0.90]
        isobins = []
        sampMan.DefineAll("isEB",   "fabs(Lep_eta) <= 1.4442")
        #sampMan.DefineAll("RelIso", "isEB ? (relIso + 0.0287 - 0.0478) : (relIso + 0.0445 - 0.0658)")
        sampMan.DefineAll("RelIso", "(pfCombIso/lep.Pt())")
        for isobin in range(4, 13):
            sampMan.DefineAll(f"w_iso{isobin}", f"(RelIso > {isoCuts[isobin-4]} && RelIso < {isoCuts[isobin-3]})")
            isobins.append(f"iso{isobin}")

    sampMan.DefineMC("met_pt", "metVars[1]", excludes=['Data_aiso'])
    sampMan.DefineMC("met_phi", "TVector2::Phi_mpi_pi(metVarsPhi[1])", excludes=['Data_aiso'])
    # Data does not run any recoil corrections
    DataAisoSamp.Define("met_pt", "metVars[0]")
    DataAisoSamp.Define("met_phi", "TVector2::Phi_mpi_pi(metVarsPhi[0])")

    # recoil variables
    sampMan.DefineAll("V2W", "UVec(lep.Pt(), lep.Phi(), met_pt, met_phi)")
    sampMan.DefineAll("WpT", "V2W.Mod()")
    sampMan.DefineAll("Wphi", "TVector2::Phi_mpi_pi(V2W.Phi())")

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
        wptbins = ["WpT_bin0", "WpT_bin1", "WpT_bin2", "WpT_bin3", "WpT_bin4", "WpT_bin5", "WpT_bin6", "WpT_bin7", "WpT_bin8", "WpT_bin9"]

    # eta bins for electrons: barral and endcap
    sampMan.DefineAll("lepEta_bin0", "1.0")
    sampMan.DefineAll("lepEta_bin1", "abs(lep.Eta()) <= 1.4442")
    sampMan.DefineAll("lepEta_bin2", "abs(lep.Eta()) > 1.4442")

    if doMuon:
        etabins = ["lepEta_bin0"]
    else:
        etabins = ["lepEta_bin0", "lepEta_bin1", "lepEta_bin2"]

    # draw the lepton isolation distribution
    nbins = 50
    xmin = 0.
    xmax = 0.75
    ymax = 5e7 if doMuon else 1e6
    for wpt in wptbins:
        for lepeta in etabins:
            strname = f"weight_{wpt}_{lepeta}"
            sampMan.DefineAll(strname, f"weight_WoVpt * {wpt} * {lepeta}")
            sampMan.cacheDraw("RelIso", f"histo_wjets_{lepname}_RelIso_{lepeta}_{wpt}", 100, 0, 0.72, DrawConfig(xmin=xmin, xmax=xmax, xlabel="Relative Isolation", ylabel=f"Events / {(xmax-xmin)/nbins:.2f}", dology=True, ymax=ymax, donormalizebin=False, addOverflow=True, addUnderflow=True, showratio=False, legendPos=[0.94, 0.88, 0.70, 0.68]), weightname = strname)

    # do the fine binning first; then rebin in the processHists
    mass_bins = np.array([0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0, 95.0, 100.0, 105.0, 110.0, 115.0, 120.0, 125.0, 130.0, 135.0, 140.])
    nbins = 12
    xmin = 0
    xmax = 120
    idx = 0
    for iso in isobins:
        for wpt in wptbins:
            for lepeta in etabins:
                for chg in chgbins:
                    strname = "weight_{}_{}_{}_{}".format(chg, iso, wpt, lepeta)

                    sampMan.DefineAll(strname, "w_{} * weight_WoVpt * {} * {} * {}".format(iso, wpt, lepeta, chg))
                    legends[1] = label_plus if "plus" in chg else label_minus

                    outputname = "histo_wjetsAntiIso_mtcorr_" + strname 
                    sampMan.cacheDraw("mtCorr", outputname, mass_bins, DrawConfig(xmin=xmin, xmax=xmax, xlabel="m_{T} [GeV]", ylabel=f"Events / {(xmax-xmin)/nbins:.0f} GeV", dology=True, ymax=6e5, donormalizebin=False, addOverflow=True, addUnderflow=True, showratio=False, lheader = f"{isoCuts[idx]} < I < {isoCuts[idx + 1]}", legendPos=[0.94, 0.88, 0.70, 0.68], legends = legends.copy()), weightname = strname)
                    #outputname = "histo_wjetsAntiIso_lepEta_" + strname
                    #sampMan.cacheDraw("Lep_eta", outputname, 24, -2.4, 2.4, DrawConfig(xmin=-2.4, xmax=2.4, xlabel="Lepton Eta", dology=True, ymax=2e6, donormalizebin=False, addOverflow=False, addUnderflow=False, showratio=True), weightname = strname)

                    # draw the isolation in different bins, in order to calculate the mean in each bin
                    sampMan.cacheDraw("RelIso", f"histo_wjetsAntiIso_RelIso_{strname}", 100, isoCuts[idx], isoCuts[idx + 1], DrawConfig(xmin=xmin, xmax=xmax, xlabel="Relative Isolation", ylabel=f"Events / {(xmax-xmin)/nbins:.2f}", dology=True, ymax=ymax, donormalizebin=False, addOverflow=True, addUnderflow=True, showratio=False, legendPos=[0.94, 0.88, 0.70, 0.68]), weightname = strname)
        idx += 1

    sampMan.launchDraw()


    hmts_comp = OrderedDict()
    hIsos = OrderedDict()

    #hetas_mtCut_comp = OrderedDict()
    for iso in isobins:
        for wpt in wptbins:
            for lepeta in etabins:
                for chg in chgbins:
                    strname = "weight_{}_{}_{}_{}".format(chg, iso, wpt, lepeta)
                    
                    # for mT
                    outputname = "histo_wjetsAntiIso_mtcorr_" + strname
                    hstacked = THStack2TH1(sampMan.hsmcs[outputname])
                    for ibin in range(hstacked.GetNbinsX()+1):
                        # hstacked should always be above 0
                        hstacked.SetBinContent(ibin, max(hstacked.GetBinContent(ibin), 0))
                    sampMan.hdatas[outputname].Add( hstacked, -1.0 )
                    hmts_comp[strname] = sampMan.hdatas[outputname]
                    hmts_comp[strname].SetName(outputname)

                    ## for lepton eta with mT cuts
                    #strname += "_mtCut"
                    #outputname = "histo_wjetsAntiIso_lepEta_" + strname
                    #hstacked = THStack2TH1(sampMan.hsmcs[outputname])
                    #for ibin in xrange(hstacked.GetNbinsX()+1):
                    #    # hstacked should always be above 0
                    #    hstacked.SetBinContent(ibin, max(hstacked.GetBinContent(ibin), 0))
                    #sampMan.hdatas[outputname].Add( hstacked, -1.0 )
                    #hetas_mtCut_comp[strname] = sampMan.hdatas[outputname]
                    #hetas_mtCut_comp[strname].SetName(outputname)

                    # for isolation
                    outputname = f"histo_wjetsAntiIso_RelIso_{strname}"
                    hstacked = THStack2TH1(sampMan.hsmcs[outputname])
                    for ibin in range(hstacked.GetNbinsX()+1):
                        # hstacked should always be above 0
                        hstacked.SetBinContent(ibin, max(hstacked.GetBinContent(ibin), 0))
                    sampMan.hdatas[outputname].Add( hstacked, -1.0 )
                    hIsos[strname] = sampMan.hdatas[outputname]
                    hIsos[strname].SetName(outputname)

    postfix = lepname + "nu"
    if applyScaling:
        postfix += "_applyScaling"
    if do5TeV:
        postfix += "_5TeV"
    postfix += ".root"
    outfile = ROOT.TFile.Open("root/output_qcdshape_"+postfix, "recreate")

    for wpt in wptbins:
        #odir = outfile.mkdir(wpt)
        #outfile.cd(wpt)
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
                    strname = "weight_{}_{}_{}_{}".format(chg, iso, wpt, lepeta)
                    strname_next = "weight_{}_{}_{}_{}".format(chg, iso_next, wpt, lepeta)

                    outputname = "histo_wjetsAntiIso_mtcorr_" + strname

                    hcenter = hmts_comp[strname]
                    # shape uncertaintis are using the shape difference from the neighboring bins
                    # used in the original qcd background predictions
                    hup = hmts_comp[strname_next].Clone(outputname+"_shapeUp")
                    hdown = hmts_comp[strname_next].Clone(outputname+"_shapeDown")

                    hup.Scale(hcenter.Integral() / (hup.Integral()+1e-6))

                    for ibin in range(1, hcenter.GetNbinsX()+1):
                        center = hcenter.GetBinContent(ibin)
                        up = hup.GetBinContent(ibin)
                        hdown.SetBinContent(ibin, max(2*center - up, 0))

                        hcenter.SetBinContent(ibin, max(center, 0))
                        hup.SetBinContent(ibin, max(up, 0))

                    #hcenter.SetDirectory(odir)
                    hcenter.SetDirectory(outfile)
                    hcenter.Write()
                    #hup.SetDirectory(odir)
                    hup.SetDirectory(outfile)
                    hup.Write()
                    #hdown.SetDirectory(odir)
                    hdown.SetDirectory(outfile)
                    hdown.Write()

    outfile.Close()


    sampMan.dumpCounts()

    outfile = ROOT.TFile.Open("root/output_qcdIso_"+postfix, "recreate")
    for isobin, h in hIsos.items():
        print(f"{isobin} mean: {h.GetMean():.3f}")
        h.SetDirectory(outfile)
        h.Write()
    outfile.Close()

    print("Program end...")

    input()
    
    return 

if __name__ == "__main__":
   main()

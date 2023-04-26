import ROOT
import re
import numpy as np
import os
from CMSPLOTS.myFunction import AddOverflowsTH1, RebinHisto
from modules.qcdExtrapolater import ExtrapolateQCD, LoadQCDNorms, InclusiveQCD
from modules.cardMaker import MakeWJetsCards, MakeZJetsCards, GenerateRunCommand, MakeXSecCard
from modules.histProcessor import ProcessHists, CopyandMergeTau
from modules.Binnings import mass_bins_w, mass_bins_z, mass_bins_test
from modules.Utils import GetValues
from collections import OrderedDict
import argparse

ROOT.gROOT.SetBatch(True)

def RunPreparations(fwsig_input, fwsig_rebin, fwsig_mergeTau, fqcd_input, fqcd_rebin, fqcd_input_scaled, fqcd_rebin_scaled, fqcd_output, lepname="mu", is5TeV=False, mass_bins_w=mass_bins_w, outdir_card="cards", fzsig_input=None, fzsig_output=None, mass_bins_z=mass_bins_z, applyLFU=False, doInclusive=False):
    """
    rebin the mT hists for sig and qcd bkg,
    merge tau (wx) process to sig (wlnu),
    extrapolate the qcd shape from anti-isolated to isolated
    genereate datacard for W;
    rebin the mll hists for z's;
    generate datacard for Z;
    return the list of datacard names
    """
    includeUnderflow = False
    includeOverflow = True

    ProcessHists(fwsig_input, fwsig_rebin, mass_bins_w,
                 includeUnderflow, includeOverflow)
    # merge the tau processes into the signal process
    CopyandMergeTau(fwsig_rebin, fwsig_mergeTau)

    # for QCD
    #ProcessHists(fqcd_input, fqcd_rebin, mass_bins_w,
    #             includeUnderflow, includeOverflow)
    #if fqcd_input_scaled and fqcd_rebin_scaled:
    #    ProcessHists(fqcd_input_scaled, fqcd_rebin_scaled,
    #                 mass_bins_w, includeUnderflow, includeOverflow)

    ## extrapolate the QCD template from anti-isolated region to isolated region
    ## manually set the QCD normalization factors for the prefit impacts
    ## very hacky, not ideal at all; only temporary solution
    #fqcdnorm = "data/QCDNorms_inc.json"
    #qcdnorms = LoadQCDNorms(fqcdnorm)
    #print("QCD Norms: ", qcdnorms)
    ## qcdnorms = {}
    #ExtrapolateQCD(fqcd_rebin, fqcd_output, [lepname+"plus", lepname+"minus"], "WpT_bin0", [
    #               "lepEta_bin0"], fname_scaled=fqcd_rebin_scaled, rebinned=True, is5TeV=is5TeV, qcdnorms=qcdnorms)
    ##InclusiveQCD(fqcd_rebin, fqcd_output, [lepname+"plus", lepname+"minus"], "WpT_bin0", [
    ##                "lepEta_bin0"], fname_scaled=fqcd_rebin_scaled, rebinned=True, is5TeV=is5TeV, qcdnorms=qcdnorms)
    
    # copy fqcd_input to fqcd_output
    os.system("cp "+fqcd_input+" "+fqcd_output)

    # generate card based on the signal and qcd templates
    card_plus = MakeWJetsCards(fwsig_mergeTau, fqcd_output, lepname+"plus", "WpT_bin0", "lepEta_bin0",
                               rebinned=True, nMTBins=len(mass_bins_w)-1, is5TeV=is5TeV, outdir=outdir_card, applyLFU=applyLFU, doInclusive=doInclusive)
    card_minus = MakeWJetsCards(fwsig_mergeTau, fqcd_output, lepname+"minus", "WpT_bin0", "lepEta_bin0",
                                rebinned=True, nMTBins=len(mass_bins_w)-1, is5TeV=is5TeV, outdir=outdir_card, applyLFU=applyLFU, doInclusive=doInclusive)

    # generate xsec hists and cards
    outdir_xsec_root = fwsig_mergeTau.rpartition("/")[0]
    card_xsec_plus = MakeXSecCard(
        lepname+"plus", is5TeV, outdir_root=outdir_xsec_root, outdir_card=outdir_card, applyLFU=applyLFU)
    card_xsec_minus = MakeXSecCard(
        lepname+"minus", is5TeV, outdir_root=outdir_xsec_root, outdir_card=outdir_card, applyLFU=applyLFU)

    card_z = None
    if fzsig_input:
        # rebin the z histograms
        # for Z's should not include the overflow or underflows, since there is the mll cut
        ProcessHists(fzsig_input, fzsig_output, mass_bins_z, False, False)
        # generate z card
        card_z = MakeZJetsCards(fzsig_output, lepname+lepname, rebinned=True,
                                is5TeV=is5TeV, outdir=outdir_card, applyLFU=applyLFU, doInclusive=doInclusive)
        card_xsec_z = MakeXSecCard(
            lepname+lepname, is5TeV, outdir_root=outdir_xsec_root, outdir_card=outdir_card, applyLFU=applyLFU)

    return card_plus, card_minus, card_z, card_xsec_plus, card_xsec_minus, card_xsec_z


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Make postfit plots and tables")
    parser.add_argument("--doInclusive", action="store_true", dest="doInclusive",help="Run on inclusive channel; default is fiducial")
    parser.add_argument("--do13TeV", action="store_true", dest="do13TeV",help="Run on 13 TeV channel;")
    parser.add_argument("--do5TeV", action="store_true", dest="do5TeV",help="Run on 5 TeV channel;")
    
    args = parser.parse_args()
    doInclusive = args.doInclusive
    do13TeV = args.do13TeV
    do5TeV = args.do5TeV
    
    combineSqrtS = (do13TeV and do5TeV)
    sqrtSs = []
    if do13TeV:
        sqrtSs.append("13TeV")
    if do5TeV:
        sqrtSs.append("5TeV")
    
    # apply lepton flavor university
    applyLFU = True

    # scan fit range
    for key, val in mass_bins_test.items():
        if key != 0:
            continue

        cards = OrderedDict()
        channels = OrderedDict()
        cards_xsec = OrderedDict()

        for sqrtS in sqrtSs:
            card_muplus = None
            card_muminus = None
            card_zmumu = None
            card_eplus = None
            card_eminus = None
            card_zee = None

            suffix = f"_{sqrtS}"
            if doInclusive:
                suffix += "_noTheoryNorm"
            else:
                suffix += "_TheoryNormed"
            suffix += ".root"
            suffix_qcd = f"_{sqrtS}.root"
            is5TeV = (sqrtS == "5TeV")

            fzsig_mumu_input = "root/output_shapes_mumu" + suffix
            fzsig_ee_input = "root/output_shapes_ee" + suffix

            fwsig_munu_input = "root/output_shapes_munu" + suffix
            fwsig_enu_input = "root/output_shapes_enu" + suffix

            odir = "Inclusive" if doInclusive else "Fiducial"

            # muon channel
            fwsig_rebin = f"root/{odir}/test{key}/{sqrtS}/output_shapes_munu_Rebin" + suffix
            fwsig_mergeTau = f"forCombine/{odir}/test{key}/root/{sqrtS}/output_shapes_munu_mergeTau" + suffix

            fqcd_input = "root/output_qcdshape_munu" + suffix_qcd
            fqcd_rebin = f"root/{odir}/test{key}/{sqrtS}/output_qcdshape_munu_Rebin" + suffix_qcd
            fqcd_input_scaled = "root/output_qcdshape_munu_applyScaling" + suffix_qcd
            fqcd_rebin_scaled = f"root/{odir}/test{key}/{sqrtS}/output_qcdshape_munu_Rebin_applyScaling" + suffix_qcd
            fqcd_input_scaled = None
            fqcd_rebin_scaled = None
            fqcd_input = f"plots/{sqrtS}/FRAndClosure/relIso/mu_pt_vs_eta_TwoPhis/histos_qcdFR_mu_{sqrtS}.root"
            fqcd_output = f"forCombine/{odir}/test{key}/root/{sqrtS}/qcdshape_extrapolated_munu" + suffix_qcd
            fqcd_output = f"forCombine/{odir}/test{key}/root/{sqrtS}/histos_qcdFR_mu_{sqrtS}.root"

            fzsig_output = f"forCombine/{odir}/test{key}/root/{sqrtS}/output_shapes_mumu_Rebin" + suffix

            card_muplus, card_muminus, card_zmumu, card_xsec_muplus, card_xsec_muminus, card_xsec_mumu = RunPreparations(
                fwsig_munu_input, fwsig_rebin, fwsig_mergeTau, fqcd_input, fqcd_rebin, fqcd_input_scaled, fqcd_rebin_scaled, fqcd_output, "mu", outdir_card=f"forCombine/{odir}/test{key}/cards/{sqrtS}", mass_bins_w=val, fzsig_input=fzsig_mumu_input, fzsig_output=fzsig_output, applyLFU=applyLFU, is5TeV=is5TeV, doInclusive=doInclusive)

            cards[f"mu_{sqrtS}"] = [card_muplus, card_muminus, card_zmumu]
            cards_xsec[f"mu_{sqrtS}"] = [
                card_xsec_muplus, card_xsec_muminus, card_xsec_mumu]
            channels[f"mu_{sqrtS}"] = [f"muplus_{sqrtS}",
                                       f"muminus_{sqrtS}", f"mumu_{sqrtS}"]

            # electron channel
            fwsig_rebin = f"root/{odir}/test{key}/{sqrtS}/output_shapes_enu_Rebin" + suffix
            fwsig_mergeTau = f"forCombine/{odir}/test{key}/root/{sqrtS}/output_shapes_enu_mergeTau" + suffix

            fqcd_input = "root/output_qcdshape_enu" + suffix_qcd
            fqcd_rebin = f"root/{odir}/test{key}/{sqrtS}/output_qcdshape_enu_Rebin" + suffix_qcd
            fqcd_input_scaled = "root/output_qcdshape_enu_applyScaling" + suffix_qcd
            fqcd_rebin_scaled = f"root/{odir}/test{key}/{sqrtS}/output_qcdshape_enu_Rebin_applyScaling" + suffix_qcd
            fqcd_input_scaled = None
            fqcd_rebin_scaled = None
            fqcd_output = f"forCombine/{odir}/test{key}/root/{sqrtS}/qcdshape_extrapolated_enu" + suffix_qcd
            fqcd_input = f"plots/{sqrtS}/FRAndClosure/relIso/e_pt_vs_eta_TwoPhis/histos_qcdFR_e_{sqrtS}.root"
            fqcd_output = f"forCombine/{odir}/test{key}/root/{sqrtS}/histos_qcdFR_e_{sqrtS}.root"

            fzsig_output = f"forCombine/{odir}/test{key}/root/{sqrtS}/output_shapes_ee_Rebin" + suffix

            card_eplus, card_eminus, card_zee, card_xsec_eplus, card_xsec_eminus, card_xsec_ee = RunPreparations(
                fwsig_enu_input, fwsig_rebin, fwsig_mergeTau, fqcd_input, fqcd_rebin, fqcd_input_scaled, fqcd_rebin_scaled, fqcd_output, "e", outdir_card=f"forCombine/{odir}/test{key}/cards/{sqrtS}", mass_bins_w=val, fzsig_input=fzsig_ee_input, fzsig_output=fzsig_output, applyLFU=applyLFU, is5TeV=is5TeV, doInclusive=doInclusive)

            cards[f"e_{sqrtS}"] = [card_eplus, card_eminus, card_zee]
            cards_xsec[f"e_{sqrtS}"] = [
                card_xsec_eplus, card_xsec_eminus, card_xsec_ee]
            channels[f"e_{sqrtS}"] = [f"eplus_{sqrtS}",
                                      f"eminus_{sqrtS}", f"ee_{sqrtS}"]

            # combined fit
            GenerateRunCommand(f"card_combined_{sqrtS}", GetValues(cards, f".*_{sqrtS}"), GetValues(channels, f".*_{sqrtS}"), GetValues(
                cards_xsec, f".*_{sqrtS}"), output_dir=f"forCombine/{odir}/test{key}/commands/scripts/", applyLFU=applyLFU, is5TeV=is5TeV)

            # mu channel fit
            GenerateRunCommand(f"card_mu_{sqrtS}", GetValues(cards, f"mu_{sqrtS}"), GetValues(channels, f"mu_{sqrtS}"), GetValues(
                cards_xsec, f"mu_{sqrtS}"), output_dir=f"forCombine/{odir}/test{key}/commands/scripts/", applyLFU=applyLFU, is5TeV=is5TeV)

            # e channel fit
            GenerateRunCommand(f"card_e_{sqrtS}", GetValues(cards, f"e_{sqrtS}"), GetValues(channels, f"e_{sqrtS}"), GetValues(
                cards_xsec, f"e_{sqrtS}"), output_dir=f"forCombine/{odir}/test{key}/commands/scripts/", applyLFU=applyLFU, is5TeV=is5TeV)

        if combineSqrtS:
            # combine 13 and 5 TeV
            GenerateRunCommand("card_combined", GetValues(cards, ".*"), GetValues(channels, ".*"), GetValues(cards_xsec, ".*"),
                               output_dir=f"forCombine/{odir}/test{key}/commands/scripts/", applyLFU=applyLFU, combineAll=True)

            # mu channel combined fit
            GenerateRunCommand("card_mu", GetValues(cards, "mu.*"), GetValues(channels, "mu.*"), GetValues(cards_xsec, "mu.*"),
                               output_dir=f"forCombine/{odir}/test{key}/commands/scripts/", applyLFU=applyLFU, combineAll=True)

            # e channel combined fit
            GenerateRunCommand("card_e", GetValues(cards, "e.*"), GetValues(channels, "e.*"), GetValues(cards_xsec, "e.*"),
                               output_dir=f"forCombine/{odir}/test{key}/commands/scripts/", applyLFU=applyLFU, combineAll=True)

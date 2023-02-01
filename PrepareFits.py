import ROOT
import re
import numpy as np
from CMSPLOTS.myFunction import AddOverflowsTH1, RebinHisto
from modules.qcdExtrapolater import ExtrapolateQCD, LoadQCDNorms
from modules.cardMaker import MakeWJetsCards, MakeZJetsCards, GenerateRunCommand, MakeXSecCard
from modules.histProcessor import ProcessHists, CopyandMergeTau
from modules.Binnings import mass_bins_w, mass_bins_z, mass_bins_test
from modules.Utils import GetValues
from collections import OrderedDict

ROOT.gROOT.SetBatch(True)

def RunPreparations(fwsig_input, fwsig_rebin, fwsig_mergeTau, fqcd_input, fqcd_rebin, fqcd_input_scaled, fqcd_rebin_scaled, fqcd_output, lepname = "mu", is5TeV = False, mass_bins_w = mass_bins_w, outdir_card = "cards", fzsig_input = None, fzsig_output = None, mass_bins_z = mass_bins_z, applyLFU = False):
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

    ProcessHists(fwsig_input, fwsig_rebin, mass_bins_w, includeUnderflow, includeOverflow)
    # merge the tau processes into the signal process
    CopyandMergeTau(fwsig_rebin, fwsig_mergeTau)

    # for QCD
    ProcessHists(fqcd_input, fqcd_rebin, mass_bins_w, includeUnderflow, includeOverflow)
    ProcessHists(fqcd_input_scaled, fqcd_rebin_scaled, mass_bins_w, includeUnderflow, includeOverflow)

    # extrapolate the QCD template from anti-isolated region to isolated region
    # manually set the QCD normalization factors for the prefit impacts
    # very hacky, not ideal at all; only temporary solution
    fqcdnorm = "data/QCDNorms.json"
    qcdnorms = LoadQCDNorms(fqcdnorm)
    print("QCD Norms: ", qcdnorms)
    #qcdnorms = {}
    ExtrapolateQCD(fqcd_rebin, fqcd_output, [lepname+"plus", lepname+"minus"], "WpT_bin0", ["lepEta_bin0"], fname_scaled=fqcd_rebin_scaled, rebinned = True, is5TeV = is5TeV, qcdnorms=qcdnorms)

    # generate card based on the signal and qcd templates
    card_plus = MakeWJetsCards(fwsig_mergeTau, fqcd_output, lepname+"plus", "WpT_bin0", "lepEta_bin0", rebinned = True, nMTBins = len(mass_bins_w)-1, is5TeV = is5TeV, outdir = outdir_card, applyLFU=applyLFU)
    card_minus = MakeWJetsCards(fwsig_mergeTau, fqcd_output, lepname+"minus", "WpT_bin0", "lepEta_bin0", rebinned = True, nMTBins = len(mass_bins_w)-1, is5TeV = is5TeV, outdir = outdir_card, applyLFU=applyLFU)

    # generate xsec hists and cards
    outdir_xsec_root = fwsig_mergeTau.rpartition("/")[0]
    card_xsec_plus = MakeXSecCard(lepname+"plus", is5TeV, outdir_root = outdir_xsec_root, outdir_card = outdir_card, applyLFU=applyLFU)
    card_xsec_minus = MakeXSecCard(lepname+"minus", is5TeV, outdir_root = outdir_xsec_root, outdir_card = outdir_card, applyLFU=applyLFU)

    card_z = None
    if fzsig_input:
        # rebin the z histograms
        # for Z's should not include the overflow or underflows, since there is the mll cut
        ProcessHists(fzsig_input, fzsig_output, mass_bins_z, False, False)
        # generate z card
        card_z = MakeZJetsCards(fzsig_output, lepname+lepname, rebinned = True, is5TeV = is5TeV, outdir = outdir_card, applyLFU=applyLFU)
        card_xsec_z = MakeXSecCard(lepname+lepname, is5TeV, outdir_root = outdir_xsec_root, outdir_card = outdir_card, applyLFU=applyLFU)

    return card_plus, card_minus, card_z, card_xsec_plus, card_xsec_minus, card_xsec_z

if __name__  == "__main__":
    # apply lepton flavor university
    applyLFU = True
    # inclusive xsec or fiducial xsec
    doInclusive = True

    # scan fit range
    for key, val in mass_bins_test.items():
        if key != 0:
            continue

        cards = OrderedDict()
        channels = OrderedDict()
        cards_xsec = OrderedDict()

        for sqrtS in ["13TeV", "5TeV"]:
            #if sqrtS != "13TeV":
            #    continue
            card_muplus = None
            card_muminus = None
            card_zmumu = None
            card_eplus = None
            card_eminus = None
            card_zee = None

            suffix = f"_{sqrtS}"
            if doInclusive:
                suffix += "_noTheoryNorm"
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
            fqcd_output = f"forCombine/{odir}/test{key}/root/{sqrtS}/qcdshape_extrapolated_munu" + suffix_qcd

            fzsig_output = f"forCombine/{odir}/test{key}/root/{sqrtS}/output_shapes_mumu_Rebin" + suffix

            card_muplus, card_muminus, card_zmumu, card_xsec_muplus, card_xsec_muminus, card_xsec_mumu = RunPreparations(fwsig_munu_input, fwsig_rebin, fwsig_mergeTau, fqcd_input, fqcd_rebin, fqcd_input_scaled, fqcd_rebin_scaled, fqcd_output, "mu", outdir_card = f"forCombine/{odir}/test{key}/cards/{sqrtS}", mass_bins_w = val, fzsig_input=fzsig_mumu_input, fzsig_output = fzsig_output, applyLFU=applyLFU, is5TeV=is5TeV)

            cards[f"mu_{sqrtS}"] = [card_muplus, card_muminus, card_zmumu]
            cards_xsec[f"mu_{sqrtS}"] = [card_xsec_muplus, card_xsec_muminus, card_xsec_mumu]
            channels[f"mu_{sqrtS}"] = [f"muplus_{sqrtS}", f"muminus_{sqrtS}", f"mumu_{sqrtS}"]

            # electron channel
            fwsig_rebin = f"root/{odir}/test{key}/{sqrtS}/output_shapes_enu_Rebin" + suffix
            fwsig_mergeTau = f"forCombine/{odir}/test{key}/root/{sqrtS}/output_shapes_enu_mergeTau" + suffix

            fqcd_input = "root/output_qcdshape_enu" + suffix_qcd
            fqcd_rebin = f"root/{odir}/test{key}/{sqrtS}/output_qcdshape_enu_Rebin" + suffix_qcd
            fqcd_input_scaled = "root/output_qcdshape_enu_applyScaling" + suffix_qcd
            fqcd_rebin_scaled = f"root/{odir}/test{key}/{sqrtS}/output_qcdshape_enu_Rebin_applyScaling" + suffix_qcd
            fqcd_output = f"forCombine/{odir}/test{key}/root/{sqrtS}/qcdshape_extrapolated_enu" + suffix_qcd

            fzsig_output = f"forCombine/{odir}/test{key}/root/{sqrtS}/output_shapes_ee_Rebin" + suffix

            card_eplus, card_eminus, card_zee, card_xsec_eplus, card_xsec_eminus, card_xsec_ee = RunPreparations(fwsig_enu_input, fwsig_rebin, fwsig_mergeTau, fqcd_input, fqcd_rebin, fqcd_input_scaled, fqcd_rebin_scaled, fqcd_output, "e", outdir_card = f"forCombine/{odir}/test{key}/cards/{sqrtS}", mass_bins_w = val, fzsig_input=fzsig_ee_input, fzsig_output=fzsig_output, applyLFU=applyLFU, is5TeV=is5TeV)

            cards[f"e_{sqrtS}"] = [card_eplus, card_eminus, card_zee]
            cards_xsec[f"e_{sqrtS}"] = [card_xsec_eplus, card_xsec_eminus, card_xsec_ee]
            channels[f"e_{sqrtS}"] = [f"eplus_{sqrtS}", f"eminus_{sqrtS}", f"ee_{sqrtS}"]

            # combined fit
            GenerateRunCommand(f"card_combined_{sqrtS}", GetValues(cards, f".*_{sqrtS}"), GetValues(channels, f".*_{sqrtS}"), GetValues(cards_xsec, f".*_{sqrtS}"), output_dir = f"forCombine/{odir}/test{key}/commands/scripts/", applyLFU=applyLFU, is5TeV=is5TeV)

            # mu channel fit
            GenerateRunCommand(f"card_mu_{sqrtS}", GetValues(cards, f"mu_{sqrtS}"), GetValues(channels, f"mu_{sqrtS}"), GetValues(cards_xsec, f"mu_{sqrtS}"), output_dir = f"forCombine/{odir}/test{key}/commands/scripts/", applyLFU=applyLFU, is5TeV=is5TeV)

            # e channel fit
            GenerateRunCommand(f"card_e_{sqrtS}", GetValues(cards, f"e_{sqrtS}"), GetValues(channels, f"e_{sqrtS}"), GetValues(cards_xsec, f"e_{sqrtS}"), output_dir = f"forCombine/{odir}/test{key}/commands/scripts/", applyLFU=applyLFU, is5TeV=is5TeV)

        # combine 13 and 5 TeV
        GenerateRunCommand("card_combined", GetValues(cards, ".*"), GetValues(channels, ".*"), GetValues(cards_xsec, ".*"), output_dir = f"forCombine/{odir}/test{key}/commands/scripts/", applyLFU=applyLFU, combineAll=True)

        # mu channel combined fit
        GenerateRunCommand("card_mu", GetValues(cards, "mu.*"), GetValues(channels, "mu.*"), GetValues(cards_xsec, "mu.*"), output_dir = f"forCombine/{odir}/test{key}/commands/scripts/", applyLFU=applyLFU, combineAll=True)

        # e channel combined fit
        GenerateRunCommand("card_e", GetValues(cards, "e.*"), GetValues(channels, "e.*"), GetValues(cards_xsec, "e.*"), output_dir = f"forCombine/{odir}/test{key}/commands/scripts/", applyLFU=applyLFU, combineAll=True)

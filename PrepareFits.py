import ROOT
import re
import numpy as np
from CMSPLOTS.myFunction import AddOverflowsTH1, RebinHisto
from modules.qcdExtrapolater import ExtrapolateQCD
from modules.cardMaker import MakeWJetsCards, MakeZJetsCards, GenerateRunCommand, MakeXSecCard
from modules.histProcessor import ProcessHists, CopyandMergeTau
from modules.mass_bins import mass_bins_w, mass_bins_z
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
    ExtrapolateQCD(fqcd_rebin, fqcd_output, [lepname+"plus", lepname+"minus"], "WpT_bin0", ["lepEta_bin0"], fname_scaled=fqcd_rebin_scaled, rebinned = True, is5TeV = is5TeV)

    # generate card based on the signal and qcd templates
    card_plus = MakeWJetsCards(fwsig_mergeTau, fqcd_output, lepname+"plus", "WpT_bin0", "lepEta_bin0", rebinned = True, nMTBins = len(mass_bins_w)-1, is5TeV = is5TeV, outdir = outdir_card, applyLFU=applyLFU)
    card_minus = MakeWJetsCards(fwsig_mergeTau, fqcd_output, lepname+"minus", "WpT_bin0", "lepEta_bin0", rebinned = True, nMTBins = len(mass_bins_w)-1, is5TeV = is5TeV, outdir = outdir_card, applyLFU=applyLFU)

    # generate xsec hists and cards
    card_xsec_plus = MakeXSecCard(lepname+"plus", is5TeV, outdir = outdir_card, applyLFU=applyLFU)
    card_xsec_minus = MakeXSecCard(lepname+"minus", is5TeV, outdir = outdir_card, applyLFU=applyLFU)

    card_z = None
    if fzsig_input:
        # rebin the z histograms
        # for Z's should not include the overflow or underflows, since there is the mll cut
        ProcessHists(fzsig_input, fzsig_output, mass_bins_z, False, False)
        # generate z card
        card_z = MakeZJetsCards(fzsig_output, lepname+lepname, rebinned = True, is5TeV = is5TeV, outdir = outdir_card, applyLFU=applyLFU)
        card_xsec_z = MakeXSecCard(lepname+lepname, is5TeV, outdir = outdir_card, applyLFU=applyLFU)

    return card_plus, card_minus, card_z, card_xsec_plus, card_xsec_minus, card_xsec_z

if __name__  == "__main__":
    applyLFU = True

    # scan fit range
    mass_bins_test = OrderedDict()
    mass_bins_test[0] = np.array([0., 20, 30, 40, 50, 60, 70, 80, 90, 120])
    mass_bins_test[1] = np.array([5., 25, 35, 45, 55, 65, 75, 85, 95, 120])
    mass_bins_test[2] = np.array([10., 20, 30, 40, 50, 60, 70, 80, 90, 120])
    mass_bins_test[3] = np.array([15., 25, 35, 45, 55, 65, 75, 85, 95, 120])
    mass_bins_test[4] = np.array([20., 30, 40, 50, 60, 70, 80, 90, 120])
    mass_bins_test[5] = np.array([25., 35, 45, 55, 65, 75, 85, 95, 120])
    mass_bins_test[6] = np.array([30., 40, 50, 60, 70, 80, 90, 100, 120])
    mass_bins_test[7] = np.array([35., 45, 55, 65, 75, 85, 95, 120])
    mass_bins_test[8] = np.array([40., 50, 60, 70, 80, 90, 120])
    mass_bins_test[9] = np.array([45., 55, 65, 75, 85, 95, 120])
    mass_bins_test[10] = np.array([50., 60, 70, 80, 90, 120])

    mass_bin_z = np.array([60, 70, 80, 84, 88, 92, 96, 100, 104, 112, 120])

    for sqrtS in ["13TeV", "5TeV"]:
        card_muplus = None
        card_muminus = None
        card_zmumu = None
        card_eplus = None
        card_eminus = None
        card_zee = None

        suffix = ".root" if sqrtS == "13TeV" else "_5TeV.root"
        do5TeV = (sqrtS == "5TeV")

        fzsig_mumu_input = "root/output_shapes_mumu" + suffix
        fzsig_ee_input = "root/output_shapes_ee" + suffix

        for key, val in mass_bins_test.items():
            #if key != 4:
            #    continue
            # muon channel
            fwsig_input = "root/output_shapes_munu" + suffix
            fwsig_rebin = f"root/test/{sqrtS}/test{key}/output_shapes_munu_Rebin" + suffix
            fwsig_mergeTau = f"root/test/{sqrtS}/test{key}/output_shapes_munu_mergeTau" + suffix

            fqcd_input = "root/output_qcdshape_munu" + suffix
            fqcd_rebin = f"root/test/{sqrtS}/test{key}/output_qcdshape_munu_Rebin" + suffix
            fqcd_input_scaled = "root/output_qcdshape_munu_applyScaling" + suffix
            fqcd_rebin_scaled = f"root/test/{sqrtS}/test{key}/output_qcdshape_munu_Rebin_applyScaling" + suffix
            fqcd_output = f"root/test/{sqrtS}/test{key}/qcdshape_extrapolated_munu" + suffix

            fzsig_output = f"root/test/{sqrtS}/test{key}/output_shapes_mumu_Rebin" + suffix

            card_muplus, card_muminus, card_zmumu, card_xsec_muplus, card_xsec_muminus, card_xsec_mumu = RunPreparations(fwsig_input, fwsig_rebin, fwsig_mergeTau, fqcd_input, fqcd_rebin, fqcd_input_scaled, fqcd_rebin_scaled, fqcd_output, "mu", outdir_card = f"cards/test/{sqrtS}/test{key}", mass_bins_w = val, fzsig_input=fzsig_mumu_input, fzsig_output = fzsig_output, applyLFU=applyLFU, is5TeV=do5TeV)

            # electron channel
            fwsig_input = "root/output_shapes_enu" + suffix
            fwsig_rebin = f"root/test/{sqrtS}/test{key}/output_shapes_enu_Rebin" + suffix
            fwsig_mergeTau = f"root/test/{sqrtS}/test{key}/output_shapes_enu_mergeTau" + suffix

            fqcd_input = "root/output_qcdshape_enu" + suffix
            fqcd_rebin = f"root/test/{sqrtS}/test{key}/output_qcdshape_enu_Rebin" + suffix
            fqcd_input_scaled = "root/output_qcdshape_enu_applyScaling" + suffix
            fqcd_rebin_scaled = f"root/test/{sqrtS}/test{key}/output_qcdshape_enu_Rebin_applyScaling" + suffix
            fqcd_output = f"root/test/{sqrtS}/test{key}/qcdshape_extrapolated_enu" + suffix

            fzsig_output = f"root/test/{sqrtS}/test{key}/output_shapes_ee_Rebin" + suffix

            card_eplus, card_eminus, card_zee, card_xsec_eplus, card_xsec_eminus, card_xsec_ee = RunPreparations(fwsig_input, fwsig_rebin, fwsig_mergeTau, fqcd_input, fqcd_rebin, fqcd_input_scaled, fqcd_rebin_scaled, fqcd_output, "e", outdir_card = f"cards/test/{sqrtS}/test{key}", mass_bins_w = val, fzsig_input=fzsig_ee_input, fzsig_output=fzsig_output, applyLFU=applyLFU, is5TeV=do5TeV)

            # combined fit
            GenerateRunCommand("card_combined", [card_muplus, card_muminus, card_zmumu, card_eplus, card_eminus, card_zee], ["muplus", "muminus", "mumu", "eplus", "eminus", "ee"], [card_xsec_muplus, card_xsec_muminus, card_xsec_mumu, card_xsec_eplus, card_xsec_eminus, card_xsec_ee], applyLFU=applyLFU)

            # mu channel fit
            GenerateRunCommand("card_mu", [card_muplus, card_muminus, card_zmumu], ["muplus", "muminus", "mumu"], [card_xsec_muplus, card_xsec_muminus, card_xsec_mumu], applyLFU=applyLFU)

            # e channel fit
            GenerateRunCommand("card_e", [card_eplus, card_eminus, card_zee], ["eplus", "eminus", "ee"], [card_xsec_eplus, card_xsec_eminus, card_xsec_ee], applyLFU=applyLFU)


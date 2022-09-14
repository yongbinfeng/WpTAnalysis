import ROOT
import re
import numpy as np
from CMSPLOTS.myFunction import AddOverflowsTH1, RebinHisto
from modules.qcdExtrapolater import ExtrapolateQCD
from modules.cardMaker import MakeCards
from modules.histProcessor import ProcessHists, CopyandMergeTau
from modules.mass_bins import mass_bins

ROOT.gROOT.SetBatch(True)


def RunPreparations(fsig_input, fsig_rebin, fsig_mergeTau, fqcd_input, fqcd_rebin, fqcd_input_scaled, fqcd_rebin_scaled, fqcd_output, lepname = "mu", is5TeV = False):
    """
    rebin the mT hists for sig and qcd bkg,
    merge tau (wx) process to sig (wlnu),
    exptrapolate the qcd shape from anti-isolated to isolated
    genereate datacard
    """
    includeUnderflow = False
    includeOverflow = True

    ProcessHists(fsig_input, fsig_rebin, mass_bins, includeUnderflow, includeOverflow)
    # merge the tau processes into the signal process
    CopyandMergeTau(fsig_rebin, fsig_mergeTau)

    # for QCD
    ProcessHists(fqcd_input, fqcd_rebin, mass_bins, includeUnderflow, includeOverflow)
    ProcessHists(fqcd_input_scaled, fqcd_rebin_scaled, mass_bins, includeUnderflow, includeOverflow)

    # extrapolate the QCD template from anti-isolated region to isolated region
    ExtrapolateQCD(fqcd_rebin, fqcd_output, [lepname+"plus", lepname+"minus"], "WpT_bin0", ["lepEta_bin0"], fname_scaled=fqcd_rebin_scaled, rebinned = True, is5TeV = is5TeV)

    # generate card based on the signal and qcd templates
    MakeCards(fsig_mergeTau, fqcd_output, lepname+"plus", "WpT_bin0", "lepEta_bin0", rebinned = True, nMTBins = len(mass_bins)-1, is5TeV = is5TeV)
    MakeCards(fsig_mergeTau, fqcd_output, lepname+"minus", "WpT_bin0", "lepEta_bin0", rebinned = True, nMTBins = len(mass_bins)-1, is5TeV = is5TeV)


if "__main__" == __name__:
    do13TeV = False
    do5TeV = True

    if do13TeV:
        # 13TeV muon
        fsig_input = "root/output_shapes_munu.root"
        fsig_rebin = "root/output_shapes_munu_Rebin.root"
        fsig_mergeTau = "root/output_shapes_munu_mergeTau.root"

        fqcd_input = "root/output_qcdshape_munu.root"
        fqcd_rebin = "root/output_qcdshape_munu_Rebin.root"
        fqcd_input_scaled = "root/output_qcdshape_munu_applyScaling.root"
        fqcd_rebin_scaled = "root/output_qcdshape_munu_Rebin_applyScaling.root"
        fqcd_output = "root/qcdshape_extrapolated_munu.root"

        RunPreparations(fsig_input, fsig_rebin, fsig_mergeTau, fqcd_input, fqcd_rebin, fqcd_input_scaled, fqcd_rebin_scaled, fqcd_output, "mu")

        # 13TeV electron
        fsig_input = "root/output_shapes_enu.root"
        fsig_rebin = "root/output_shapes_enu_Rebin.root"
        fsig_mergeTau = "root/output_shapes_enu_mergeTau.root"

        fqcd_input = "root/output_qcdshape_enu.root"
        fqcd_rebin = "root/output_qcdshape_enu_Rebin.root"
        fqcd_input_scaled = "root/output_qcdshape_enu_applyScaling.root"
        fqcd_rebin_scaled = "root/output_qcdshape_enu_Rebin_applyScaling.root"
        fqcd_output = "root/qcdshape_extrapolated_enu.root"

        RunPreparations(fsig_input, fsig_rebin, fsig_mergeTau, fqcd_input, fqcd_rebin, fqcd_input_scaled, fqcd_rebin_scaled, fqcd_output, "e")

    if do5TeV:
        # 5TeV muon channel
        fsig_input = "root/output_shapes_munu_5TeV.root"
        fsig_rebin = "root/output_shapes_munu_Rebin_5TeV.root"
        fsig_mergeTau = "root/output_shapes_munu_mergeTau_5TeV.root"

        fqcd_input = "root/output_qcdshape_munu_5TeV.root"
        fqcd_rebin = "root/output_qcdshape_munu_Rebin_5TeV.root"
        fqcd_input_scaled = "root/output_qcdshape_munu_applyScaling_5TeV.root"
        fqcd_rebin_scaled = "root/output_qcdshape_munu_Rebin_applyScaling_5TeV.root"
        fqcd_output = "root/qcdshape_extrapolated_munu_5TeV.root"

        RunPreparations(fsig_input, fsig_rebin, fsig_mergeTau, fqcd_input, fqcd_rebin, fqcd_input_scaled, fqcd_rebin_scaled, fqcd_output, "mu", is5TeV = True)


        # 5TeV electron channel
        fsig_input = "root/output_shapes_enu_5TeV.root"
        fsig_rebin = "root/output_shapes_enu_Rebin_5TeV.root"
        fsig_mergeTau = "root/output_shapes_enu_mergeTau_5TeV.root"

        fqcd_input = "root/output_qcdshape_enu_5TeV.root"
        fqcd_rebin = "root/output_qcdshape_enu_Rebin_5TeV.root"
        fqcd_input_scaled = "root/output_qcdshape_enu_applyScaling_5TeV.root"
        fqcd_rebin_scaled = "root/output_qcdshape_enu_Rebin_applyScaling_5TeV.root"
        fqcd_output = "root/qcdshape_extrapolated_enu_5TeV.root"

        RunPreparations(fsig_input, fsig_rebin, fsig_mergeTau, fqcd_input, fqcd_rebin, fqcd_input_scaled, fqcd_rebin_scaled, fqcd_output, "e", is5TeV = True)


"""
script to make postfit comparisons
"""
from modules.postFitScripts import MakePostPlot, result2json, MakeWpTPostFitPlots, GetPOIValue, ComparePOIs, DumpGroupImpacts
from modules.CombineHarvester.plotImpacts import plotImpacts
from modules.mass_bins import mass_bins
from modules.Utils import FormatTable
import ROOT
import numpy as np
from collections import OrderedDict
ROOT.gROOT.SetBatch(True)

doInclusive = True
doWpT = False
doRebin = False

# boolean flag to config if the pull
# distribution should be included in the plots
showPULL = True
if doInclusive:
    mass_bins_test = OrderedDict()
    mass_bins_test[0] = np.array([0., 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 120])
    mass_bins_test[1] = np.array([5., 15, 25, 35, 45, 55, 65, 75, 85, 95, 105, 120])
    mass_bins_test[2] = np.array([10., 20, 30, 40, 50, 60, 70, 80, 90, 100, 120])
    mass_bins_test[3] = np.array([15., 25, 35, 45, 55, 65, 75, 85, 95, 105, 120])
    mass_bins_test[4] = np.array([20., 30, 40, 50, 60, 70, 80, 90, 100, 120])
    mass_bins_test[5] = np.array([25., 35, 45, 55, 65, 75, 85, 95, 105, 120])
    mass_bins_test[6] = np.array([30., 40, 50, 60, 70, 80, 90, 100, 120])
    mass_bins_test[7] = np.array([35., 45, 55, 65, 75, 85, 95, 105, 120])
    mass_bins_test[8] = np.array([40., 50, 60, 70, 80, 90, 100, 120])
    mass_bins_test[9] = np.array([45., 55, 65, 75, 85, 95, 105, 120])
    mass_bins_test[10] = np.array([50., 60, 70, 80, 90, 100, 120])
    z_mass_bins = np.arange(60.0, 124.0, 4)

    vals_x = np.array([0, 5, 10, 15, 20, 25, 30, 35, 40, 45.0, 50.0], dtype='d')
    ntests = len(vals_x)
    vals_lep_pos = np.zeros(ntests)
    vals_lep_neg = np.zeros(ntests)
    errs_lep_pos = np.zeros(ntests)
    errs_lep_neg = np.zeros(ntests)

    vals_mu_pos = np.zeros(ntests)
    errs_mu_pos = np.zeros(ntests)
    vals_mu_neg = np.zeros(ntests)
    errs_mu_neg = np.zeros(ntests)

    vals_e_pos = np.zeros(ntests)
    errs_e_pos = np.zeros(ntests)
    vals_e_neg = np.zeros(ntests)
    errs_e_neg = np.zeros(ntests)

    lumi_unc = 0.017

    #for idx in range(ntests):
    #    if idx != 0:
    #        continue
    for idx in range(4,5):
        sqrtS = "5TeV"
        do5TeV = (sqrtS == "5TeV")
        mass_bins = mass_bins_test[idx]

        #filename = f"cards/test{idx}/datacard_mucombined_lepEta_bin0_WpT_bin0.root"
        filename = f"cards/{sqrtS}/test{idx}/card_combined.root"
        filename_e = f"cards/{sqrtS}/test{idx}/card_combined_e.root"
        filename_mu = f"cards/{sqrtS}/test{idx}/card_combined_mu.root"

        vals_lep_pos[idx], errs_lep_pos[idx] = GetPOIValue(filename, "lepplus_sig_mu")
        vals_lep_neg[idx], errs_lep_neg[idx] = GetPOIValue(filename, "lepminus_sig_mu")
        #vals_mu_pos[idx], errs_mu_pos[idx] = GetPOIValue(filename_mu, "w_lepplus_sig_mu")
        #vals_mu_neg[idx], errs_mu_neg[idx] = GetPOIValue(filename_mu, "w_lepminus_sig_mu")

        nevts = OrderedDict()

        nevts['muplus']  = MakePostPlot(filename, "muplus",  mass_bins, f"_mT{idx}_{sqrtS}", showPULL, is5TeV=do5TeV)
        nevts['muminus'] = MakePostPlot(filename, "muminus", mass_bins, f"_mT{idx}_{sqrtS}", showPULL, startbin = len(mass_bins), is5TeV=do5TeV)
        nevts['mumu']    = MakePostPlot(filename, "mumu", z_mass_bins,  f"_mT{idx}_{sqrtS}", showPULL, startbin = len(mass_bins)*2 - 1, is5TeV=do5TeV)
        nevts['eplus']   = MakePostPlot(filename, "eplus",  mass_bins,  f"_mT{idx}_{sqrtS}", showPULL, startbin = len(mass_bins)*2 + len(z_mass_bins) - 2, is5TeV=do5TeV)
        nevts['eminus']  = MakePostPlot(filename, "eminus", mass_bins,  f"_mT{idx}_{sqrtS}", showPULL, startbin = len(mass_bins)*3 + len(z_mass_bins) - 3, is5TeV=do5TeV)
        nevts['ee']      = MakePostPlot(filename, "ee",   z_mass_bins,  f"_mT{idx}_{sqrtS}", showPULL, startbin = len(mass_bins)*4 + len(z_mass_bins) - 4, is5TeV=do5TeV)

        result2json(filename, "lepplus_sig_mu",  f"cards/{sqrtS}/test{idx}/impacts_lepplus.json")
        result2json(filename, "lepminus_sig_mu", f"cards/{sqrtS}/test{idx}/impacts_lepminus.json")
        result2json(filename, "leplep_sig_mu",   f"cards/{sqrtS}/test{idx}/impacts_leplep.json")

        result2json(filename, "WZRatio_ratiometaratio",   f"cards/{sqrtS}/test{idx}/impacts_WZRatio.json",   "nuisance_impact_ratiometapois")
        result2json(filename, "WchgRatio_ratiometaratio", f"cards/{sqrtS}/test{idx}/impacts_WchgRatio.json", "nuisance_impact_ratiometapois")
        result2json(filename, "WchgAsym_chargemetaasym",  f"cards/{sqrtS}/test{idx}/impacts_WchgAsym.json",  "nuisance_impact_chargemetapois")

        plotImpacts(f"cards/{sqrtS}/test{idx}/impacts_lepplus.json",  f"impacts_lepplus_mT{idx}_{sqrtS}")
        plotImpacts(f"cards/{sqrtS}/test{idx}/impacts_lepminus.json", f"impacts_lepminus_mT{idx}_{sqrtS}")
        plotImpacts(f"cards/{sqrtS}/test{idx}/impacts_leplep.json",   f"impacts_leplep_mll_mT{idx}_{sqrtS}")

        plotImpacts(f"cards/{sqrtS}/test{idx}/impacts_WZRatio.json",   f"impacts_WZRatio_mT{idx}_{sqrtS}")
        plotImpacts(f"cards/{sqrtS}/test{idx}/impacts_WchgRatio.json", f"impacts_WchgRatio_mT{idx}_{sqrtS}")
        plotImpacts(f"cards/{sqrtS}/test{idx}/impacts_WchgAsym.json",  f"impacts_WchgAsym_mT{idx}_{sqrtS}")

        #vals_e_pos[idx], errs_e_pos[idx] = GetPOIValue(filename, "w_eplus_sig_mu")
        #vals_e_neg[idx], errs_e_neg[idx] = GetPOIValue(filename, "w_eminus_sig_mu")

        #MakePostPlot(filename, "eplus",  mass_bins, f"_mT{idx}", showPULL)
        #MakePostPlot(filename, "eminus", mass_bins, f"_mT{idx}", showPULL, startbin = len(mass_bins))

        #result2json(filename, "w_eplus_sig_mu", f"cards/test{idx}/impacts_eplus.json")
        #result2json(filename, "w_eminus_sig_mu",  f"cards/test{idx}/impacts_eminus.json")

        #plotImpacts(f"cards/test{idx}/impacts_eplus.json",  f"impacts_eplus_mT{idx}")
        #plotImpacts(f"cards/test{idx}/impacts_eminus.json", f"impacts_eminus_mT{idx}")


        impacts = OrderedDict()

        impacts['lepplus']  = DumpGroupImpacts(filename, "lepplus_sig_mu")
        impacts['lepminus'] = DumpGroupImpacts(filename, "lepminus_sig_mu")
        impacts['leplep']   = DumpGroupImpacts(filename, "leplep_sig_mu")

        impacts['WOverZ']   = DumpGroupImpacts(filename, "WZRatio_ratiometaratio",   "nuisance_group_impact_ratiometapois")
        impacts['WAsym']    = DumpGroupImpacts(filename, "WchgAsym_chargemetaasym",  "nuisance_group_impact_chargemetapois")


    #errs_mu_pos = np.sqrt(errs_mu_pos**2 - lumi_unc**2)
    #errs_mu_neg = np.sqrt(errs_mu_neg**2 - lumi_unc**2)
    #errs_e_pos = np.sqrt(errs_e_pos**2 - lumi_unc**2)
    #errs_e_neg = np.sqrt(errs_e_neg**2 - lumi_unc**2)
    errs_lep_pos = np.sqrt(errs_lep_pos**2 - lumi_unc**2)
    errs_lep_neg = np.sqrt(errs_lep_neg**2 - lumi_unc**2)

    labels = [
        "W^{+} #rightarrow #mu^{+}#nu",
        "W^{-} #rightarrow #mu^{-}#bar{#nu}",
        "W^{+} #rightarrow e^{+}#nu",
        "W^{-} #rightarrow e^{-}#bar{#nu}",
        "W^{+} #rightarrow #ell^{+}#nu",
        "W^{-} #rightarrow #ell^{-}#bar{#nu}",
    ]
    markers = [22, 23, 26, 32]
    colors = [2, 4, 6, 8]
    #ComparePOIs(vals_x, [vals_mu_pos, vals_mu_neg], [errs_mu_pos, errs_mu_neg], labels[:2], colors[:2], markers[:2], "test_mu")
    #ComparePOIs(vals_x, [vals_e_pos, vals_e_neg], [errs_e_pos, errs_e_neg], labels[2:], colors[2:], markers[2:], "test_e")
    #ComparePOIs(vals_x, [vals_mu_pos, vals_mu_neg, vals_e_pos, vals_e_neg], [errs_mu_pos, errs_mu_neg, errs_e_pos, errs_e_neg], labels, colors, markers, "test")
    #ComparePOIs(vals_x, [vals_lep_pos, vals_lep_neg], [errs_lep_pos, errs_lep_neg], labels[4:], colors[:2], markers[:2], "test_lep")

    # plot the impacts on QCD
    #result2json("cards/datacard_muplus_lepEta_bin0_WpT_bin0", "QCD_muplus_lepEta_bin0_WpT_bin0_mu", "cards/impacts_muplus_lepEta_bin0_WpT_bin0_qcd.json")
    #result2json("cards/datacard_eplus",       "QCD_eplus_lepEta_bin1_WpT_bin0_mu",  "cards/impacts_eplus_qcd.json")

    #plotImpacts("cards/impacts_muplus_lepEta_bin0_WpT_bin0_qcd.json", "impacts_muplus_lepEta_bin0_qcd")
    #plotImpacts("cards/impacts_eplus_qcd.json",       "impacts_eplus_lepEta_qcd")

    ## print out the nevts information
    #print(FormatTable(nevts, ['muplus', 'muminus', 'mumu']))
    print(FormatTable(nevts, caption="Event yield at 13 TeV", label = f"tab:nevts_{sqrtS}"))
    print("\n\n\n\n")
    print(FormatTable(impacts, caption=f"Systematic uncertainties in percentage at {sqrtS}", label = f"tab:impacts_{sqrtS}", precision=2))

if doWpT:
    MakePostPlot("cards/datacard_muplus_lepEta_bin0_WpT.root", "muplus", "WpT", showPULL)

    wpttruthbins  = ["WpT_truth_bin1", "WpT_truth_bin2", "WpT_truth_bin3", "WpT_truth_bin4", "WpT_truth_bin5", "WpT_truth_bin6", "WpT_truth_bin7", "WpT_truth_bin8", "WpT_truth_bin9"]
   
    jsonNames = []
    for wpttruthbin in wpttruthbins:
        result2json("cards/datacard_muplus_lepEta_bin0_WpT.root", "w_muplus_{}_sig_mu".format(wpttruthbin), "cards/impacts_muplus_{}_sig_mu.json".format(wpttruthbin))
        jsonNames.append("cards/impacts_muplus_{}_sig_mu.json".format(wpttruthbin))
        plotImpacts("cards/impacts_muplus_{}_sig_mu.json".format(wpttruthbin), "impacts_muplus_{}".format(wpttruthbin))

    MakeWpTPostFitPlots(jsonNames)


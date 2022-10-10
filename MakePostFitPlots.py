"""
script to make postfit comparisons
"""
from modules.postFitScripts import MakePostPlot, result2json, MakeWpTPostFitPlots, GetPOIValue, ComparePOIs, DumpGroupImpacts
from modules.CombineHarvester.plotImpacts import plotImpacts
from modules.mass_bins import mass_bins
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
    MakePostPlot("cards/test_mu_withXsecSys.root", "muplus", mass_bins, "_withXsecSys", showPULL)
    MakePostPlot("cards/test_mu_withXsecSys.root", "muminus", mass_bins, "_withXsecSys", showPULL, startbin = len(mass_bins)-1)

    result2json("cards/test_mu_withXsecSys.root", "w_muminus_sig_mu", "cards/impacts_muminus_lepEta_bin0_WpT_bin0_withXsecSys.json")
    result2json("cards/test_mu_withXsecSys.root", "w_muplus_sig_mu", "cards/impacts_muplus_lepEta_bin0_WpT_bin0_withXsecSys.json")

    plotImpacts("cards/impacts_muplus_lepEta_bin0_WpT_bin0_withXsecSys.json", "impacts_muplus_lepEta_bin0_withXsecSys")
    plotImpacts("cards/impacts_muminus_lepEta_bin0_WpT_bin0_withXsecSys.json", "impacts_muminus_lepEta_bin0_withXsecSys")

    DumpGroupImpacts("cards/test_mu_withXsecSys.root", "w_muminus_sig_mu")

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

    vals_x = np.array([0, 5, 10, 15, 20, 25, 30, 35, 40, 45.0, 50.0], dtype='d')
    ntests = len(vals_x)
    vals_mu_pos = np.zeros(ntests)
    errs_mu_pos = np.zeros(ntests)
    vals_mu_neg = np.zeros(ntests)
    errs_mu_neg = np.zeros(ntests)
    vals_e_pos = np.zeros(ntests)
    errs_e_pos = np.zeros(ntests)
    vals_e_neg = np.zeros(ntests)
    errs_e_neg = np.zeros(ntests)

    lumi_unc = 0.017

    for idx in range(ntests):
        mass_bins = mass_bins_test[idx]
        filename = f"cards/test{idx}/datacard_mucombined_lepEta_bin0_WpT_bin0.root"
        vals_mu_pos[idx], errs_mu_pos[idx] = GetPOIValue(filename, "w_muplus_sig_mu")
        vals_mu_neg[idx], errs_mu_neg[idx] = GetPOIValue(filename, "w_muminus_sig_mu")

        MakePostPlot(filename, "muplus",  mass_bins, f"_mT{idx}", showPULL)
        MakePostPlot(filename, "muminus", mass_bins, f"_mT{idx}", showPULL, startbin = len(mass_bins)-1)

        result2json(filename, "w_muplus_sig_mu", f"cards/test{idx}/impacts_muplus.json")
        result2json(filename, "w_muminus_sig_mu",  f"cards/test{idx}/impacts_muminus.json")

        plotImpacts(f"cards/test{idx}/impacts_muplus.json",  f"impacts_muplus_mT{idx}")
        plotImpacts(f"cards/test{idx}/impacts_muminus.json", f"impacts_muminus_mT{idx}")

        filename = f"cards/test{idx}/datacard_ecombined_lepEta_bin0_WpT_bin0.root"
        vals_e_pos[idx], errs_e_pos[idx] = GetPOIValue(filename, "w_eplus_sig_mu")
        vals_e_neg[idx], errs_e_neg[idx] = GetPOIValue(filename, "w_eminus_sig_mu")

        MakePostPlot(filename, "eplus",  mass_bins, f"_mT{idx}", showPULL)
        MakePostPlot(filename, "eminus", mass_bins, f"_mT{idx}", showPULL, startbin = len(mass_bins)-1)

        result2json(filename, "w_eplus_sig_mu", f"cards/test{idx}/impacts_eplus.json")
        result2json(filename, "w_eminus_sig_mu",  f"cards/test{idx}/impacts_eminus.json")

        plotImpacts(f"cards/test{idx}/impacts_eplus.json",  f"impacts_eplus_mT{idx}")
        plotImpacts(f"cards/test{idx}/impacts_eminus.json", f"impacts_eminus_mT{idx}")

        DumpGroupImpacts(filename, "w_eplus_sig_mu")

    errs_mu_pos = np.sqrt(errs_mu_pos**2 - lumi_unc**2)
    errs_mu_neg = np.sqrt(errs_mu_neg**2 - lumi_unc**2)
    errs_e_pos = np.sqrt(errs_e_pos**2 - lumi_unc**2)
    errs_e_neg = np.sqrt(errs_e_neg**2 - lumi_unc**2)

    labels = [
        "W^{+} #rightarrow #mu^{+}#nu",
        "W^{-} #rightarrow #mu^{-}#bar{#nu}",
        "W^{+} #rightarrow e^{+}#nu",
        "W^{-} #rightarrow e^{-}#bar{#nu}",
    ]
    markers = [22, 23, 26, 32]
    colors = [2, 4, 6, 8]
    ComparePOIs(vals_x, [vals_mu_pos, vals_mu_neg], [errs_mu_pos, errs_mu_neg], labels[:2], colors[:2], markers[:2], "test_mu")
    ComparePOIs(vals_x, [vals_e_pos, vals_e_neg], [errs_e_pos, errs_e_neg], labels[2:], colors[2:], markers[2:], "test_e")
    ComparePOIs(vals_x, [vals_mu_pos, vals_mu_neg, vals_e_pos, vals_e_neg], [errs_mu_pos, errs_mu_neg, errs_e_pos, errs_e_neg], labels, colors, markers, "test")

    #MakePostPlot("cards/test_mu.root", "muplus", mass_bins, "", showPULL)
    #MakePostPlot("cards/test_mu.root", "muminus", mass_bins, "", showPULL, startbin = len(mass_bins)-1)
    #MakePostPlot("cards/test_e.root",  "eplus", mass_bins, "", showPULL)
    #MakePostPlot("cards/test_e.root",  "eminus", mass_bins, "", showPULL, startbin = len(mass_bins)-1)
    #MakePostPlot("cards/test_5TeV.root", "muplus", mass_bins, "_5TeV", showPULL, is5TeV=True)
    #MakePostPlot("cards/test_5TeV.root", "muminus", mass_bins, "_5TeV", showPULL, is5TeV=True, startbin = len(mass_bins)-1)
    #MakePostPlot("cards/test_e_5TeV.root", "eplus", mass_bins, "_5TeV", showPULL, is5TeV=True)
    #MakePostPlot("cards/test_e_5TeV.root", "eminus", mass_bins, "_5TeV", showPULL, is5TeV=True, startbin = len(mass_bins)-1)
    ##
    #result2json("cards/test_mu.root", "w_muminus_sig_mu", "cards/impacts_muminus_lepEta_bin0_WpT_bin0.json")
    #result2json("cards/test_mu.root", "w_muplus_sig_mu", "cards/impacts_muplus_lepEta_bin0_WpT_bin0.json")
    #result2json("cards/test_e.root",  "w_eplus_sig_mu",  "cards/impacts_eplus_lepEta_bin0_WpT_bin0.json")
    #result2json("cards/test_e.root",  "w_eminus_sig_mu",  "cards/impacts_eminus_lepEta_bin0_WpT_bin0.json")

    #result2json("cards/test_5TeV.root", "w_muminus_sig_mu", "cards/impacts_muminus_lepEta_bin0_WpT_bin0_5TeV.json")
    #result2json("cards/test_5TeV.root", "w_muplus_sig_mu", "cards/impacts_muplus_lepEta_bin0_WpT_bin0_5TeV.json")
    #result2json("cards/test_e_5TeV.root", "w_eminus_sig_mu", "cards/impacts_eminus_lepEta_bin0_WpT_bin0_5TeV.json")
    #result2json("cards/test_e_5TeV.root", "w_eplus_sig_mu",  "cards/impacts_eplus_lepEta_bin0_WpT_bin0_5TeV.json")

    #plotImpacts("cards/impacts_muplus_lepEta_bin0_WpT_bin0.json", "impacts_muplus_lepEta_bin0")
    #plotImpacts("cards/impacts_muminus_lepEta_bin0_WpT_bin0.json", "impacts_muminus_lepEta_bin0")
    #plotImpacts("cards/impacts_eplus_lepEta_bin0_WpT_bin0.json",       "impacts_eplus_lepEta_bin0")
    #plotImpacts("cards/impacts_eminus_lepEta_bin0_WpT_bin0.json",       "impacts_eminus_lepEta_bin0")

    #plotImpacts("cards/impacts_muplus_lepEta_bin0_WpT_bin0_5TeV.json", "impacts_muplus_lepEta_bin0_5TeV")
    #plotImpacts("cards/impacts_muminus_lepEta_bin0_WpT_bin0_5TeV.json", "impacts_muminus_lepEta_bin0_5TeV")
    #plotImpacts("cards/impacts_eplus_lepEta_bin0_WpT_bin0_5TeV.json", "impacts_eplus_lepEta_bin0_5TeV")
    #plotImpacts("cards/impacts_eminus_lepEta_bin0_WpT_bin0_5TeV.json", "impacts_eminus_lepEta_bin0_5TeV")

    # plot the impacts on QCD
    #result2json("cards/datacard_muplus_lepEta_bin0_WpT_bin0", "QCD_muplus_lepEta_bin0_WpT_bin0_mu", "cards/impacts_muplus_lepEta_bin0_WpT_bin0_qcd.json")
    #result2json("cards/datacard_eplus",       "QCD_eplus_lepEta_bin1_WpT_bin0_mu",  "cards/impacts_eplus_qcd.json")

    #plotImpacts("cards/impacts_muplus_lepEta_bin0_WpT_bin0_qcd.json", "impacts_muplus_lepEta_bin0_qcd")
    #plotImpacts("cards/impacts_eplus_qcd.json",       "impacts_eplus_lepEta_qcd")

if doWpT:
    MakePostPlot("cards/datacard_muplus_lepEta_bin0_WpT.root", "muplus", "WpT", showPULL)

    wpttruthbins  = ["WpT_truth_bin1", "WpT_truth_bin2", "WpT_truth_bin3", "WpT_truth_bin4", "WpT_truth_bin5", "WpT_truth_bin6", "WpT_truth_bin7", "WpT_truth_bin8", "WpT_truth_bin9"]
   
    jsonNames = []
    for wpttruthbin in wpttruthbins:
        result2json("cards/datacard_muplus_lepEta_bin0_WpT.root", "w_muplus_{}_sig_mu".format(wpttruthbin), "cards/impacts_muplus_{}_sig_mu.json".format(wpttruthbin))
        jsonNames.append("cards/impacts_muplus_{}_sig_mu.json".format(wpttruthbin))
        plotImpacts("cards/impacts_muplus_{}_sig_mu.json".format(wpttruthbin), "impacts_muplus_{}".format(wpttruthbin))

    MakeWpTPostFitPlots(jsonNames)


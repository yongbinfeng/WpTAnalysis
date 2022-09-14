"""
script to make postfit comparisons
"""
from modules.postFitScripts import MakePostPlot, result2json, MakeWpTPostFitPlots
from modules.CombineHarvester.plotImpacts import plotImpacts
from modules.mass_bins import mass_bins
import ROOT
ROOT.gROOT.SetBatch(True)

doInclusive = True
doWpT = False
doRebin = False

# boolean flag to config if the pull
# distribution should be included in the plots
showPULL = True
if doInclusive:
    MakePostPlot("cards/datacard_muminus_lepEta_bin0_WpT_bin0.root", "muminus", mass_bins, "", showPULL)
    MakePostPlot("cards/datacard_muplus_lepEta_bin0_WpT_bin0.root", "muplus", mass_bins, "", showPULL)
    MakePostPlot("cards/test_e.root",  "eplus", mass_bins, "", showPULL)
    MakePostPlot("cards/test_5TeV.root", "muplus", mass_bins, "_5TeV", showPULL, is5TeV=True)
    MakePostPlot("cards/test_5TeV.root", "muminus", mass_bins, "_5TeV", showPULL, is5TeV=True, startbin = len(mass_bins)-1)
    MakePostPlot("cards/test_e_5TeV.root", "eplus", mass_bins, "_5TeV", showPULL, is5TeV=True)
    MakePostPlot("cards/test_e_5TeV.root", "eminus", mass_bins, "_5TeV", showPULL, is5TeV=True, startbin = len(mass_bins)-1)
    #
    #result2json("cards/datacard_mu_lepEta_bin0_WpT_bin0.root", "w_muminus_sig_mu", "cards/impacts_muminus_lepEta_bin0_WpT_bin0.json")
    #result2json("cards/datacard_mu_lepEta_bin0_WpT_bin0.root", "w_muplus_sig_mu", "cards/impacts_muplus_lepEta_bin0_WpT_bin0.json")
    #result2json("cards/test_e.root",       "w_eplus_sig_mu",  "cards/impacts_eplus.json")
    #result2json("cards/test_e.root",       "w_eminus_sig_mu",  "cards/impacts_eminus.json")

    #result2json("cards/test_5TeV.root", "w_muminus_sig_mu", "cards/impacts_muminus_lepEta_bin0_WpT_bin0_5TeV.json")
    #result2json("cards/test_5TeV.root", "w_muplus_sig_mu", "cards/impacts_muplus_lepEta_bin0_WpT_bin0_5TeV.json")

    #result2json("cards/test_e_5TeV.root", "w_eminus_sig_mu", "cards/impacts_eminus_lepEta_bin0_WpT_bin0_5TeV.json")
    #result2json("cards/test_e_5TeV.root", "w_eplus_sig_mu",  "cards/impacts_eplus_lepEta_bin0_WpT_bin0_5TeV.json")

    #plotImpacts("cards/impacts_muminus_lepEta_bin0_WpT_bin0.json", "impacts_muminus_lepEta_bin0")
    #plotImpacts("cards/impacts_muplus_lepEta_bin0_WpT_bin0.json", "impacts_muplus_lepEta_bin0")
    #plotImpacts("cards/impacts_eplus.json",       "impacts_eplus_lepEta")
    #plotImpacts("cards/impacts_eminus.json",       "impacts_eminus_lepEta")

    #plotImpacts("cards/impacts_muminus_lepEta_bin0_WpT_bin0_5TeV.json", "impacts_muminus_lepEta_bin0_5TeV")
    #plotImpacts("cards/impacts_muplus_lepEta_bin0_WpT_bin0_5TeV.json", "impacts_muplus_lepEta_bin0_5TeV")
    #
    #plotImpacts("cards/impacts_eminus_lepEta_bin0_WpT_bin0_5TeV.json", "impacts_eminus_lepEta_bin0_5TeV")
    #plotImpacts("cards/impacts_eplus_lepEta_bin0_WpT_bin0_5TeV.json", "impacts_eplus_lepEta_bin0_5TeV")


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


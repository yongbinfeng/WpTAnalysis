"""
script to make postfit comparisons
"""
from modules.postFitScripts import MakePostPlot, result2json, MakeWpTPostFitPlots
from modules.CombineHarvester.plotImpacts import plotImpacts
import ROOT
ROOT.gROOT.SetBatch(True)

doInclusive = True
doWpT = False
doRebin = False

# boolean flag to config if the pull
# distribution should be included in the plots
showPULL = True
if doInclusive:
    suffix = "_Rebin.root" if doRebin else ".root"
    MakePostPlot("cards/datacard_muminus_lepEta_bin0_WpT_bin0" + suffix, "muminus", "", showPULL)
    MakePostPlot("cards/datacard_muplus_lepEta_bin0_WpT_bin0" + suffix, "muplus", "", showPULL)
    MakePostPlot("cards/test_e" + suffix,  "eplus", "", showPULL)
    MakePostPlot("cards/test_5TeV" + suffix, "muplus", "", showPULL)
    
    #result2json("cards/datacard_muminus_lepEta_bin0_WpT_bin0" + suffix, "w_muminus_sig_mu", "cards/impacts_muminus_lepEta_bin0_WpT_bin0.json")
    #result2json("cards/datacard_muplus_lepEta_bin0_WpT_bin0" + suffix, "w_muplus_sig_mu", "cards/impacts_muplus_lepEta_bin0_WpT_bin0.json")
    result2json("cards/datacard_mu_lepEta_bin0_WpT_bin0" + suffix, "w_muminus_sig_mu", "cards/impacts_muminus_lepEta_bin0_WpT_bin0.json")
    result2json("cards/datacard_mu_lepEta_bin0_WpT_bin0" + suffix, "w_muplus_sig_mu", "cards/impacts_muplus_lepEta_bin0_WpT_bin0.json")
    result2json("cards/test_e" + suffix,       "w_eplus_sig_mu",  "cards/impacts_eplus.json")
    result2json("cards/test_e" + suffix,       "w_eminus_sig_mu",  "cards/impacts_eminus.json")

    result2json("cards/test_5TeV" + suffix, "w_muminus_sig_mu", "cards/impacts_muminus_lepEta_bin0_WpT_bin0_5TeV.json")
    result2json("cards/test_5TeV" + suffix, "w_muplus_sig_mu", "cards/impacts_muplus_lepEta_bin0_WpT_bin0_5TeV.json")

    result2json("cards/test_e_5TeV" + suffix, "w_eminus_sig_mu", "cards/impacts_eminus_lepEta_bin0_WpT_bin0_5TeV.json")
    result2json("cards/test_e_5TeV" + suffix, "w_eplus_sig_mu",  "cards/impacts_eplus_lepEta_bin0_WpT_bin0_5TeV.json")

    plotImpacts("cards/impacts_muminus_lepEta_bin0_WpT_bin0.json", "impacts_muminus_lepEta_bin0")
    plotImpacts("cards/impacts_muplus_lepEta_bin0_WpT_bin0.json", "impacts_muplus_lepEta_bin0")
    plotImpacts("cards/impacts_eplus.json",       "impacts_eplus_lepEta")
    plotImpacts("cards/impacts_eminus.json",       "impacts_eminus_lepEta")

    plotImpacts("cards/impacts_muminus_lepEta_bin0_WpT_bin0_5TeV.json", "impacts_muminus_lepEta_bin0_5TeV")
    plotImpacts("cards/impacts_muplus_lepEta_bin0_WpT_bin0_5TeV.json", "impacts_muplus_lepEta_bin0_5TeV")
    
    plotImpacts("cards/impacts_eminus_lepEta_bin0_WpT_bin0_5TeV.json", "impacts_eminus_lepEta_bin0_5TeV")
    plotImpacts("cards/impacts_eplus_lepEta_bin0_WpT_bin0_5TeV.json", "impacts_eplus_lepEta_bin0_5TeV")


    # plot the impacts on QCD
    #result2json("cards/datacard_muplus_lepEta_bin0_WpT_bin0" + suffix, "QCD_muplus_lepEta_bin0_WpT_bin0_mu", "cards/impacts_muplus_lepEta_bin0_WpT_bin0_qcd.json")
    #result2json("cards/datacard_eplus" + suffix,       "QCD_eplus_lepEta_bin1_WpT_bin0_mu",  "cards/impacts_eplus_qcd.json")

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


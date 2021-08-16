"""
script to make postfit comparisons
"""
from PostFits.postFitScripts import MakePostPlot, result2json
from PostFits.CombineHarvester.plotImpacts import plotImpacts
import ROOT
ROOT.gROOT.SetBatch(True)

doWpT = True

# boolean flag to config if the pull
# distribution should be included in the plots
showPULL = False
MakePostPlot("cards/datacard_muplus_lepEta_bin0_WpT_bin0.root", "muplus", "", showPULL)
MakePostPlot("cards/datacard_eplus.root",  "eplus", "", showPULL)
MakePostPlot("cards/datacard_muplus_lepEta_bin0_WpT.root", "muplus", "WpT", showPULL)

result2json("cards/datacard_muplus_lepEta_bin0_WpT_bin0.root", "w_muplus_sig_mu", "cards/impacts_muplus_lepEta_bin0_WpT_bin0.json")
result2json("cards/datacard_eplus.root",       "w_eplus_sig_mu",  "cards/impacts_eplus.json")

plotImpacts("cards/impacts_muplus_lepEta_bin0_WpT_bin0.json", "impacts_muplus_lepEta_bin0")
plotImpacts("cards/impacts_eplus.json",       "impacts_eplus_lepEta")

if doWpT:
    MakePostPlot("cards/datacard_muplus_lepEta_bin0_WpT.root", "muplus", "WpT", showPULL)

    wpttruthbins  = ["WpT_truth_bin1", "WpT_truth_bin2", "WpT_truth_bin3", "WpT_truth_bin4", "WpT_truth_bin5", "WpT_truth_bin6", "WpT_truth_bin7", "WpT_truth_bin8", "WpT_truth_bin9"]
    
    for wpttruthbin in wpttruthbins:
        result2json("cards/datacard_muplus_lepEta_bin0_WpT.root", "w_muplus_{}_sig_mu".format(wpttruthbin), "cards/impacts_muplus_{}_sig_mu.json".format(wpttruthbin))
        plotImpacts("cards/impacts_muplus_{}_sig_mu.json".format(wpttruthbin), "impacts_muplus_{}".format(wpttruthbin))


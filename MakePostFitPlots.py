"""
script to make postfit comparisons
"""
from PostFits.postFitScripts import MakePostPlot, result2json
from PostFits.CombineHarvester.plotImpacts import plotImpacts
import ROOT
ROOT.gROOT.SetBatch(True)

# boolean flag to config if the pull
# distribution should be included in the plots
showPULL = False
MakePostPlot("cards/datacard_muplus_lepEta_bin0_WpT_bin0.root", "muplus", ["lepEta_bin0"], "", showPULL)
MakePostPlot("cards/datacard_eplus.root",  "eplus", ["lepEta_bin1", "lepEta_bin2"], "", showPULL)

result2json("cards/datacard_muplus_lepEta_bin0_WpT_bin0.root", "w_muplus_sig_mu", "cards/impacts_muplus_lepEta_bin0_WpT_bin0.json")
result2json("cards/datacard_eplus.root",       "w_eplus_sig_mu",  "cards/impacts_eplus.json")

plotImpacts("cards/impacts_muplus_lepEta_bin0_WpT_bin0.json", "impacts_muplus_lepEta_bin0")
plotImpacts("cards/impacts_eplus.json",       "impacts_eplus_lepEta")

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
MakePostPlot("Cards/datacard_muplus_lepEta_bin0.root", "muplus", ["lepEta_bin0"], "", showPULL)
MakePostPlot("Cards/datacard_eplus_lepEta.root",  "eplus", ["lepEta_bin1", "lepEta_bin2"], "", showPULL)

result2json("Cards/datacard_muplus_lepEta_bin0.root", "w_muplus_sig_mu", "Cards/impacts_muplus_lepEta_bin0.json")
result2json("Cards/datacard_eplus_lepEta.root",       "w_eplus_sig_mu",  "Cards/impacts_eplus_lepEta.json")

plotImpacts("Cards/impacts_muplus_lepEta_bin0.json", "impacts_muplus_lepEta_bin0")
plotImpacts("Cards/impacts_eplus_lepEta.json",       "impacts_eplus_lepEta")

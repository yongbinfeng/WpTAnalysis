#
# extract the QCD normalization from the postfit distributions,
# in order to normalize QCD properly for the prefit impacts and plots
# A bit hacky, but it works
#

import ROOT
import json

fileName = "forCombine/test/test4/commands/scripts/card_combined.root"
ofileName = "data/QCDNorms.json"

results = {}
QCDNorms = [
    "QCD_muminus_lepEta_bin0_WpT_bin0_13TeV",
    "QCD_muplus_lepEta_bin0_WpT_bin0_13TeV",
    "QCD_eminus_lepEta_bin0_WpT_bin0_13TeV",
    "QCD_eplus_lepEta_bin0_WpT_bin0_13TeV",
    "QCD_muminus_lepEta_bin0_WpT_bin0_5TeV",
    "QCD_muplus_lepEta_bin0_WpT_bin0_5TeV",
    "QCD_eminus_lepEta_bin0_WpT_bin0_5TeV",
    "QCD_eplus_lepEta_bin0_WpT_bin0_5TeV",
]

f = ROOT.TFile(fileName)
t = f.Get("fitresults")
t.GetEntry(0)

for QCDNorm in QCDNorms:
    print(QCDNorm)
    results[QCDNorm] = t.GetLeaf(f"{QCDNorm}_mu").GetValue()

# print results
for key in results:
    print(key, results[key])

# save to json
with open(ofileName, 'w') as f:
    json.dump(results, f, indent = 2, separators=(',', ': '))


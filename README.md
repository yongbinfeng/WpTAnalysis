# W cross section and W pT Analysis

Some analyzers to make plots and run statistical analysis for the inclusive W and Z cross section, and W pT analysis. Still under development.

```
git clone git@github.com:yongbinfeng/WpTAnalysis.git
git checkout master
```

The analyzer is based on the ntuples from this repository: https://github.com/MiT-HEP/MitEwk13TeV/tree/CMSSW_94X


### Make Histograms
The scripts are based on RDataFrame. Tested using `CMSSW_10_6_0` environment and it runs fine.

This histogram makers reply on some functions in `Functions.cc`. Compile it first before using:
```
root -l
root [0] .L Functions.cc+
```
This only needs to be done once.

To produce the histograms in the signal region, just run
```
python MakePlots_Wlnu.py
```
which will produce all the needed histograms for data and MC templates in the muon channel. 

This might take a while since it needs to loop over a large number of events with different systematic variations. Change `doTest` flag to true will run in the test mode, with only a subset of samples and variations. Suitable for debugging.

Change `doMuon` flag to false will run in the electron channel. 

Change `doWpT` flag to true will add some W pT bins at reco level, and also some W pT truth bins for the signal MC. *Note* making the W pT bins in the signal region would take a lot of time. Under investigation.

To produce the histograms (QCD background templates) in the lepton anti-isolated region (CR), run
```
python MakePlots_Wlnu_AntiIso.py
```

### QCD Extrapolation

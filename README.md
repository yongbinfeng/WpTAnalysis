# W cross section and W pT Analysis

Some analyzers to make plots and run statistical analysis for the inclusive W and Z cross section, and W pT analysis. Under development.

```
git clone git@github.com:yongbinfeng/WpTAnalysis.git
git checkout master
```

The analyzer is based on the ntuples from this repository: https://github.com/MiT-HEP/MitEwk13TeV/tree/CMSSW_94X. Some plotting scripts are based on RDataFrame. The combine command needs to be run with tfcombine (https://github.com/bendavid/HiggsAnalysis-CombinedLimit/tree/tensorflowfit). (I'm on branch `tensorflowfit`.). Tested with the `CMSSW_10_6_0` environment and everything works well.


## Make Histograms
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
which will produce all the needed histograms for data and MC templates in the muon channel, saved in the `root` directory. 

This might take a while since it needs to loop over a large number of events with different systematic variations. Change `doTest` flag to true will run in the test mode, with only a subset of samples and variations. Suitable for debugging.

Change `doMuon` flag to false will run in the electron channel. 

Change `doWpT` flag to true will add some W pT bins at reco level, and also some W pT truth bins for the signal MC. *Note* making the W pT bins in the signal region would take a lot of time. Under investigation.

To produce the histograms (QCD background templates) in the lepton anti-isolated region (CR), run
```
python MakePlots_Wlnu_AntiIso.py
```
Some of the preprocessed root files have already been saved in the root directory


## QCD Extrapolation
The `ExtrapolateQCD.py` script is used to extrapolate the QCD template from a set of antiIsolated control regions to the signal region. Run 
```
python ExtrapolateQCD.py
```
would run the extrapolation for the muon plus channel. Change `doMuon` to false would run the eplus channel. Note currently in the electron channel we do the barrel and endcap extrapolation separately.

## Make DataCards
`MakeCards.py` script is used to generate the datacards needed for combine. It takes the data and MC templates from the 1st step, and the QCD extrapolated templates from the 2nd step, and add some systematic uncertainties. 

Then one can run the tfcombine in the `cards` directory, e.g.:
```
text2hdf5.py datacard_muplus_lepEta_bin0_WpT_bin0.txt
combinetf.py datacard_muplus_lepEta_bin0_WpT_bin0.hdf5 --saveHists --doImpacts --output datacard_muplus_lepEta_bin0_WpT_bin0.root
```

Similarly, in the electron channel,
```
text2hdf5.py datacard_eplus.txt 
combinetf.py datacard_eplus.hdf5 --saveHists --doImpacts --output datacard_eplus.root
```

In the W to muplus nu channel for W pT,
```
text2hdf5.py datacard_muplus_lepEta_bin0_WpT.txt
combinetf.py datacard_muplus_lepEta_bin0_WpT.hdf5 --saveHists --doImpacts --output datacard_muplus_lepEta_bin0_WpT.root
```

These would produce root files with `datacard_muplus_lepEta_bin0_WpT_bin0.root`, `datacard_eplus.root`, and `datacard_muplus_lepEta_bin0_WpT.root`, with the pre/postfit distributions and nuisance impacts information saved.

## Make Postfit plots
One can compare the data-Template postfit distributions and the impacts of nuisance parameters, etc by running
```
python MakePostFitPlots.py
```

## ToDO
- Understand better the QCD background and its impact
- Add electron scale uncertainties, and possibly the muon rochester correction uncertainties

# W cross section and W pT Analysis

Collection of the analysis code to make plots and run statistical analysis for the inclusive W and Z cross section, and W pT analysis, ideally end to end. Under development.

```
git clone git@github.com:yongbinfeng/WpTAnalysis.git
git checkout master
```

The code is based on python3, and ROOT is required. On nodes mounted with `/cvmfs` one can set up the environment by
```
source setenv.sh
```

The analysis code takes the ntuples from this repository: https://github.com/MiT-HEP/MitEwk13TeV/tree/CMSSW_94X. Some plotting scripts are based on RDataFrame. The combine command needs to be run with tfcombine (https://github.com/bendavid/HiggsAnalysis-CombinedLimit/tree/tensorflowfit). (I'm on branch `tensorflowfit`.)


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

Change `doElectron` flag will run in the electron channel, e.g.,
```
python MakePlots_Wlnu.py --doElectron
```

Change `doWpT` flag to true will add some W pT bins at reco level, and also some W pT truth bins for the signal MC. *Note* making the W pT bins in the signal region would take a lot of time. Under investigation.

To produce the histograms (QCD background templates) in the lepton anti-isolated region (CR), run
```
python MakePlots_Wlnu_AntiIso.py
```
Do this for electron, and muon channels separately. `--applyScaling` flag will scale the MC xsec by 30\% and the shape difference will be one source of systematic uncertainty of the QCD shapes.

## PrepareFits.py

The histogram processings and data card generations, etc, are all included in `PrepareFits.py`. 
```
python PrepareFits.py
```
It includes processing histograms (rebinning, deal with overflow/underflow, etc), merging tau contributions to the signal process, extrapolating QCD templates from anti-isolated region to signal region, and generate the datacards and commands to run combine. Each of the steps are wrapped into one function, and one can also run each function separately.

## Run combine fits
The commands to run combine, the datacards, and also the root files with histograms and xsecs are all generated in the `PrepareFits.py` step. The commands are saved in `run_combine.sh`. The combine fits has to run under the tfcombine environment (different from the analyzer environment as it is python2.)
```
bash card_combined.sh # for combination of electron and muon channel
bash card_e.sh # for electron channel only
bash card_mu.sh # for muon channel only
```

This will run all fits and produce post-fit histograms, impacts, and grouped systematic uncertainties.

## Make Postfit plots
One can compare the data-Template postfit distributions and the impacts of nuisance parameters, etc by running
```
python MakePostFitPlots.py
```

## ToDO
- Add electron scale uncertainties, and possibly the muon rochester correction uncertainties

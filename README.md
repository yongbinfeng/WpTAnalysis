# W pT Analysis

Some code to make plots for the W pT analysis. Still under development.

The scripts needs RDataFrame. Tested under `CMSSW_10_6_0` environment and it works fine.

```
git clone git@github.com:yongbinfeng/WpTAnalysis.git
git checkout master
```

And compile the `Functions.cc` for the first time:
```
root -l
root [0] .L Functions.cc+
```

The plot making scripts for Zee, Zmm, and Wmnu should work fine. But the Wenu channel script has not been updated yet. To make the plots, simply run
```
python MakePlots_Zmumu.py
```

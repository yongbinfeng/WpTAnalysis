
python FRAndClusure.py --doPtVsEta --useTwoPhis
python FRAndClusure.py --doPtVsEta --useTwoPhis --doElectron

python PrepareFits.py --do13TeV --do5TeV --reweightZpt

python MakePostFitPlots.py --do13TeV --do5TeV --combineSqrtS --doElectron --doMuon --reweightZpt


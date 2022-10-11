from modules.cardMaker import MakeWJetsCards, combineCards
import os
    
doMuon = True
doElectron = False
doWpT = False

if doWpT:
    wptbins = ["WpT_bin1", "WpT_bin2", "WpT_bin3", "WpT_bin4", "WpT_bin5", "WpT_bin6", "WpT_bin7", "WpT_bin8", "WpT_bin9"]
else:
    wptbins = ["WpT_bin0"]


channel = "muplus"
etabin = "lepEta_bin0"
fnames = {
    "muplus" : "output_shapes_mu_mergeTau",
    "muminus": "output_shapes_mu_mergeTau",
    "eplus" :  "output_shapes_e",
    "eminus":  "output_shapes_e",
}

fqcds = {
    "muplus": "QCD/qcdshape_extrapolated_mu",
    "muminus": "QCD/qcdshape_extrapolated_mu",
    "eplus":  "QCD/qcdshape_extrapolated_e",
    "eplus":  "QCD/qcdshape_extrapolated_e",
}

pwd = os.getcwd()
if doMuon:
    channels = ["muplus", "muminus"]
    for channel in channels:
        etabin = "lepEta_bin0"
        pwd = os.getcwd()
        cards = []
        for wptbin in wptbins:
            postfix = "" if not doWpT else "_WpT"
            card = MakeWJetsCards(pwd+"/root/"+fnames[channel]+postfix+".root", pwd+"/root/"+fqcds[channel]+"_"+wptbin+".root", channel, wptbin, etabin, doWpT=doWpT)
            cards.append(card)

        if doWpT:
            # combine the cards in different w pt bins
            combineCards(wptbins, cards, pwd+"/cards/datacard_{}_{}_WpT.txt".format(channel, etabin))
if doElectron:
    # do electron
    channel = "eplus"
    pwd = os.getcwd()
    for wptbin in ["WpT_bin0"]:
        card1 = MakeWJetsCards(pwd+"/root/"+fnames[channel]+".root", pwd+"/root/"+fqcds[channel]+"_"+wptbin+".root", channel, wptbin, "lepEta_bin1")
        card2 = MakeWJetsCards(pwd+"/root/"+fnames[channel]+".root", pwd+"/root/"+fqcds[channel]+"_"+wptbin+".root", channel, wptbin, "lepEta_bin2")
    
    # combine the two cards in 2 eta bins
    combineCards(["Barrel", "Endcap"], [card1, card2], pwd+"/cards/datacard_{}.txt".format(channel))

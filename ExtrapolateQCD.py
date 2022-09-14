from modules.qcdExtrapolater import ExtrapolateQCD

doMuon = True
doElectron = False
doWpT = False
do5TeV = False

if doWpT:
    wptbins = ["WpT_bin1", "WpT_bin2", "WpT_bin3", "WpT_bin4", "WpT_bin5", "WpT_bin6", "WpT_bin7", "WpT_bin8", "WpT_bin9"]
else:
    wptbins = ["WpT_bin0"]

if doMuon:
    fname = "root/output_qcdshape_munu"
    oname = "root/qcdshape_extrapolated_mu"
    for wptbin in wptbins:
        if do5TeV:
            ExtrapolateQCD(fname+"_5TeV.root", oname+"_"+wptbin+"_5TeV.root", ["muplus", "muminus"], wptbin, ["lepEta_bin0"], fname_scaled = fname+"_applyScaling_5TeV.root", is5TeV=True)
        else:
            ExtrapolateQCD(fname+".root", oname+"_"+wptbin+".root", ["muplus", "muminus"], wptbin, ["lepEta_bin0"], fname_scaled = fname+"_applyScaling.root")

if doElectron:
    # on electron channel
    # currently separate barrel and endcap for electrons
    fname = "root/output_qcdshape_enu"
    oname = "root/qcdshape_extrapolated_e"
    for wptbin in ["WpT_bin0"]:
        # for electrons, only implement the inclusive version for now
        if do5TeV:
            ExtrapolateQCD(fname+"_5TeV.root", oname+"_"+wptbin+"_5TeV.root", ["eplus", "eminus"], wptbin, ["lepEta_bin0"], fname_scaled = fname+"_applyScaling_5TeV.root", is5TeV = True)
        else:
            ExtrapolateQCD(fname+".root", oname+"_"+wptbin+".root", ["eplus", "eminus"], wptbin, ["lepEta_bin0"], fname_scaled = fname+"_applyScaling.root")


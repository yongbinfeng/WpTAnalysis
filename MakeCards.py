"""
script to generate the W(lnv) cards for tfCombine
"""
from collections import OrderedDict
import os

doMuon = True
doElectron = False
doWpT = False
doQCDSys = True

if doWpT:
    wptbins = ["WpT_bin1", "WpT_bin2", "WpT_bin3", "WpT_bin4", "WpT_bin5", "WpT_bin6", "WpT_bin7", "WpT_bin8", "WpT_bin9"]
else:
    wptbins = ["WpT_bin0"]


class Process(object):
    """
    define one process and related information for the cards
    """
    def __init__(self,**kwargs):
        self.name = kwargs.get('name', 'sig')
        self.fname = kwargs.get('fname', 'test.root')
        self.hname = kwargs.get('hname', 'htest')
        self.hsys  = kwargs.get('hsys', 'htest_sys')
        self.isObs = kwargs.get('isObs', False) # if it is the observation
        self.isSignal = kwargs.get('isSignal', False) # no contraint on the strength if true
        self.rate = kwargs.get('rate', -1.0)
        self.isMC = kwargs.get('isMC', False) # contraint on lumi if true
        self.isV = kwargs.get('isV', False) # recoil sys if true
        self.isQCD = kwargs.get('isQCD', False) # fake QCD
        self.xsecUnc = kwargs.get('xsecUnc', "0.10") # uncertainty on the cross section

        #  regulation
        if self.isObs:
            self.isSignal = False
            self.isMC = False
            self.isV = False
            self.isQCD = False
            self.xsecUnc = "-1.0"

        if self.isSignal:
            self.isMC = True # signal process is always uing MC template
            self.isV = True
            self.isQCD = False
            self.xsecUnc = "-1.0"

        if self.isQCD:
            self.isMC = False # QCD is always using data driven
            self.isV = False
            self.xsecUnc = "-1.0"


def MakeCards(fname_mc, fname_qcd, channel, wptbin, etabin):
    # data
    data = Process(name = "data_obs", fname = fname_mc,
                   hname = "Data/histo_wjets_{}_mT_1_{}_{}_data".format(channel, wptbin, etabin),
                   isObs = True,
                   )
    # sig processes
    sigs = []
    if not doWpT:
        sig = Process(name = "w_"+channel+"_sig", fname = fname_mc, 
                     hname = "MCTemplates/histo_wjets_{}_mT_1_{}_{}_grouped_wlnu0".format( channel, wptbin, etabin),
                     hsys  = "MCTemplates/histo_wjets_{}_mT_1_{}_{}_grouped_wlnu0_".format(channel, wptbin, etabin),
                     isSignal = True,
                     isMC = True,
                     isV = True,
                     isQCD = False
                     ) 
        sigs.append(sig)
    else:
        wpttruthbins  = ["WpT_truth_bin1", "WpT_truth_bin2", "WpT_truth_bin3", "WpT_truth_bin4", "WpT_truth_bin5", "WpT_truth_bin6", "WpT_truth_bin7", "WpT_truth_bin8", "WpT_truth_bin9"]
        for wpttruthbin in wpttruthbins:
            sig = Process(name = "w_"+channel+"_"+wpttruthbin+"_sig", fname = fname_mc,
                     hname = "{}/histo_wjets_{}_mT_1_{}_{}_{}_signalMC".format( wpttruthbin, channel, wptbin, etabin, wpttruthbin),
                     hsys  = "{}/histo_wjets_{}_mT_1_{}_{}_{}_signalMC_".format(wpttruthbin, channel, wptbin, etabin, wpttruthbin),
                     isSignal = True,
                     isMC = True,
                     isV = True,
                     isQCD = False
                     )
            sigs.append(sig)

    # ttbar bkg
    ttbar = Process(name = "tt", fname = fname_mc,
                    hname = "MCTemplates/histo_wjets_{}_mT_1_{}_{}_grouped_ttbar3".format(channel, wptbin, etabin),
                    hsys = "MCTemplates/histo_wjets_{}_mT_1_{}_{}_grouped_ttbar3_".format(channel, wptbin, etabin),
                    isSignal = False,
                    isMC = True,
                    isV = False,
                    isQCD = False,
                    xsecUnc = "1.10"
                )
    # Z bkg
    zxx = Process(name = "zxx", fname = fname_mc,
                  hname = "MCTemplates/histo_wjets_{}_mT_1_{}_{}_grouped_zxx9".format(channel, wptbin, etabin),
                  hsys = "MCTemplates/histo_wjets_{}_mT_1_{}_{}_grouped_zxx9_".format(channel, wptbin, etabin),
                  isSignal = False,
                  isMC = True,
                  isV = True,
                  isQCD = False,
                  xsecUnc = "1.03"
                  )
    # W->tau + nu bkg
    wtau = Process(name = "taunu", fname = fname_mc,
                   hname = "MCTemplates/histo_wjets_{}_mT_1_{}_{}_grouped_wx10".format(channel, wptbin, etabin),
                   hsys = "MCTemplates/histo_wjets_{}_mT_1_{}_{}_grouped_wx10_".format(channel, wptbin, etabin),
                   isSignal = False,
                   isMC = True,
                   isV = True,
                   isQCD = False,
                   xsecUnc = "1.10"
                )
    # VV bkg
    vv = Process(name = "VV", fname = fname_mc,
                 hname = "MCTemplates/histo_wjets_{}_mT_1_{}_{}_grouped_vv6".format(channel, wptbin, etabin),
                 hsys = "MCTemplates/histo_wjets_{}_mT_1_{}_{}_grouped_vv6_".format(channel, wptbin, etabin),
                 isSignal = False,
                 isMC = True,
                 isV = False,
                 isQCD = False,
                 xsecUnc = "1.10"
                 )
    # QCD bkg
    qcd = Process(name = "QCD_"+channel+"_"+etabin+"_"+wptbin, fname = fname_qcd,
                  hname = "h_QCD_Extrapolated_{}_{}_{}".format(channel, etabin, wptbin),
                  #hsys = "h_QCD_Extrapolated_{}_lepEta_".format(channel),
                  hsys = "h_QCD_Extrapolated_",
                  isSignal = False,
                  isMC = False,
                  isV = False,
                  isQCD = True,
                )

    processes = sigs + [ttbar, zxx, wtau, vv, qcd]
    
    lepname = "mu" if "mu" in channel else "e"

    # correction systematics
    # in Aram's ntuples, defined here: https://github.com/MiT-HEP/MitEwk13TeV/blob/CMSSW_94X/NtupleMod/eleNtupleMod.C#L71
    # and here: https://github.com/MiT-HEP/MitEwk13TeV/blob/CMSSW_94X/NtupleMod/muonNtupleMod.C#L81
    # 1 is MC, 2 is FSR, 3 is Bkg, 4 is tagpt, 
    # 6 is prefire
    sfsys = ["SysWeight1", lepname+"_SysWeight2", lepname+"_SysWeight3", lepname+"_SysWeight4", "SysWeight8", "SysWeight10"]

    # recoil correction systematics
    # in Aram's ntuples, defined here: https://github.com/MiT-HEP/MitEwk13TeV/blob/CMSSW_94X/NtupleMod/eleNtupleMod.C#L65
    # 1 is central correction, 2 is correction with different eta bins, 3 is correction using Gaussian kernels, 
    # 6-15 are corrections with different statistical uncertainties
    recoilsys = ['SysRecoil2', 'SysRecoil3', 'SysRecoil6', 'SysRecoil7', 'SysRecoil8', 'SysRecoil9', 'SysRecoil10', 'SysRecoil11',
                 'SysRecoil12', 'SysRecoil13', 'SysRecoil14', 'SysRecoil15']

    # qcd stat
    # this is hard coded for now. Will be improved later
    #prefix = etabin.split("_")[1]+"_"+wptbin
    prefix = channel + "_" + etabin + "_" + wptbin
    nbins = 12
    qcdstats = [prefix+"_bin"+str(i)+"shape" for i in range(1, nbins+1)]
    if doQCDSys:
        #qcdstats += [prefix+"_Pol2shape"]
        qcdstats += [prefix+"_ScaledMCshape"]


    #
    # writing datacards
    #
    if not os.path.exists("cards"):
        os.makedirs("cards")
    ofile = open("cards/datacard_{}_{}_{}.txt".format(channel, etabin, wptbin), "w")
    ofile.write("imax 1 number of channels\n")

    ofile.write("jmax {} number of processes -1\n".format(len(processes)-1))
    ofile.write("kmax * number of nuisance parameters\n\n")
    # observation
    ofile.write("Observation -1\n\n")
    ofile.write("shapes {pname:<30} * {fname:<30} {hname:<50}\n".format(pname = "data_obs", fname = data.fname, hname = data.hname))

    for proc in processes:
        ofile.write("shapes {pname:<30} * {fname:<30} {hname:<50} {hsys}$SYSTEMATIC\n".format(pname = proc.name, fname = proc.fname, hname = proc.hname, hsys = proc.hsys))
    ofile.write("\n")

    # systematic uncertainties
    lines = OrderedDict()
    lines["bin"]        = "{:<40}".format("bin")
    lines["proc_name"]  = "{:<40}".format("process")
    lines["proc_index"] = "{:<40}".format("process")
    lines["rate"]       = "{:<40}".format("rate")

    lines['lumi'] = "{:<30} {:<10}".format("lumi_13TeV", "lnN")
    for proc in processes:
        if not proc.isSignal and proc.isMC:
            lines["norm_"+proc.name] = "{:<30} {:<10}".format("norm_"+proc.name, "lnN")
    
    for sys in sfsys:
        lines[sys] = "{:<30} {:<10}".format(sys, "shape")
    lines["effstat"] = "{:<30} {:<10}".format("effstat", "lnN")

    for sys in recoilsys:
        lines[sys] = "{:<30} {:<10}".format(sys, "shape")

    for sys in qcdstats:
        lines[sys] = "{:<30} {:<10}".format(sys, "shape")

    # fill in the per-process information
    isig = 0
    ibkg = 1
    for proc in processes:
        lines["bin"] += " {binname:<10}".format(binname="bin1")
        lines["proc_name"] += " {pname:<10}".format(pname = proc.name)
        if proc.isSignal or proc.isQCD:
            lines["proc_index"] += " {pindex:<10}".format(pindex = isig)
            isig -= 1
        else:
            lines["proc_index"] += " {pindex:<10}".format(pindex = ibkg)
            ibkg += 1
        lines["rate"] += " {:<10}".format("-1.0")

        temp = "1.017 "if proc.isMC else "-"
        lines['lumi'] += " {:<10}".format(temp)
        for proc2 in processes:
            if not proc2.isSignal and proc2.isMC:
                temp = proc.xsecUnc if proc.name == proc2.name else "-"
                lines["norm_"+proc2.name] += " {:<10}".format(temp)

        for sys in sfsys:
            temp = "1.0" if not proc.isQCD else "-"
            lines[sys] += " {:<10}".format(temp)

        temp = "-"
        if proc.isSignal:
            temp = "1.0049" if channel.startswith("e") else "1.0029" 
        lines['effstat'] += " {:<10}".format(temp)

        for sys in recoilsys:
            temp = "1.0" if proc.isV else "-"
            lines[sys] += " {:<10}".format(temp)

        for sys in qcdstats:
            temp = "1.0" if proc.isQCD else "-"
            lines[sys] += " {:<10}".format(temp)

    # add the QCD group 
    qcdgroup = "QCD group ="
    for sys in qcdstats:
        qcdgroup += " " + sys

    # add the recoil group
    recoilgroup = "Recoil group ="
    for sys in recoilsys:
        recoilgroup += " " + sys

    # write these systematic lines
    for line in list(lines.values()):
        ofile.write(line+"\n")
    ofile.write(qcdgroup+"\n")
    ofile.write(recoilgroup+"\n")

    ofile.close()
    return ofile.name.split('/')[-1]


def combineCards(labels, cards, oname):
    """
    genrate and run the command to combine cards in different bins
    """
    cmd = "cd cards; combineCards.py"
    for label, card in zip(labels, cards):
        cmd += " {}={}".format(label, card)
    cmd += " > {}; cd ../".format(oname)
    print(("combine cards with {}".format(cmd)))
    os.system(cmd)
    
    
if __name__ == "__main__":
    channel = "muplus"
    etabin = "lepEta_bin0"
    fnames = {
        "muplus" : "output_shapes_mu",
        "muminus": "output_shapes_mu",
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
                card = MakeCards(pwd+"/root/"+fnames[channel]+postfix+".root", pwd+"/root/"+fqcds[channel]+"_"+wptbin+".root", channel, wptbin, etabin)
                cards.append(card)

            if doWpT:
                # combine the cards in different w pt bins
                combineCards(wptbins, cards, pwd+"/cards/datacard_{}_{}_WpT.txt".format(channel, etabin))
    if doElectron:
        # do electron
        channel = "eplus"
        pwd = os.getcwd()
        for wptbin in ["WpT_bin0"]:
            card1 = MakeCards(pwd+"/root/"+fnames[channel]+".root", pwd+"/root/"+fqcds[channel]+"_"+wptbin+".root", channel, wptbin, "lepEta_bin1")
            card2 = MakeCards(pwd+"/root/"+fnames[channel]+".root", pwd+"/root/"+fqcds[channel]+"_"+wptbin+".root", channel, wptbin, "lepEta_bin2")
        
        # combine the two cards in 2 eta bins
        combineCards(["Barrel", "Endcap"], [card1, card2], pwd+"/cards/datacard_{}.txt".format(channel))

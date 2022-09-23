"""
code to generate the W(lnv) cards for tfCombine. 
The desig can be improved, for example: different processes and 
different systematics uncertainties, with maps between i process and j systematic.
Then the processes and systematics can be naturally grouped into different subcategories
"""

from collections import OrderedDict
import os

class Process(object):
    """
    define one process and related information for the cards
    """
    def __init__(self,**kwargs):
        self.name = kwargs.get('name', 'sig')
        self.hname = kwargs.get('hname', 'htest')
        self.hsys  = kwargs.get('hsys', 'htest_sys')
        self.isObs = kwargs.get('isObs', False) # if it is the observation
        self.isSignal = kwargs.get('isSignal', False) # no contraint on the strength if true
        self.rate = kwargs.get('rate', -1.0)
        self.isMC = kwargs.get('isMC', False) # contraint on lumi if true
        self.isV = kwargs.get('isV', False) # recoil sys if true
        self.isQCD = kwargs.get('isQCD', False) # fake QCD
        self.xsecUnc = kwargs.get('xsecUnc', "0.10") # uncertainty on the cross section

        # a bit hacky, to point to the correct file location
        # - by default the fits run on lpc
        prefix_LPC = "/uscms/home/yfeng/nobackup/WpT/Cards/TestCode/WpTAnalysis/"
        self.fname = prefix_LPC + kwargs.get('fname', 'test.root')

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

class Nuisance(object):
    """
    define one nuisance and related information for making cards
    """
    def __init__(self, **kwargs):
        self.name = kwargs.get('name', 'lumi')
        self.type = kwargs.get('type', 'lnN')
        assert self.type in ['lnN', 'shape'], "Nuisance type can only be lnN or shape"
        # value is only useful for lnN variation
        self.value = kwargs.get('value', '1.026')

def MakeCards(fname_mc, fname_qcd, channel, wptbin, etabin, doWpT = False, rebinned = False, is5TeV = False, nMTBins = 9, outdir = "cards"):
    prefix = ""
    if rebinned:
        prefix = "Rebinned_"

    # from lumi    
    unc_lumi = {}
    unc_lumi['lumi_13TeV'] = 1.017
    unc_lumi['lumi_5TeV'] = 1.019

    # from eff calculations
    unc_effstat = {}
    unc_effstat['effstat_muplus_13TeV'] = 1.0022
    unc_effstat['effstat_muminus_13TeV'] = 1.0021
    unc_effstat['effstat_eplus_13TeV'] = 1.0059
    unc_effstat['effstat_eminus_13TeV'] = 1.0053
    unc_effstat['effstat_muplus_13TeV'] = 1.0023
    unc_effstat['effstat_muminus_13TeV'] = 1.0021
    unc_effstat['effstat_eplus_13TeV'] = 1.0080
    unc_effstat['effstat_eminus_13TeV'] = 1.0077

    # data
    data = Process(name = "data_obs", fname = fname_mc,
                   hname = prefix + "histo_wjets_{}_mT_1_{}_{}_data".format(channel, wptbin, etabin),
                   isObs = True,
                   )
    # sig processes
    sigs = []
    if not doWpT:
        sig = Process(name = "w_"+channel+"_sig", fname = fname_mc, 
                     hname = prefix + "histo_wjets_{}_mT_1_{}_{}_grouped_wsig".format( channel, wptbin, etabin),
                     hsys  = prefix + "histo_wjets_{}_mT_1_{}_{}_grouped_wsig_".format(channel, wptbin, etabin),
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
                     hname = prefix + "histo_wjets_{}_mT_1_{}_{}_{}_signalMC".format( wpttruthbin, channel, wptbin, etabin, wpttruthbin),
                     hsys  = prefix + "histo_wjets_{}_mT_1_{}_{}_{}_signalMC_".format(wpttruthbin, channel, wptbin, etabin, wpttruthbin),
                     isSignal = True,
                     isMC = True,
                     isV = True,
                     isQCD = False
                     )
            sigs.append(sig)

    # ttbar bkg
    sname = "ttbar3" if not is5TeV else "ttbar1"
    ttbar = Process(name = "tt", fname = fname_mc,
                    hname = prefix + "histo_wjets_{}_mT_1_{}_{}_grouped_{}".format(channel, wptbin, etabin, sname),
                    hsys  = prefix + "histo_wjets_{}_mT_1_{}_{}_grouped_{}_".format(channel, wptbin, etabin, sname),
                    isSignal = False,
                    isMC = True,
                    isV = False,
                    isQCD = False,
                    xsecUnc = "1.10"
                )
    # Z bkg
    sname = "zxx9" if not is5TeV else "zxx6"
    zxx = Process(name = "zxx", fname = fname_mc,
                  hname = prefix + "histo_wjets_{}_mT_1_{}_{}_grouped_{}".format(channel, wptbin, etabin, sname),
                  hsys  = prefix + "histo_wjets_{}_mT_1_{}_{}_grouped_{}_".format(channel, wptbin, etabin, sname),
                  isSignal = False,
                  isMC = True,
                  isV = True,
                  isQCD = False,
                  xsecUnc = "1.03"
                  )
    # W->tau + nu bkg
    sname = "wx10" if not is5TeV else "wx7"
    wtau = Process(name = "taunu", fname = fname_mc,
                   hname = prefix + "histo_wjets_{}_mT_1_{}_{}_grouped_{}".format(channel, wptbin, etabin, sname),
                   hsys  = prefix + "histo_wjets_{}_mT_1_{}_{}_grouped_{}_".format(channel, wptbin, etabin, sname),
                   isSignal = False,
                   isMC = True,
                   isV = True,
                   isQCD = False,
                   xsecUnc = "1.10"
                )
    # VV bkg
    sname = "vv6" if not is5TeV else "vv2"
    vv = Process(name = "VV", fname = fname_mc,
                 hname = prefix + "histo_wjets_{}_mT_1_{}_{}_grouped_{}".format(channel, wptbin, etabin, sname),
                 hsys  = prefix + "histo_wjets_{}_mT_1_{}_{}_grouped_{}_".format(channel, wptbin, etabin, sname),
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

    #processes = sigs + [ttbar, zxx, wtau, vv, qcd]
    processes = sigs + [ttbar, zxx, vv, qcd]
    
    lepname = "mu" if "mu" in channel else "e"
    era = "13TeV" if not is5TeV else "5TeV"

    sysgroups = OrderedDict()

    sysgroups["lumisys"] = ["lumi_" + era]

    sysgroups["mcsecsys"] = []
    for proc in processes:
        if not proc.isSignal and proc.isMC:
            sysgroups["mcsecsys"].append("norm_" + proc.name)

    # correction systematics
    # in Aram's ntuples, defined here: https://github.com/MiT-HEP/MitEwk13TeV/blob/CMSSW_94X/NtupleMod/eleNtupleMod.C#L71
    # and here: https://github.com/MiT-HEP/MitEwk13TeV/blob/CMSSW_94X/NtupleMod/muonNtupleMod.C#L81
    # 1 is MC, 2 is FSR, 3 is Bkg, 4 is tagpt, 
    # 8 is ecal prefire
    # 10 is muon prefire
    sysgroups["sfsys"] = ["SysWeight1", lepname+"_SysWeight2", lepname+"_SysWeight3", lepname+"_SysWeight4", "SysWeight8", "SysWeight10"]

    sysgroups["sfstatsys"] = ["effstat_" + channel + "_" + era]

    # recoil correction systematics
    # in Aram's ntuples, defined here: https://github.com/MiT-HEP/MitEwk13TeV/blob/CMSSW_94X/NtupleMod/eleNtupleMod.C#L65
    # 1 is central correction, 2 is correction with different eta bins, 3 is correction using Gaussian kernels, 
    # 6-15 are corrections with different statistical uncertainties
    sysgroups["recoilsys"] = ['SysRecoil2', 'SysRecoil3', 'SysRecoil6', 'SysRecoil7', 'SysRecoil8', 'SysRecoil9', 'SysRecoil10', 'SysRecoil11', 'SysRecoil12', 'SysRecoil13', 'SysRecoil14', 'SysRecoil15']

    # qcd stat
    # this is hard coded for now. Will be improved later
    #prefix = etabin.split("_")[1]+"_"+wptbin
    prefix = channel + "_" + etabin + "_" + wptbin
    nbins = nMTBins
    sysgroups["qcdstats"] = [prefix+"_bin"+str(i)+"shape" for i in range(1, nbins+1)]
    if True:
        #qcdstats += [prefix+"_Pol2shape"]
        sysgroups["qcdstats"] += [prefix+"_ScaledMCshape"]

    # theory systematics
    qcdscale_indices = [1, 2, 3, 4, 6, 8]
    sysgroups["qcdscalesys"] = ["TheoryUnc" + str(idx) for idx in qcdscale_indices]

    pdf_indices = list(range(9, 109))
    sysgroups["pdfsys"] = ["TheoryUnc" + str(idx) for idx in pdf_indices]

    sysgroups["othersys"] = ['SysTauFrac']

    #
    # writing datacards
    #
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    cardname = f"{outdir}/datacard_{channel}_{etabin}_{wptbin}.txt"
    if is5TeV:
        cardname = f"{outdir}/datacard_{channel}_{etabin}_{wptbin}_5TeV.txt"
    ofile = open(cardname, "w")
    ofile.write("imax 1 number of channels\n")

    ofile.write("jmax {} number of processes -1\n".format(len(processes)-1))
    ofile.write("kmax * number of nuisance parameters\n\n")
    # observation
    ofile.write("Observation -1\n\n")
    ofile.write("shapes {pname:<30} * {fname:<40} {hname:<50}\n".format(pname = "data_obs", fname = data.fname, hname = data.hname))

    for proc in processes:
        ofile.write("shapes {pname:<30} * {fname:<40} {hname:<50} {hsys}$SYSTEMATIC\n".format(pname = proc.name, fname = proc.fname, hname = proc.hname, hsys = proc.hsys))
    ofile.write("\n")

    lines = OrderedDict()
    lines["bin"]        = "{:<40}".format("bin")
    lines["proc_name"]  = "{:<40}".format("process")
    lines["proc_index"] = "{:<40}".format("process")
    lines["rate"]       = "{:<40}".format("rate")

    # systematic uncertainties
    for sys in sysgroups["lumisys"]:
        lines[sys] = "{:<30} {:<10}".format(sys, "lnN")

    for sys in sysgroups["mcsecsys"]:
        lines[sys] = "{:<30} {:<10}".format(sys, "lnN")
    
    for sys in sysgroups["sfsys"]:
        lines[sys] = "{:<30} {:<10}".format(sys, "shape")
    
    for sys in sysgroups["sfstatsys"]:
        lines[sys] = "{:<30} {:<10}".format(sys, "lnN")

    for sys in sysgroups["recoilsys"]:
        lines[sys] = "{:<30} {:<10}".format(sys, "shape")

    for sys in sysgroups["qcdstats"]:
        lines[sys] = "{:<30} {:<10}".format(sys, "shape")
    
    for sys in sysgroups["qcdscalesys"]:
        lines[sys] = "{:<30} {:<10}".format(sys, "shape")
    
    for sys in sysgroups["pdfsys"]:
        lines[sys] = "{:<30} {:<10}".format(sys, "shape")

    for sys in sysgroups["othersys"]:
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

        temp = "-"
        for lumi in sysgroups["lumisys"]:
            if proc.isMC:
                temp = unc_lumi[lumi]
            lines[lumi] += " {:<10}".format(temp)

        for proc2 in processes:
            if not proc2.isSignal and proc2.isMC:
                temp = proc.xsecUnc if proc.name == proc2.name else "-"
                lines["norm_"+proc2.name] += " {:<10}".format(temp)

        for sys in sysgroups["sfsys"]:
            temp = "1.0" if not proc.isQCD else "-"
            lines[sys] += " {:<10}".format(temp)

        for sys in sysgroups["sfstatsys"]:
            temp = "-"
            if proc.isSignal:
                temp = unc_effstat[sys]
            lines[sys] += " {:<10}".format(temp)

        for sys in sysgroups["recoilsys"]:
            temp = "1.0" if proc.isV else "-"
            lines[sys] += " {:<10}".format(temp)

        for sys in sysgroups["qcdstats"]:
            temp = "1.0" if proc.isQCD else "-"
            lines[sys] += " {:<10}".format(temp)

        for sys in sysgroups["qcdscalesys"]:
            temp = "1.0" if proc.isSignal else "-"
            lines[sys] += " {:<10}".format(temp)
        
        for sys in sysgroups["pdfsys"]:
            temp = "1.0" if proc.isSignal else "-"
            lines[sys] += " {:<10}".format(temp)

        for sys in sysgroups["othersys"]:
            temp = "1.0" if proc.isSignal else "-"
            lines[sys] += " {:<10}".format(temp)

    # write these systematic lines
    for line in list(lines.values()):
        ofile.write(line+"\n")

    for sysgroup, systematics in sysgroups.items():
        ofile.write(sysgroup+" group =")
        for sys in systematics:
            ofile.write(" " + sys)
        ofile.write("\n")
    ofile.close()
    return cardname

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

def GenerateRunCommand(card_plus: str, card_minus: str, prefix: str = "./"):
    """
    generate the script with commands to run combine.
    inputs can probably be improved here
    """
    workdir = prefix + card_plus.rpartition('/')[0]
    card_plus = card_plus.split('/')[-1]
    card_minus = card_minus.split('/')[-1]
    output = card_plus.replace("plus", "combined").replace(".txt", "")

    cmd = ""
    cmd += "#!/bin/bash\n\n"
    cmd += f"cd {workdir}\n"
    cmd += f"combineCards.py plus={card_plus} minus={card_minus} > {output}.txt\n"
    cmd += f"text2hdf5.py {output}.txt\n"
    cmd += f"combinetf.py {output}.hdf5 --binByBinStat --computeHistErrors --saveHists --doImpacts --output {output}.root\n"

    # write outputs to scripts
    with open(f"{workdir}/{output}.sh", "w") as f:
        f.write(cmd)
    os.system(f"chmod +x {workdir}/{output}.sh")

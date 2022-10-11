"""
code to generate the W(lnv) cards for tfCombine. 
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
        # valuemap saves how much one process is affected by the nuisance parameter
        # key is the process name, value is the uncertainty
        self.valuemap = kwargs.get('valuemap', {})
    
    def __getitem__(self, key):
        return self.valuemap.get(key, '-')
    
    def __setitem__(self, key, val):
        self.valuemap[key] = val


def WriteCard(data, processes, nuisgroups, cardname):
    """
    write the datacard provided with data, processes, and nuisances
    """
    dirpath = cardname.rpartition('/')[0]
    if not os.path.exists(dirpath):
        print(f"Make the directory {dirpath}")
        os.makedirs(dirpath)

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

    ofile.write("{:<40}".format("bin"))
    for proc in processes:
        # one data card is one bin
        ofile.write(" {binname:<10}".format(binname="bin1"))
    ofile.write("\n")

    ofile.write("{:<40}".format("process"))
    for proc in processes:
        ofile.write(" {pname:<10}".format(pname=proc.name))
    ofile.write("\n")

    ofile.write("{:<40}".format("process"))
    isig = 0
    ibkg = 1
    for proc in processes:
        if proc.isSignal or proc.isQCD:
            ofile.write(f" {isig:<10}")
            isig -= 1
        else:
            ofile.write(f" {ibkg:<10}")
            ibkg += 1
    ofile.write("\n")

    ofile.write("{:<40}".format("rate"))
    for proc in processes:
        ofile.write(" {:<10}".format("-1.0"))
    ofile.write("\n\n")

    ## write out all systematics
    for nuisgroup in nuisgroups.values():
        for nuis in nuisgroup:
            ofile.write(f"{nuis.name:<30} {nuis.type:<10}")
            for proc in processes:
                ofile.write(" {:<10}".format(nuis[proc.name]))
            ofile.write("\n")
    ofile.write("\n")

    # write out nuisance groups
    for gname, gnuisances in nuisgroups.items():
        ofile.write(f"{gname:<30} group =")
        for nuis in gnuisances:
            ofile.write(f" {nuis.name}")
        ofile.write("\n")

    ofile.close()
    

def MakeWJetsCards(fname_mc, fname_qcd, channel, wptbin, etabin, doWpT = False, rebinned = False, is5TeV = False, nMTBins = 9, outdir = "cards"):
    # prefix of all histo names
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
    unc_effstat['effstat_muplus_5TeV'] = 1.0023
    unc_effstat['effstat_muminus_5TeV'] = 1.0021
    unc_effstat['effstat_eplus_5TeV'] = 1.0080
    unc_effstat['effstat_eminus_5TeV'] = 1.0077

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

    # list of all processes
    #processes = sigs + [ttbar, zxx, wtau, vv, qcd]
    processes = sigs + [ttbar, zxx, vv, qcd]
    
    lepname = "mu" if "mu" in channel else "e"
    era = "13TeV" if not is5TeV else "5TeV"

    # define different nuisance parameters, their groups,
    # and the impacts on each process
    nuisgroups = OrderedDict()

    nuis_lumi = Nuisance(name = "lumi_" + era, type = "lnN")
    for proc in processes:
        if proc.isMC:
            nuis_lumi[proc.name] = unc_lumi[nuis_lumi.name]
    nuisgroups["lumisys"] = [nuis_lumi]

    nuisgroups["mcsecsys"] = []
    for proc in processes:
        if not proc.isSignal and proc.isMC:
            nuis_norm = Nuisance(name = "norm_" + proc.name, type = "lnN")
            nuis_norm[proc.name] = proc.xsecUnc
            nuisgroups["mcsecsys"].append(nuis_norm)

    # correction systematics
    # in Aram's ntuples, defined here: https://github.com/MiT-HEP/MitEwk13TeV/blob/CMSSW_94X/NtupleMod/eleNtupleMod.C#L71
    # and here: https://github.com/MiT-HEP/MitEwk13TeV/blob/CMSSW_94X/NtupleMod/muonNtupleMod.C#L81
    # 1 is MC, 2 is FSR, 3 is Bkg, 4 is tagpt, 
    # 8 is ecal prefire
    # 10 is muon prefire
    # correlate MC systematic, ECAL and muon prefire 
    # between electron and muon channel
    # uncorrelate FSR, bkg, tagpt in electron ahd muon channel
    nuis_SysWeight1 = Nuisance(name = "SysWeight1", type = "shape")
    nuis_SysWeight2 = Nuisance(name = lepname + "_SysWeight2", type = "shape")
    nuis_SysWeight3 = Nuisance(name = lepname + "_SysWeight3", type = "shape")
    nuis_SysWeight4 = Nuisance(name = lepname + "_SysWeight4", type = "shape")
    nuis_SysWeight8 = Nuisance(name = "SysWeight8", type = "shape")
    nuis_SysWeight10 = Nuisance(name = "SysWeight10", type = "shape")
    nuisgroups["sfsys"] = [nuis_SysWeight1, nuis_SysWeight2, nuis_SysWeight3, nuis_SysWeight4, nuis_SysWeight8, nuis_SysWeight10]
    for proc in processes:
        if not proc.isQCD:
            # all the samples except the QCD apply the corrections
            # so they should be affected by sf systematics
            for sysweight in nuisgroups["sfsys"]:
                sysweight[proc.name] = 1.0

    nuis_effstat = Nuisance(name = "effstat_" + channel + "_" + era, type = "lnN")
    for proc in processes:
        if proc.isSignal:
            # only apply eff/sf stat uncertainty to signals for now
            nuis_effstat[proc.name] = unc_effstat[nuis_effstat.name]
    nuisgroups["sfstatsys"] = [nuis_effstat]

    # recoil correction systematics
    # in Aram's ntuples, defined here: https://github.com/MiT-HEP/MitEwk13TeV/blob/CMSSW_94X/NtupleMod/eleNtupleMod.C#L65
    # 1 is central correction, 2 is correction with different eta bins, 3 is correction using Gaussian kernels, 
    # 6-15 are corrections with different statistical uncertainties
    nuisgroups["recoilsys"] = []
    recoil_indices = ["2", "3", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15"]
    for iuvar in recoil_indices:
        nuis_SysRecoil = Nuisance(name = "SysRecoil" + iuvar, type = "shape")
        for proc in processes:
            if proc.isV:
                # all the V processes (including W and Z) apply the recoil corrections
                nuis_SysRecoil[proc.name] = 1.0
        nuisgroups["recoilsys"].append(nuis_SysRecoil)

    # qcd stat
    # this is hard coded for now. Will be improved later
    nuisgroups["qcdstats"] = []
    prefix = channel + "_" + etabin + "_" + wptbin
    nbins = nMTBins
    for ibin in range(1, nbins+1):
        nuis_QCDStat = Nuisance(name = prefix + f"_bin{ibin}shape", type = "shape")
        for proc in processes:
            if proc.isQCD:
                nuis_QCDStat[proc.name] = 1.0
        nuisgroups["qcdstats"].append(nuis_QCDStat)

    if True:
        nuis_QCDMCCont = Nuisance(name = prefix + "_ScaledMCshape", type = "shape")
        for proc in processes:
            if proc.isQCD:
                nuis_QCDMCCont[proc.name] = 1.0
        nuisgroups["qcdstats"].append(nuis_QCDMCCont)

    # theory systematics
    # qcd scale
    nuisgroups["qcdscalesys"] = []
    for par in ["MuF", "MuR", "MuFMuR"]:
        nuis_QCDScale = Nuisance(name = par, type = "shape")
        for proc in processes:
            if proc.isSignal:
                nuis_QCDScale[proc.name] = 1.0
        nuisgroups["qcdscalesys"].append(nuis_QCDScale)

    # pdf variations
    nuisgroups["pdfsys"] = []
    pdf_indices = list(range(1, 101))
    for ipdf in pdf_indices:
        nuis_PDF = Nuisance(name = f"PDF{ipdf}", type = "shape")
        for proc in processes:
            if proc.isSignal:
                nuis_PDF[proc.name] = 1.0
        nuisgroups["pdfsys"].append(nuis_PDF)
    
    # alphaS variations
    nuisgroups["alphaS"] = []
    nuis_alphaS = Nuisance(name = "alphaS", type = "shape")
    for proc in processes:
        if proc.isSignal:
            nuis_alphaS[proc.name] = 1.0
    nuisgroups["alphaS"].append(nuis_alphaS)

    # tau fraction variation in the signal process
    nuis_TauFrac = Nuisance(name = "SysTauFrac", type = "shape")
    for proc in processes:
        if proc.isSignal:
            nuis_TauFrac[proc.name] = 1.0
    nuisgroups["othersys"] = [nuis_TauFrac]

    #
    # writing datacards
    #
    cardname = f"{outdir}/datacard_{channel}_{etabin}_{wptbin}.txt"
    if is5TeV:
        cardname = f"{outdir}/datacard_{channel}_{etabin}_{wptbin}_5TeV.txt"
    WriteCard(data, processes, nuisgroups, cardname)

    return cardname


def MakeZJetsCards(fname, channel, rebinned = False, is5TeV = False, outdir = "cards"):
    """
    Generate the combine datacard for Z+jets signal region
    """
    # prefix of all histo names
    prefix = ""
    if rebinned:
        prefix = "Rebinned_"

    # from lumi    
    unc_lumi = {}
    unc_lumi['lumi_13TeV'] = 1.017
    unc_lumi['lumi_5TeV'] = 1.019

    # from eff calculations
    unc_effstat = {}
    unc_effstat['effstat_muplus_13TeV'] = 1.0020
    unc_effstat['effstat_muminus_13TeV'] = 1.0020
    unc_effstat['effstat_eplus_13TeV'] = 1.0041
    unc_effstat['effstat_eminus_13TeV'] = 1.0041
    unc_effstat['effstat_muplus_5TeV'] = 1.0019
    unc_effstat['effstat_muminus_5TeV'] = 1.0019
    unc_effstat['effstat_eplus_5TeV'] = 1.0061
    unc_effstat['effstat_eminus_13TeV'] = 1.0061

    # data
    data = Process(name = "data_obs", fname = fname,
                   hname = prefix + f"histo_zjets_zmass_{channel}_weight_0_data",
                   isObs = True,
                   )
    # sig processes
    sigs = []
    sig = Process(name = "z_"+channel+"_sig", fname = fname, 
                 hname = prefix + f"histo_zjets_zmass_{channel}_weight_0_grouped_DY0",
                 hsys = prefix + f"histo_zjets_zmass_{channel}_weight_0_grouped_DY0_",
                 isSignal = True,
                 isMC = True,
                 isV = True,
                 isQCD = False
                 ) 
    sigs.append(sig)

    # ttbar bkg
    sname = "ttbar1" if not is5TeV else "ttbar1"
    ttbar = Process(name = "tt", fname = fname,
                    hname = prefix + f"histo_zjets_zmass_{channel}_weight_0_grouped_ttbar1",
                    hsys = prefix + f"histo_zjets_zmass_{channel}_weight_0_grouped_ttbar1_",
                    isSignal = False,
                    isMC = True,
                    isV = False,
                    isQCD = False,
                    xsecUnc = "1.10"
                )
    # EWK bkg
    sname = "EWK2" if not is5TeV else "EWK2"
    ewk = Process(name = "EWK", fname = fname,
                 hname = prefix + f"histo_zjets_zmass_{channel}_weight_0_grouped_EWK2",
                 hsys = prefix + f"histo_zjets_zmass_{channel}_weight_0_grouped_EWK2_",
                 isSignal = False,
                 isMC = True,
                 isV = False,
                 isQCD = False,
                 xsecUnc = "1.10"
                 )

    # list of all processes
    #processes = sigs + [ttbar, zxx, wtau, vv, qcd]
    processes = sigs + [ttbar, ewk]
    
    lepname = "mu" if "mu" in channel else "e"
    era = "13TeV" if not is5TeV else "5TeV"

    # define different nuisance parameters, their groups,
    # and the impacts on each process
    nuisgroups = OrderedDict()

    nuis_lumi = Nuisance(name = "lumi_" + era, type = "lnN")
    for proc in processes:
        if proc.isMC:
            nuis_lumi[proc.name] = unc_lumi[nuis_lumi.name]
    nuisgroups["lumisys"] = [nuis_lumi]

    nuisgroups["mcsecsys"] = []
    for proc in processes:
        if not proc.isSignal and proc.isMC:
            nuis_norm = Nuisance(name = "norm_" + proc.name, type = "lnN")
            nuis_norm[proc.name] = proc.xsecUnc
            nuisgroups["mcsecsys"].append(nuis_norm)

    # correction systematics
    # in Aram's ntuples, defined here: https://github.com/MiT-HEP/MitEwk13TeV/blob/CMSSW_94X/NtupleMod/eleNtupleMod.C#L71
    # and here: https://github.com/MiT-HEP/MitEwk13TeV/blob/CMSSW_94X/NtupleMod/muonNtupleMod.C#L81
    # 1 is MC, 2 is FSR, 3 is Bkg, 4 is tagpt, 
    # 8 is ecal prefire
    # 10 is muon prefire
    # correlate MC systematic, ECAL and muon prefire 
    # between electron and muon channel
    # uncorrelate FSR, bkg, tagpt in electron ahd muon channel
    nuis_SysWeight1 = Nuisance(name = "SysWeight1", type = "shape")
    nuis_SysWeight2 = Nuisance(name = lepname + "_SysWeight2", type = "shape")
    nuis_SysWeight3 = Nuisance(name = lepname + "_SysWeight3", type = "shape")
    nuis_SysWeight4 = Nuisance(name = lepname + "_SysWeight4", type = "shape")
    nuis_SysWeight8 = Nuisance(name = "SysWeight8", type = "shape")
    nuis_SysWeight10 = Nuisance(name = "SysWeight10", type = "shape")
    nuisgroups["sfsys"] = [nuis_SysWeight1, nuis_SysWeight2, nuis_SysWeight3, nuis_SysWeight4, nuis_SysWeight8, nuis_SysWeight10]
    for proc in processes:
        if not proc.isQCD:
            # all the samples except the QCD apply the corrections
            # so they should be affected by sf systematics
            for sysweight in nuisgroups["sfsys"]:
                sysweight[proc.name] = 1.0

    nuis_effstat_plus = Nuisance(name = "effstat_" + lepname + "plus_" + era, type = "lnN")
    nuis_effstat_minus = Nuisance(name = "effstat_" + lepname + "minus_" + era, type = "lnN")
    for proc in processes:
        if proc.isSignal:
            # only apply eff/sf stat uncertainty to signals for now
            nuis_effstat_plus[proc.name] = unc_effstat[nuis_effstat_plus.name]
            nuis_effstat_minus[proc.name] = unc_effstat[nuis_effstat_minus.name]
    nuisgroups["sfstatsys"] = [nuis_effstat_plus, nuis_effstat_minus]

    # theory systematics
    # qcd scale
    nuisgroups["qcdscalesys"] = []
    for par in ["MuF", "MuR", "MuFMuR"]:
        nuis_QCDScale = Nuisance(name = par, type = "shape")
        for proc in processes:
            if proc.isSignal:
                nuis_QCDScale[proc.name] = 1.0
        nuisgroups["qcdscalesys"].append(nuis_QCDScale)

    # pdf variations
    nuisgroups["pdfsys"] = []
    pdf_indices = list(range(1, 101))
    for ipdf in pdf_indices:
        nuis_PDF = Nuisance(name = f"PDF{ipdf}", type = "shape")
        for proc in processes:
            if proc.isSignal:
                nuis_PDF[proc.name] = 1.0
        nuisgroups["pdfsys"].append(nuis_PDF)
    
    # alphaS variations
    nuisgroups["alphaS"] = []
    nuis_alphaS = Nuisance(name = "alphaS", type = "shape")
    for proc in processes:
        if proc.isSignal:
            nuis_alphaS[proc.name] = 1.0
    nuisgroups["alphaS"].append(nuis_alphaS)

    #
    # writing datacards
    #
    cardname = f"{outdir}/datacard_{channel}.txt"
    if is5TeV:
        cardname = f"{outdir}/datacard_{channel}_5TeV.txt"
    WriteCard(data, processes, nuisgroups, cardname)

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

def GenerateRunCommand(card_plus: str, card_minus: str, prefix: str = "./", card_z = None):
    """
    generate the script with commands to run combine.
    inputs can probably be improved here
    """
    workdir = prefix + card_plus.rpartition('/')[0]
    card_plus = card_plus.split('/')[-1]
    card_minus = card_minus.split('/')[-1]
    if card_z:
        card_z = card_z.split('/')[-1]
    output = card_plus.replace("plus", "combined").replace(".txt", "")

    cmd = ""
    cmd += "#!/bin/bash\n\n"
    cmd += f"cd {workdir}\n"
    cmd += f"combineCards.py plus={card_plus} minus={card_minus}"
    if card_z:
        # include the z card in the combine
        cmd += f" z={card_z}"
    cmd += f"> {output}.txt\n"
    cmd += f"text2hdf5.py {output}.txt\n"
    cmd += f"combinetf.py {output}.hdf5 --binByBinStat --computeHistErrors --saveHists --doImpacts --output {output}.root\n"

    # write outputs to scripts
    with open(f"{workdir}/{output}.sh", "w") as f:
        f.write(cmd)
    os.system(f"chmod +x {workdir}/{output}.sh")

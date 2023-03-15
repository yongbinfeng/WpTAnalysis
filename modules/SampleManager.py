'''
basic class to do the sample loading, filtering, defining new 
variables and plotting.
Currently customized to work best for the W/Z case. can be easily
extened to general usage.
'''

import ROOT
import numpy as np
import sys
from collections import OrderedDict
from CMSPLOTS.myFunction import DrawHistos
from typing import List
from termcolor import colored

MINMASS = 60
MAXMASS = 120
LEPPTMIN = 25.0
LEPETA = 2.4
# LUMI = 199.27
LUMI = 200.87
LUMI_5TeV = 298.0

ROOT.ROOT.EnableImplicitMT()
ROOT.gSystem.Load("./modules/Functions_cc.so")


class DrawConfig(object):
    """
    try to sync with the arguments in DrawHistos
    in /afs/cern.ch/work/y/yofeng/public/CMSPLOTS
    """

    def __init__(self, **kwargs):
        self.xmin = kwargs.get('xmin', 0.0)
        self.xmax = kwargs.get('xmax', 0.0)
        self.xlabel = kwargs.get('xlabel', 'p_{T} [GeV]')

        self.ymin = kwargs.get('ymin', 0.5)
        self.ymax = kwargs.get('ymax', 5e5)
        self.ylabel = kwargs.get('ylabel', 'Events / GeV')

        self.dologx = kwargs.get('dologx', False)
        self.dology = kwargs.get('dology', True)
        self.dologz = kwargs.get('dologz', False)
        self.donormalize = kwargs.get('donormalize', False)
        self.donormalizebin = kwargs.get('donormalizebin', True)

        self.addOverflow = kwargs.get('addOverflow', True)
        self.addUnderflow = kwargs.get('addUnderflow', True)

        self.showratio = kwargs.get('showratio', True)
        self.ratiobase = kwargs.get('ratiobase', 1)
        # self.yrmin = kwargs.get('yrmin', 0.71)
        # self.yrmax = kwargs.get('yrmax', 1.29)
        self.yrmin = kwargs.get('yrmin', 0.86)
        self.yrmax = kwargs.get('yrmax', 1.14)
        self.yrlabel = kwargs.get('yrlabel', 'Data / MC ')

        self.legends = kwargs.get('legends', [])
        self.legendPos = kwargs.get('legendPos', [0.92, 0.88, 0.70, 0.62])
        self.lheader = kwargs.get('lheader', None)
        self.colors = kwargs.get('colors', [])

        self.redrawihist = kwargs.get('redrawihist', 0)
        self.extraText = kwargs.get('extraText', None)
        self.noCMS = kwargs.get('noCMS', True)
        self.drawoptions = kwargs.get('drawoptions', [])

        self.nMaxDigits = kwargs.get('nMaxDigits', 3)

        self.outputname = kwargs.get('outputname', 'test')
        self.savepdf = kwargs.get('savepdf', True)


class Sample(object):
    def __init__(self, inputfiles, isMC=True, xsec=1.0, color=1, legend="", name="",
                 isZSR=True, isWSR=False, nmcevt=-1, additionalnorm=1.0, is5TeV=False,
                 doTheoryVariation=True, reweightZpt=False):
        self.name = name
        self.mcxsec = xsec
        self.color = color
        self.is5TeV = is5TeV
        self.legend = legend
        self.isZSR = isZSR
        self.reweightZpt = reweightZpt
        if isWSR:
            # can not be ZSR and WSR at the same time
            self.isZSR = False
        self.isWSR = isWSR
        # assert not (self.isZSR and self.isWSR), "can not be ZSR and WSR at the same time !!"
        self.additionalnorm = additionalnorm
        self.isMC = isMC
        # this variable is used in cases where the xsec is allowed to
        # free-float, or normalized to the rest of (data-MC_sum)
        self.renormalize = False
        if isinstance(inputfiles, str):
            inputfiles = [inputfiles]
        self.initRDF(inputfiles)

        if isMC:
            if nmcevt > 0:
                self.nmcevt = nmcevt
            else:
                self.nmcevt = self.getNMCEvt(inputfiles)

        # get normalization factor
        # for data, always 1; for MC, lumi * sec / N
        self.fnorm = self.calculateNorm()

        if self.isMC and doTheoryVariation:
            # get sum of weights before selection
            # to normalize varied weights
            # 111 is the current saved number of variations
            self.nmcevts_varied = []
            self.fnorms_varied = []
            for ivar in range(111):
                nevts = self.getNMCEvtWithTheoryVariations(inputfiles, ivar)
                # sanity check
                # sometimes the varied weights are dropped in the middle
                # causing the sum of varied is zero
                if nevts == 0:
                    nevts = self.nmcevt
                self.nmcevts_varied.append(nevts)
                self.fnorms_varied.append(self.calculateNorm(nevts))

        # do preprocessing
        self.prepareVars()
        self.select()

        self.DoZptReweighting()

        self._garbagerdfs = []

    def initRDF(self, inputfiles):
        # this is a bit weird in ROOT.
        # the tree must be held by some vars, or else it would be
        # closed after this function and the rdf would crash
        self.tree = ROOT.TChain("Events")
        for ifile in inputfiles:
            for line in open(ifile, "r"):
                fname = line.rstrip()
                if fname.startswith('#'):
                    continue
                print(fname)
                self.tree.Add(fname)

        # self.rdf_org = ROOT.ROOT.RDataFrame( self.tree )
        self.rdf_temp = ROOT.ROOT.RDataFrame(self.tree)
        self.rdf_org = self.rdf_temp.Filter(
            'if(tdfentry_ == 0) {cout << "Running evtloop" << endl; } return true; ')
        # print self.rdf_org.Count().GetValue()

    def getNMCEvt(self, inputfiles):
        # print("count total number of MC events from:")
        Nevt = 0.
        for ifile in inputfiles:
            for line in open(ifile, "r"):
                fname = line.rstrip()
                if fname.startswith('#'):
                    continue
                # print(fname)
                ifile = ROOT.TFile.Open(fname)
                hist = ifile.Get("hGenWeights")
                Nevt += hist.Integral()
                ifile.Close()
                # if Nevt > 0:
                #    # the histograms saved in other files of the same dataset
                #    # is copied from the same mother histogram, so only needs the
                #    # first histogram, to get the total number of events before selection
                #    break;
        print("total number of events: {}".format(Nevt))
        return Nevt

    def getNMCEvtFromBranch(self):
        hname = "h_evtWeight_"+self.name
        self.rdf_org = self.rdf_org.Define(
            "EvtWeight", "((genWeight > 0) ? 1 : -1)")
        h_evts = self.rdf_org.Histo1D((hname, hname, 2, -2, 2), "EvtWeight")
        Npos = h_evts.GetBinContent(2)
        Nneg = h_evts.GetBinContent(1)
        print("total number of events: {}, in which {} are positive, {} are negative".format(
            Npos+Nneg, Npos, Nneg))
        return Npos-Nneg

    def getNMCEvtWithTheoryVariations(self, inputfiles: List[str], ivariation: int) -> int:
        """
        hLHEWeightSum saves the sum of weights for each theory variation
        before the selections. 
        can be used to normalize the MC histograms
        with different theory variations.
        """
        # print(f"count total number of MC events with theory variation : {ivariation}")
        Nevt = 0.
        for ifile in inputfiles:
            for line in open(ifile, "r"):
                fname = line.rstrip()
                if fname.startswith('#'):
                    continue
                ifile = ROOT.TFile.Open(fname)
                hist = ifile.Get("hLHEWeightSum")
                Nevt += hist.GetBinContent(ivariation + 1)
                ifile.Close()
        # print("total number of events: {}".format(Nevt))
        return Nevt

    def select(self):
        """ 
        Event selections.
        Not sure why but the rdf seems to have to 
        be different before and after Filter.
        So use rdf_org to hold the pre-selected one.
        """
        if self.isZSR:
            self.rdf = self.rdf_org.Filter("category==1 || category==2 || category==3") \
                .Filter("abs(lep1.Eta()) < 2.4 && abs(lep2.Eta()) < 2.4") \
                .Define("Z", "(lep1 + lep2)") \
                .Define("ZMass", "Z.M()") \
                .Filter("ZMass > {MINMASS} && ZMass < {MAXMASS}".format(MINMASS=MINMASS, MAXMASS=MAXMASS))
            # .Filter("Muon.pt[0] > {} && Muon.pt[1] > {}".format(LEPPTMIN, LEPPTMIN))  \
            # .Filter("abs(Muon.eta[0]) < {} && abs(Muon.eta[1]) < {}".format(LEPETA, LEPETA))
            # .Filter("lep_n==2")
            # .Filter("Z_pt < 40.0")
        elif self.isWSR:
            self.rdf = self.rdf_org.Filter("abs(lep.Eta())<2.4") \
                                   .Filter("lep.Pt()>25.0")
        else:
            self.rdf = self.rdf_org
        # print self.rdf.Count().GetValue()

    def ApplyCut(self, cutstring):
        rdf_postcut = self.rdf.Filter(cutstring)
        self._garbagerdfs.append(self.rdf)
        print("applied cut {} to Sample {}".format(cutstring, self.name))
        self.rdf = rdf_postcut

    def Define(self, varname, formula):
        """
        add the support to replace Sample member values (starting with self.)
        inside formula; mostly used to scale (MC) samples with their xsecs
        """
        if "self." in formula:
            pieces = formula.split()
            for i, piece in enumerate(pieces):
                if "self." not in piece:
                    continue
                try:
                    val = eval(piece)
                    pieces[i] = str(val)
                    # formula = formula.replace(piece + " ", str(val) + " ")
                except:
                    print(
                        f"Sample {self.name} does not have attribute {piece} in the string {formula}")
                    sys.exit(1)
            formula = " ".join(pieces)
        self.rdf = self.rdf.Define(varname, formula)

    def calculateNorm(self, total=-1) -> float:
        if not self.isMC:
            # data
            return 1
        lumi = LUMI_5TeV if self.is5TeV else LUMI
        if total < 0:
            total = self.nmcevt
        # print(f"lumi {lumi}, weight {lumi/total}")
        return lumi / total * self.additionalnorm

    def prepareVars(self):
        """
        sample preprocessing
        define some common variables (weights, make V, etc.)
        """
        # define weight
        # if self.isMC:
        #    self.rdf_org = self.rdf_org.Define("weight_WoVpt", "( PUWeight * NLOWeight )")
        if self.isMC:
            if self.isWSR or self.isZSR:
                self.rdf_org = self.rdf_org.Define(
                    "weight_WoVpt", f"evtWeight[0] * {self.fnorm}")
        else:
            self.rdf_org = self.rdf_org.Define("weight_WoVpt", "1.0")

        # if self.isZSR:
        #    pass
        #    ## make Z boson
        #    #self.rdf_org = self.rdf_org.Define("Zvec",  "VVec(mu_pt[0], mu_eta[0], mu_phi[0], mu_e[0], mu_pt[1], mu_eta[1], mu_phi[1], mu_e[1])")  \
        #    #                   .Define("Z_pt",  "Zvec.Pt()")    \
        #    #                   .Define("Z_eta", "Zvec.Eta()")   \
        #    #                   .Define("Z_phi", "Zvec.Phi()")   \
        #    ## make recoil vec
        #    #self.rdf_org = self.rdf_org.Define("Uvec",  "UVec(Z_pt, Z_phi, deepmet_pt, deepmet_phi)") \
        #    #                   .Define("u_pt",  "Uvec.Mod()")        \
        #    #                   .Define("u_phi", "Uvec.Phi()")        \
        #    #                   .Define("u1",    "u_pt * TMath::Cos(u_phi + TMath::Pi() - Z_phi)") \
        #    #                   .Define("u2",    "u_pt * TMath::Sin(u_phi + TMath::Pi() - Z_phi)") \
        #    ## z pt reweight
        #    #self.rdf_org = self.rdf_org.Define("ZptWeight", "ZptReWeight(Z_pt, h_zpt_ratio, isData)" if self.reweightZpt else "1.")

        #    #self.rdf_org = self.rdf_org.Define("weight", "ZptWeight * weight_WoVpt")
        # elif self.isWSR and self.isMC:
        #    ## make W boson in the transverse plane
        #    #self.rdf_org = self.rdf_org.Define("Wvec", "VVec(mu_pt[0], mu_eta[0], mu_phi[0], mu_e[0], gen_met_pt, 0., gen_met_phi, gen_met_pt)") \
        #    #                           .Define("W_pt",  "Wvec.Pt()") \
        #    #                           .Define("W_phi", "Wvec.Phi()")
        #    #self.rdf_org = self.rdf_org.Define("W_pt", "trueW_pt[0]") \
        #    #                           .Define("W_phi", "trueW_phi[0]")
        #    # make recoil vec
        #    #self.rdf_org = self.rdf_org.Define("Uvec", "UVec(W_pt, W_phi, deepmet_pt, deepmet_phi)") \
        #    #                           .Define("u_pt",  "Uvec.Mod()")   \
        #    #                           .Define("u_phi", "Uvec.Phi()")   \
        #    # weights
        #    self.rdf_org = self.rdf_org.Define("weight", "weight_WoVpt")

        # else:
        # self.rdf_org = self.rdf_org.Define("weight", "weight_WoVpt")
    def DoZptReweighting(self):
        if self.reweightZpt:
            print(
                colored(f"Apply Zpt Reweighting to Sample {self.name}", "red"))
            self.rdf = self.rdf.Define("Z_pt", "genV.Pt()")\
                               .Define("ZptWeight_Bare", "ZptReWeight(Z_pt, h_zpt_ratio, 0.)") \
                               .Define("evtWeight_Bare", "evtWeight[0]") \
                               .Define("evtWeight_Bare_WVpt", "evtWeight_Bare * ZptWeight_Bare") \
                               .Define("dumbVal", "1.0")
            # need to renormalize the weights
            self.htemp_weighted = self.rdf.Histo1D(
                ("_htemp_weighted", "_htemp_weighted", 2, 0, 2.0), "dumbVal", "evtWeight_Bare_WVpt")
            self.htemp_woweight = self.rdf.Histo1D(
                ("_htemp_woweight", "_htemp_woweight", 2, 0, 2.0), "dumbVal", "evtWeight_Bare")
            self.zptrenorm = self.htemp_woweight.Integral() / self.htemp_weighted.Integral()
            print("Without Zpt Reweighting: ", self.htemp_woweight.Integral())
            print("With Zpt Reweighting: ", self.htemp_weighted.Integral())
            print("Zpt renormalization factor: ", self.zptrenorm)
            self.rdf = self.rdf.Define(
                "ZptWeight", f"ZptWeight_Bare * {self.zptrenorm}")
        else:
            self.rdf = self.rdf.Define("ZptWeight", "1.")
        self.rdf = self.rdf.Define("weight_WVpt", "ZptWeight * weight_WoVpt")


class SampleManager(object):
    def __init__(self, data, mcs=[], is5TeV=False, outdir="plots"):
        self.data = data
        self.mcs = mcs
        self.to_draw = OrderedDict()
        self.to_draw2D = OrderedDict()
        self.initGroups()

        self.is5TeV = is5TeV
        # sanity check:
        # if samplemanager is at 5Tev, then all MC samples should be at 5Tev
        # also true for the opposite
        for mc in self.mcs:
            assert self.is5TeV == mc.is5TeV, "for MC {} is at 5TeV? {}, inconsistent with the samplemanager {}".format(
                mc.name, mc.is5TeV, self.is5TeV)

        # output directory for plots
        self.outdir = outdir

        # collection to save all histograms, stacks, etc, for further uses
        # For future improvements, some of them are sample-related,
        # so can probably go to Sample class
        self.hdatas = {}
        self.hsmcs = {}
        self.hratios = {}
        self.hmcs = {}
        
        self.hdatas2D = {}
        self.hmcs2D = {}
        self.hsubtracted2D = {}

        # count the number of events in data and MC
        self.counts = []

    def ApplyCutAll(self, cutstring):
        self.data.ApplyCut(cutstring)
        for mc in self.mcs:
            mc.ApplyCut(cutstring)

    def ApplyCutSpecificMCs(self, cutstring, sampnames=[], sampgroupnames=[]):
        for mc in self.mcs:
            if mc.name in sampnames:
                print("Apply Cut {} for sample {}".format(cutstring, mc.name))
                mc.ApplyCut(cutstring)
            elif mc.groupname in sampgroupnames:
                print("Apply Cut {} for sample {} in grouped sample {}".format(
                    cutstring, mc.name, mc.groupname))
                mc.ApplyCut(cutstring)

    def DefineAll(self, varname, formula, excludes=[], excludeGroups=[]):
        self.DefineData(varname, formula)
        self.DefineMC(varname, formula, excludes, excludeGroups)

    def DefineMC(self, varname, formula, excludes=[], excludeGroups=[]):
        for mc in self.mcs:
            if mc.name in excludes:
                print("Define {} with {}. Skip sample {}".format(
                    varname, formula, mc.name))
                continue
            if mc.groupname in excludeGroups:
                print("Define {} with {}. Skip sample {} in grouped sample {}".format(
                    varname, formula, mc.name, mc.groupname))
                continue
            mc.Define(varname, formula)

    def DefineSpecificMCs(self, varname, formula, sampnames=[], sampgroupnames=[]):
        for mc in self.mcs:
            if mc.name in sampnames:
                print("Define {} with {} for sample {}".format(
                    varname, formula, mc.name))
                mc.Define(varname, formula)
            elif mc.groupname in sampgroupnames:
                print("Define {} with {} for sample {} in grouped sample {}".format(
                    varname, formula, mc.name, mc.groupname))
                mc.Define(varname, formula)

    def DefineData(self, varname, formula):
        self.data.Define(varname, formula)

    def initGroups(self):
        # default group information is the same as the sample itself
        for mc in self.mcs:
            mc.groupname = mc.name
            mc.groupcolor = mc.color
            mc.grouplegend = mc.legend

    def groupMCs(self, mcnames_to_be_grouped, groupname, groupcolor, grouplegend, renormalize=False):
        for mc in self.mcs:
            if mc.name in mcnames_to_be_grouped:
                print("Group sample ", mc.name, " to ", groupname)
                mc.groupname = groupname
                mc.groupcolor = groupcolor
                mc.grouplegend = grouplegend
                if renormalize:
                    mc.renormalize = True
                mcnames_to_be_grouped.remove(mc.name)

        if mcnames_to_be_grouped:
            print(colored('some somples are not grouped yet, please check the names...: {}'.format(
                mcnames_to_be_grouped), 'red'))

    def getMCByName(self, mcname):
        for mc in self.mcs:
            if mc.name == mcname:
                return mc
        raise Exception('can not find the mcname {}'.format(mcname))
    
    def cacheDraw2D(self, *args, **kwds):
        if len(args) == 10 or len(args) == 11:
            self.cacheDraw2D_fb(*args, **kwds)
        elif len(args) == 6 or len(args) == 7:
            self.cacheDraw2D_vb(*args, **kwds)
        else:
            raise ValueError(
                'The argument is problematic. It has to be either 6 or 10.')
    
    def cacheDraw2D_fb(self, varname1, varname2, hname, nbins1, xmin1, xmax1, nbins2, xmin2, xmax2, drawconfigs, weightname="weight_WVpt"):
        """ 
        cache the var to be drawn.
        But do not launch the action by the 'lazy action' in RDataFrame
        """
        h_data = self.data.rdf.Histo2D(
            (hname+"_data", hname, nbins1, xmin1, xmax1, nbins2, xmin2, xmax2), varname1, varname2, weightname)
        h_mcs = []
        for imc in range(len(self.mcs)):
            mc = self.mcs[imc]
            h_mcs.append(mc.rdf.Histo2D(
                (hname+"_mc{}".format(imc), hname, nbins1, xmin1, xmax1, nbins2, xmin2, xmax2), varname1, varname2, weightname))
        self.to_draw2D[hname] = (h_data, h_mcs, drawconfigs)
        
    def cacheDraw2D_vb(self, varname1, varname2, hname, xbins, ybins, drawconfigs, weightname="weight_WVPt"):
        assert isinstance(xbins, np.ndarray) and isinstance(ybins, np.ndarray), "input must be np array"
        nbinsx = xbins.size - 1
        nbinsy = ybins.size - 1
        h_data = self.data.rdf.Histo2D(
            (hname+"_data", hname, nbinsx, xbins, nbinsy, ybins), varname1, varname2, weightname)
        h_mcs = []
        for imc in range(len(self.mcs)):
            mc = self.mcs[imc]
            h_mcs.append(mc.rdf.Histo2D(
                (hname+"_mc_"+str(imc), hname, nbinsx, xbins, nbinsy, ybins), varname1, varname2, weightname))

        self.to_draw[hname] = (h_data, h_mcs, drawconfigs)
        

    def cacheDraw(self, *args, **kwds):
        if len(args) == 6 or len(args) == 7:
            self.cacheDraw_fb(*args, **kwds)
        elif len(args) == 4 or len(args) == 5:
            self.cacheDraw_vb(*args, **kwds)
        else:
            raise ValueError(
                'The argument is problematic. It has to be either 4 or 6')

    def cacheDraw_fb(self, varname, hname, nbins, xmin, xmax, drawconfigs, weightname="weight_WVpt"):
        """ 
        cache the var to be drawn.
        But do not launch the action by the 'lazy action' in RDataFrame
        """
        h_data = self.data.rdf.Histo1D(
            (hname+"_data", hname, nbins, xmin, xmax), varname, weightname)
        h_mcs = []
        for imc in range(len(self.mcs)):
            mc = self.mcs[imc]
            h_mcs.append(mc.rdf.Histo1D(
                (hname+"_mc_"+str(imc), hname, nbins, xmin, xmax), varname, weightname))

        self.to_draw[hname] = (h_data, h_mcs, drawconfigs)

    def cacheDraw_vb(self, varname, hname, xbins, drawconfigs, weightname="weight_WVpt"):
        """
        cache the var to be drawn.
        But do not launch the action by the 'lazy action' in RDataFrame
        """
        assert isinstance(xbins, np.ndarray), "input must be np array"
        nbins = xbins.size - 1
        h_data = self.data.rdf.Histo1D(
            (hname+"_data", hname, nbins, xbins), varname, weightname)
        h_mcs = []
        for imc in range(len(self.mcs)):
            mc = self.mcs[imc]
            h_mcs.append(mc.rdf.Histo1D(
                (hname+"_mc_"+str(imc), hname, nbins, xbins), varname, weightname))

        self.to_draw[hname] = (h_data, h_mcs, drawconfigs)

    def _DrawPlot(self, h_data, h_mcs, drawconfigs, hname, is2D=False):
        legends = []

        h_data.Scale(1.0)
        h_data.SetLineColor(1)
        h_data.SetMarkerStyle(20)
        h_data.SetMarkerSize(1)
        if drawconfigs.donormalizebin:
            h_data.Scale(1.0, "width")
        legends.append("Data")

        # count the number of events only need to be done once
        docounting = (len(self.counts) == 0)
        intgoption = "width" if drawconfigs.donormalizebin == "width" else ""

        if docounting:
            self.counts.append(h_data.Integral(
                0, h_data.GetNbinsX()+1, intgoption))

        hgroupedmcs = OrderedDict()
        hmcs = {}
        group_to_renormalize = ""
        for imc in range(len(h_mcs)):
            if drawconfigs.donormalizebin:
                h_mcs[imc].Scale(1.0, "width")

            if docounting:
                self.counts.append(h_mcs[imc].Integral(
                    0, h_mcs[imc].GetNbinsX()+1, intgoption))

            groupname = self.mcs[imc].groupname
            if groupname not in hgroupedmcs:
                hgroupedmcs[groupname] = h_mcs[imc].Clone(
                    h_mcs[imc].GetName().replace("_mc_", "_grouped_"+groupname))
                hgroupedmcs[groupname].SetLineColor(self.mcs[imc].groupcolor)
                hgroupedmcs[groupname].SetFillColor(self.mcs[imc].groupcolor)

                legends.append(self.mcs[imc].grouplegend)

                if self.mcs[imc].renormalize:
                    group_to_renormalize = groupname
            else:
                # group mc already exist. Add to the histogram
                hgroupedmcs[groupname].Add(h_mcs[imc].GetValue())

            # save the histogram to collection
            hmcs[self.mcs[imc].name] = h_mcs[imc].Clone(
                h_mcs[imc].GetName().replace("_mc_", "_"+self.mcs[imc].name))

        if group_to_renormalize:
            ndata = h_data.Integral(0, h_data.GetNbinsX()+1, intgoption)
            nmc = 0
            for gname, ghisto in hgroupedmcs.items():
                if gname != group_to_renormalize:
                    nmc += ghisto.Integral(0, ghisto.GetNbinsX()+1, intgoption)
            weight = float(ndata-nmc) / hgroupedmcs[group_to_renormalize].Integral(
                0, hgroupedmcs[group_to_renormalize].GetNbinsX()+1)
            print("Renormalize group for {}, with {} data, {} MC and weight {}".format(
                group_to_renormalize, ndata, nmc, weight))
            hgroupedmcs[group_to_renormalize].Scale(weight)
            
        if not drawconfigs.legends:
            drawconfigs.legends = legends
        if drawconfigs.outputname == "test":
            # change default outputname to the histo name
            drawconfigs.outputname = hname

        if not is2D:
            hsname = "hs_" + drawconfigs.outputname
            hs_gmc = ROOT.THStack(hsname, hsname)
            n_mcs = len(hgroupedmcs)
            for h_gmc in reversed(list(hgroupedmcs.values())):
                # if drawconfigs.donormalizebin:
                #    h_gmc.Scale(1.0, "width")
                hs_gmc.Add(h_gmc)

            self.hdatas[drawconfigs.outputname] = h_data
            self.hsmcs[drawconfigs.outputname] = hs_gmc
            self.hmcs[drawconfigs.outputname] = hmcs
            if n_mcs > 0:
                h_to_draw = [h_data, hs_gmc]
            else:
                h_to_draw = [h_data]
        else:
            # subtract MC from data for 2D plot
            h_MC = None
            h_subtract = h_data.Clone(h_data.GetName() + "_DataForPlotting")
            for h_gmc in reversed(list(hgroupedmcs.values())):
                if h_MC is None:
                    h_MC = h_gmc.Clone(h_gmc.GetName() + "_AllMCCombined")
                else:
                    h_MC.Add(h_gmc)
            if h_MC:
                h_subtract.Add(h_MC, -1)
            
            self.hdatas2D[drawconfigs.outputname] = h_data
            self.hmcs2D[drawconfigs.outputname] = h_MC
            self.hsubtracted2D[drawconfigs.outputname] = h_subtract
            
            h_to_draw = [h_subtract]
            drawconfigs.legends = []
            drawconfigs.drawoptions = ["colz,text"]
            drawconfigs.dologz = True
            
        self.hratios[drawconfigs.outputname] = DrawHistos(h_to_draw, drawconfigs.legends, drawconfigs.xmin, drawconfigs.xmax, drawconfigs.xlabel, drawconfigs.ymin, drawconfigs.ymax, drawconfigs.ylabel, drawconfigs.outputname, dology=drawconfigs.dology, dologx=drawconfigs.dologx, showratio=drawconfigs.showratio, yrmax=drawconfigs.yrmax, yrmin=drawconfigs.yrmin, yrlabel=drawconfigs.yrlabel, donormalize=drawconfigs.donormalize, ratiobase=drawconfigs.ratiobase, legendPos=drawconfigs.legendPos, redrawihist=drawconfigs.redrawihist, extraText=drawconfigs.extraText, noCMS=drawconfigs.noCMS, addOverflow=drawconfigs.addOverflow, addUnderflow=drawconfigs.addUnderflow, nMaxDigits=drawconfigs.nMaxDigits, is5TeV=self.is5TeV, lheader=drawconfigs.lheader, outdir=self.outdir, drawoptions = drawconfigs.drawoptions, dologz = drawconfigs.dologz, doth2=is2D, savepdf=drawconfigs.savepdf)
            

    def launchDraw(self):
        """
        launch all the draw options in to_draw
        """
        print("starting drawing")
        for hname, plotinfo in self.to_draw.items():
            self._DrawPlot(plotinfo[0], plotinfo[1], plotinfo[2], hname)
        for hname, plotinfo in self.to_draw2D.items():
            self._DrawPlot(plotinfo[0], plotinfo[1], plotinfo[2], hname, is2D=True)
        self.to_draw.clear()
        print("finished drawing..")

    def dumpCounts(self):
        """
        dump the number of events in data and MC
        """
        print("dump the counts..")
        if len(self.counts) == 0:
            print("draw at least one distribution and then count")
            return

        print((self.counts))

        print("Data %d" % self.counts[0])
        for imc in range(len(self.mcs)):
            print("MC %s count %f" % (self.mcs[imc].name, self.counts[imc+1]))

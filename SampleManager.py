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
sys.path.append("/afs/cern.ch/work/y/yofeng/public/CMSPLOTS")
from myFunction import DrawHistos

MINMASS = 60
MAXMASS = 120
LEPPTMIN = 25.0
LEPETA = 2.4
LUMI = 199.27

ROOT.ROOT.EnableImplicitMT()
ROOT.gSystem.Load("/afs/cern.ch/work/y/yofeng/public/WpT/CMSSW_10_6_0/src/PostCorrNTuple/Functions_cc.so")


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
        self.donormalize = kwargs.get('donormalize', False)
        self.donormalizebin = kwargs.get('donormalizebin', True)

        self.addOverflow = kwargs.get('addOverflow', True)
        self.addUnderflow = kwargs.get('addUnderflow', True)

        self.showratio = kwargs.get('showratio', True)
        self.ratiobase = kwargs.get('ratiobase', 1)
        #self.yrmin = kwargs.get('yrmin', 0.71)
        #self.yrmax = kwargs.get('yrmax', 1.29)
        self.yrmin = kwargs.get('yrmin', 0.86)
        self.yrmax = kwargs.get('yrmax', 1.14)
        self.yrlabel = kwargs.get('yrlabel', 'Data / MC ')

        self.legends = kwargs.get('legends', [])
        self.legendPos = kwargs.get('legendPos', [0.92, 0.88, 0.70, 0.62])
        self.colors = kwargs.get('colors', [])

        self.redrawihist = kwargs.get('redrawihist', 0)
        self.extraText = kwargs.get('extraText', None)
        self.noCMS = kwargs.get('noCMS', True)
        
        self.nMaxDigits = kwargs.get('nMaxDigits', 3)

        self.outputname = kwargs.get('outputname', 'test')


class Sample(object):
    def __init__(self, inputfiles, isMC=True, xsec=1.0, color=1, reweightzpt = False, legend="", name="", 
                 isZSR=True, isWSR=False, applySF=True, bjetVeto=False, nmcevt=-1, additionalnorm=1.0):
        if isMC:
            if nmcevt>0:
                self.nmcevt = nmcevt
            else:
                self.nmcevt = self.getNMCEvt(inputfiles)
            self.nmcevt = self.nmcevt
            self.mcxsec = xsec
            self.color = color
            #self.normfactor = LUMI * self.mcxsec / self.nmcevt
            self.normfactor = 1.0
            self.reweightzpt = reweightzpt
        else:
            self.nmcevt = 0
            self.mcxsec = 0
            self.color = 1
            self.normfactor = 1.0
            self.reweightzpt = False
        self.additionalnorm = additionalnorm

        self.isMC = isMC
        self.legend = legend
        self.name = name
        self.isZSR = isZSR
        if isWSR:
            # can not be ZSR and WSR at the same time
            self.isZSR = False
        self.isWSR = isWSR
        #assert not (self.isZSR and self.isWSR), "can not be ZSR and WSR at the same time !!"
        self.applySF = applySF
        # this variable is used in cases where the xsec is allowed to 
        # free-float, or normalized to the rest of (data-MC_sum)
        self.renormalize = False
        self.renormalizefactor = 1.0
        self.bjetVeto = bjetVeto
        self.initRDF(inputfiles)

        self.prepareVars()
        self.select()

        self._garbagerdfs = []

    def initRDF(self, inputfiles):
        # this is a bit weird in ROOT.
        # the tree must be held by some vars, or else it would be 
        # closed after this function and the rdf would crash
        self.tree = ROOT.TChain("Events")
        for line in open( inputfiles, "r"):
            fname = line.rstrip()
            if fname.startswith('#'):
                continue
            print fname
            self.tree.Add( fname )

        self.rdf_org = ROOT.ROOT.RDataFrame( self.tree )
        #print self.rdf_org.Count().GetValue()

    def getNMCEvt(self, inputfiles):
        print "count total number of MC events from:"
        for line in open( inputfiles, "r"):
            fname = line.rstrip()
            if fname.startswith('#'):
                continue
            print fname
            ifile = ROOT.TFile(fname)
            hist = ifile.Get("hGenWeights")
            Nevt = hist.Integral()
            ifile.Close()
        print "total number of events: {}".format(Nevt)
        return Nevt

    def getNMCEvtFromBranch(self):
        hname = "h_evtWeight_"+self.name
        self.rdf_org = self.rdf_org.Define("EvtWeight", "((genWeight > 0) ? 1 : -1)")
        h_evts = self.rdf_org.Histo1D( (hname, hname, 2, -2, 2), "EvtWeight")
        Npos = h_evts.GetBinContent(2)
        Nneg = h_evts.GetBinContent(1)
        print "total number of events: {}, in which {} are positive, {} are negative".format(Npos+Nneg, Npos, Nneg)
        return Npos-Nneg

    def select(self):
        """ 
        Event selections.
        Not sure why but the rdf seems to have to 
        be different before and after Filter.
        So use rdf_org to hold the pre-selected one.
        """
        if self.isZSR:
            self.rdf_temp = self.rdf_org.Filter("category==1 || category==2 || category==3") \
                                .Filter("abs(lep1.Eta()) < 2.4 && abs(lep2.Eta()) < 2.4") \
                                .Define("lep1_corr", "VLep(lep1_corrected_pt, lep1.Eta(), lep1.Phi(), lep1.M())") \
                                .Define("lep2_corr", "VLep(lep2_corrected_pt, lep2.Eta(), lep2.Phi(), lep2.M())") \
                                .Define("Z", "(lep1_corr + lep2_corr)") \
                                .Define("ZMass", "Z.M()") \
                                .Filter("ZMass > {MINMASS} && ZMass < {MAXMASS}".format(MINMASS=MINMASS, MAXMASS=MAXMASS))
                                #.Filter("Muon.pt[0] > {} && Muon.pt[1] > {}".format(LEPPTMIN, LEPPTMIN))  \
                                #.Filter("abs(Muon.eta[0]) < {} && abs(Muon.eta[1]) < {}".format(LEPETA, LEPETA))
                                #.Filter("lep_n==2")
                                #.Filter("Z_pt < 40.0")
            if self.bjetVeto:
                self.rdf = self.rdf_temp.Filter("jet_CSVLoose_n<1")
            else:
                self.rdf = self.rdf_temp
        elif self.isWSR:
            self.rdf = self.rdf_org.Filter("abs(lep.Eta())<2.4") \
                                   .Filter("lep.Pt()>25.0") \
                                   .Filter("mtCorr>40")
        else:
            self.rdf = self.rdf_org
        #print self.rdf.Count().GetValue()

    def ApplyCut(self, cutstring):
        rdf_postcut = self.rdf.Filter( cutstring )
        self._garbagerdfs.append( self.rdf )
        print "applied cut {} to Sample {}".format(cutstring, self.name)
        self.rdf = rdf_postcut

    def Define(self, varname, formula):
        self.rdf = self.rdf.Define(varname, formula)

    def prepareVars(self):
        # define weight
        #if self.isMC and self.applySF:
        #    self.rdf_org = self.rdf_org.Define("weight_WoVpt", "( PUWeight * NLOWeight * mu_trigSF * mu_isoSF * mu_trkSF * mu_idSF )")
        #elif self.isMC:
        #    self.rdf_org = self.rdf_org.Define("weight_WoVpt", "( PUWeight * NLOWeight )")
        if self.isMC:
            print(str(LUMI/self.nmcevt))
            self.rdf_org = self.rdf_org.Define("mcnorm", str(LUMI/self.nmcevt * self.additionalnorm))
            if self.isZSR:
                self.rdf_org = self.rdf_org.Define("weight_WoVpt_WoEff", "scale1fb * prefireWeight * mcnorm") \
                                       .Define("weight_WoVpt", "scale1fb * prefireWeight * mcnorm * lepsfweight")
            elif self.isWSR:
                self.rdf_org = self.rdf_org.Define("weight_WoVpt", "evtWeight[0] * mcnorm")
        else:
            self.rdf_org = self.rdf_org.Define("weight_WoVpt", str(self.additionalnorm))

        #if self.isZSR:
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
        #    #self.rdf_org = self.rdf_org.Define("ZptWeight", "ZptReWeight(Z_pt, h_zpt_ratio, isData)" if self.reweightzpt else "1.")

        #    #self.rdf_org = self.rdf_org.Define("weight", "ZptWeight * weight_WoVpt")
        #elif self.isWSR and self.isMC:
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

        #else:
        #self.rdf_org = self.rdf_org.Define("weight", "weight_WoVpt")


class SampleManager(object):
    def __init__(self, data, mcs=[]):
        self.data = data
        self.mcs = mcs
        self.to_draw = OrderedDict()
        self.getNormFactors()
        self.initGroups()

        self.hdatas  = {}
        self.hsmcs   = {}
        self.hratios = {}

        # count the number of events in data and MC
        self.counts = []

    def ApplyCutAll(self, cutstring):
        self.data.ApplyCut(cutstring)
        for mc in self.mcs:
            mc.ApplyCut(cutstring)

    def ApplyCutSpecificMCs(self, cutstring, sampnames=[], sampgroupnames=[]):
        for mc in self.mcs:
            if mc.name in sampnames:
                print "Apply Cut {} for sample {}".format(cutstring, mc.name)
                mc.ApplyCut(cutstring)
            elif mc.groupname in sampgroupnames:
                print "Apply Cut {} for sample {} in grouped sample {}".format(cutstring, mc.name, mc.groupname)
                mc.ApplyCut(cutstring)

    def DefineAll(self, varname, formula, excludes=[], excludeGroups=[]):
        self.DefineData(varname, formula)
        self.DefineMC(  varname, formula, excludes, excludeGroups)

    def DefineMC(self, varname, formula, excludes=[], excludeGroups=[]):
        for mc in self.mcs:
            if mc.name in excludes:
                print "Define {} with {}. Skip sample {}".format(varname, formula, mc.name)
                continue
            if mc.groupname in excludeGroups:
                print "Define {} with {}. Skip sample {} in grouped sample {}".format(varname, formula, mc.name, mc.groupname)
                continue
            mc.rdf = mc.rdf.Define(varname, formula)

    def DefineSpecificMCs(self, varname, formula, sampnames=[], sampgroupnames=[]):
        for mc in self.mcs:
            if mc.name in sampnames:
                print "Define {} with {} for sample {}".format(varname, formula, mc.name)
                mc.rdf = mc.rdf.Define(varname, formula)
            elif mc.groupname in sampgroupnames:
                print "Define {} with {} for sample {} in grouped sample {}".format(varname, formula, mc.name, mc.groupname)
                mc.rdf = mc.rdf.Define(varname, formula)

    def DefineData(self, varname, formula):
        self.data.rdf = self.data.rdf.Define(varname, formula)

    def getNormFactors(self):
        self.normfactors = []
        for mc in self.mcs:
            self.normfactors.append( mc.normfactor )
        print "normalization factors ", self.normfactors

    def initGroups(self):
        # default group information is the same as the sample itself
        for mc in self.mcs:
            mc.groupname = mc.name
            mc.groupcolor = mc.color
            mc.grouplegend = mc.legend

    def groupMCs(self, mcnames_to_be_grouped, groupname, groupcolor, grouplegend, renormalize=False, renormalizefactor=1.0):
        for mc in self.mcs:
            if mc.name in mcnames_to_be_grouped:
                print "Group sample ", mc.name, " to ", groupname
                mc.groupname = groupname
                mc.groupcolor = groupcolor
                mc.grouplegend = grouplegend
                if renormalize:
                    mc.renormalize = True
                mc.renormalizefactor = renormalizefactor
                mcnames_to_be_grouped.remove(mc.name)

        if mcnames_to_be_grouped:
            from termcolor import colored
            print colored('some somples are not grouped yet, please check the names...: {}'.format(mcnames_to_be_grouped), 'red')

    def cacheDraw(self, *args, **kwds):
        if len(args) == 6:
            self.cacheDraw_fb(*args, **kwds)
        elif len(args) == 4:   
            self.cacheDraw_vb(*args, **kwds)
        else:
            raise ValueError('The argument is problematic. It has to be either 4 or 6')

    def cacheDraw_fb(self, varname, hname, nbins, xmin, xmax, drawconfigs, weightname="weight_WoVpt"):
        """ 
        cache the var to be drawn.
        But do not launch the action by the 'lazy action' in RDataFrame
        """
        h_data = self.data.rdf.Histo1D( (hname+"_data", hname, nbins, xmin, xmax), varname, weightname )
        h_mcs = []
        for imc in xrange(len(self.mcs)):
            mc = self.mcs[imc]
            h_mcs.append( mc.rdf.Histo1D( (hname+"_mc_"+str(imc), hname, nbins, xmin, xmax), varname, weightname) )

        self.to_draw[hname] = (h_data, h_mcs, drawconfigs)

    def cacheDraw_vb(self, varname, hname, xbins, drawconfigs, weightname="weight_WoVpt"):
        """
        cache the var to be drawn.
        But do not launch the action by the 'lazy action' in RDataFrame
        """
        assert isinstance(xbins, np.ndarray), "input must be np array"
        nbins = xbins.size - 1
        h_data = self.data.rdf.Histo1D( (hname+"_data", hname, nbins, xbins), varname, weightname )
        h_mcs = []
        for imc in xrange(len(self.mcs)):
            mc = self.mcs[imc]
            h_mcs.append( mc.rdf.Histo1D( (hname+"_mc_"+str(imc), hname, nbins, xbins), varname, weightname) )

        self.to_draw[hname] = (h_data, h_mcs, drawconfigs)

    def _DrawPlot(self, h_data, h_mcs, drawconfigs, hname):
        legends = []

        h_data.Scale(1.0)
        h_data.SetLineColor(1)
        h_data.SetMarkerStyle(20)
        h_data.SetMarkerSize(1)
        legends.append("Data")

        # count the number of events only need to be done once
        docounting = (len(self.counts)==0)

        if docounting:
            self.counts.append(h_data.Integral(0, h_data.GetNbinsX()+1))

        hgroupedmcs = OrderedDict()
        group_to_renormalize = ""
        for imc in xrange(len(h_mcs)):
            # scale the MC to the xsec
            h_mcs[imc].Scale( self.normfactors[imc] )
            if self.mcs[imc].renormalizefactor!=1.0:
                print "renormalize MC {} with a factor or {}".format(self.mcs[imc].name, self.mcs[imc].renormalizefactor)
                h_mcs[imc].Scale( self.mcs[imc].renormalizefactor )

            if docounting:
                self.counts.append( h_mcs[imc].Integral(0, h_mcs[imc].GetNbinsX()+1) )
    
            groupname = self.mcs[imc].groupname
            if groupname not in hgroupedmcs:
                hgroupedmcs[groupname] = h_mcs[imc].Clone(h_mcs[imc].GetName().replace("_mc_","_grouped_"+groupname))
                hgroupedmcs[groupname].SetLineColor( self.mcs[imc].groupcolor )
                hgroupedmcs[groupname].SetFillColor( self.mcs[imc].groupcolor )

                legends.append( self.mcs[imc].grouplegend )

                if self.mcs[imc].renormalize:
                    group_to_renormalize = groupname
            else:
                # group mc already exist. Add to the histogram
                hgroupedmcs[groupname].Add( h_mcs[imc].GetValue() )

        if group_to_renormalize:
            ndata = h_data.Integral(0, h_data.GetNbinsX()+1)
            nmc=0
            for gname, ghisto in hgroupedmcs.iteritems():
                if gname!=group_to_renormalize:
                    nmc += ghisto.Integral(0, ghisto.GetNbinsX()+1)
            weight = float(ndata-nmc) / hgroupedmcs[group_to_renormalize].Integral(0, hgroupedmcs[group_to_renormalize].GetNbinsX()+1)
            print "Renormalize group for {}, with {} data, {} MC and weight {}".format(group_to_renormalize, ndata, nmc, weight)
            hgroupedmcs[group_to_renormalize].Scale(weight)

        hsname = "hs_" + drawconfigs.outputname
        hs_gmc = ROOT.THStack( hsname, hsname)
        for h_gmc in reversed(hgroupedmcs.values()):
            if drawconfigs.donormalizebin:
                # scale to bin width
                h_gmc.Scale(1.0, "width")
            hs_gmc.Add( h_gmc )

        if not drawconfigs.legends:
            drawconfigs.legends = legends

        if drawconfigs.outputname == "test":
            # change default outputname to the histo name
            drawconfigs.outputname = hname

        if drawconfigs.donormalizebin:
            h_data.Scale(1.0, "width")
        
        self.hdatas[ drawconfigs.outputname] = h_data
        self.hsmcs[  drawconfigs.outputname] = hs_gmc
        self.hratios[drawconfigs.outputname] = DrawHistos( [h_data, hs_gmc], drawconfigs.legends, drawconfigs.xmin, drawconfigs.xmax, drawconfigs.xlabel, drawconfigs.ymin, drawconfigs.ymax, drawconfigs.ylabel, drawconfigs.outputname, dology=drawconfigs.dology, dologx=drawconfigs.dologx, showratio=drawconfigs.showratio, yrmax = drawconfigs.yrmax, yrmin = drawconfigs.yrmin, yrlabel = drawconfigs.yrlabel, donormalize=drawconfigs.donormalize, ratiobase=drawconfigs.ratiobase, legendPos = drawconfigs.legendPos, redrawihist = drawconfigs.redrawihist, extraText = drawconfigs.extraText, noCMS = drawconfigs.noCMS, addOverflow = drawconfigs.addOverflow, addUnderflow = drawconfigs.addUnderflow, nMaxDigits = drawconfigs.nMaxDigits)
        

    def launchDraw(self):
        """
        launch all the draw options in to_draw
        """
        print "starting drawing"
        for hname, plotinfo in self.to_draw.iteritems():
            self._DrawPlot( plotinfo[0], plotinfo[1], plotinfo[2], hname)
        self.to_draw.clear()
        print "finished drawing.."

    def dumpCounts(self):
        """
        dump the number of events in data and MC
        """
        print "dump the counts.."
        if len(self.counts)==0:
            print "draw at least one distribution and then count"
            return

        print(self.counts)

        print "Data %d"%self.counts[0]
        for imc in xrange(len(self.mcs)):
            print "MC %s count %f"%(self.mcs[imc].name, self.counts[imc+1])

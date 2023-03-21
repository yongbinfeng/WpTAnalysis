"""
script to measure FR and carry out closure test
"""
import ROOT
import numpy as np
import sys
import argparse
from collections import OrderedDict
from modules.Binnings import mass_bins_forqcd
from CMSPLOTS.myFunction import DrawHistos, MultiplyH2, PositiveProtection, IntegralAndError2D, CombineOneBin2D, GetRatioPanel, LHistos2Hist, IncludeOverflow2D, TH2ToTH1s, SymmetrizeHisto

ROOT.gROOT.SetBatch(True)

def DrawDataMCStack(hdata, hmc, hpred = None, xmin = 0, xmax = 140, xlabel = "m_{T} GeV", outputname = "test", is5TeV = False):
    hdata.SetLineColor(1)
    hdata.SetMarkerColor(1)
    hmc.SetLineColor(92)
    hmc.SetFillColor(92)
    hmc.SetMarkerColor(92)
    hsmc = ROOT.THStack(f"hsmc_{outputname}", f"hsmc_{outputname}")
    labels = ["Obs", "MC"]
    if hpred:
        hpred.SetLineColor(226)
        hpred.SetMarkerColor(226)
        hpred.SetFillColor(226)
        hsmc.Add(hpred)
        labels.append("Pred QCD")
    # add hmc later than hpred, as first added will be drawn at the bottom
    hsmc.Add(hmc)
    hratiopanel = GetRatioPanel(hsmc) 
    DrawHistos([hdata, hsmc], labels, xmin, xmax, xlabel, 0., None, "Events", f"{outputname}", dology=False, showratio=True, yrmin = 0.89, yrmax = 1.11, ratiobase=1, redrawihist=0, hratiopanel=hratiopanel, is5TeV=is5TeV)
    
def GetFR(hnums_data, hnums_mc, hdens_data, hdens_mc, projX = False, projY = False, avergeEta = False, suffix = "Clone", scaleNumMC = 1.0, scaleDenMC = 1.0):
    hnum_data = LHistos2Hist(hnums_data, f"{hnums_data[0].GetName()}_{suffix}")
    hnum_mc = LHistos2Hist(hnums_mc, f"{hnums_mc[0].GetName()}_{suffix}")
    hden_data = LHistos2Hist(hdens_data, f"{hdens_data[0].GetName()}_{suffix}")
    hden_mc = LHistos2Hist(hdens_mc, f"{hdens_mc[0].GetName()}_{suffix}")
    
    # scale up/down contributions from MC in the (anti)isolated region
    # for systematic variations
    hnum_mc.Scale(scaleNumMC)
    hden_mc.Scale(scaleDenMC)
    
    hnum = GetAbsYield(hnum_data, hnum_mc, f"num_FR_for_{suffix}")
    hden = GetAbsYield(hden_data, hden_mc, f"den_FR_for_{suffix}")
    
    if avergeEta:
        hnum = AverageEta2D(hnum)
        hden = AverageEta2D(hden)
        
    assert projX + projY < 2, "only one projection can be done at a time"
    hFR = hnum.Clone(f"{hnum.GetName()}_FR")
    
    if projX:
        # 1D histo along x axis
        hnumN = hnum.ProjectionX(f"{hnum.GetName()}_ProjX")
        hdenN = hden.ProjectionX(f"{hden.GetName()}_ProjX")
        PositiveProtection(hnumN)
        PositiveProtection(hdenN)
        hnumN.Divide(hdenN)
        for iy in range(1, hFR.GetNbinsY()+1):
            for ix in range(1, hFR.GetNbinsX()+1):
                hFR.SetBinContent(ix, iy, hnumN.GetBinContent(ix))
                hFR.SetBinError(ix, iy, hnumN.GetBinError(ix))
    elif projY:
        # 1D histo along y axis
        hnumN = hnum.ProjectionY(f"{hnum.GetName()}_ProjY")
        hdenN = hden.ProjectionY(f"{hden.GetName()}_ProjY")
        PositiveProtection(hnumN)
        PositiveProtection(hdenN)
        hnumN.Divide(hdenN)
        for ix in range(1, hFR.GetNbinsX()+1):
            for iy in range(1, hFR.GetNbinsY()+1):
                hFR.SetBinContent(ix, iy, hnumN.GetBinContent(iy))
                hFR.SetBinError(ix, iy, hnumN.GetBinError(iy))
    else:
        PositiveProtection(hnum)
        PositiveProtection(hden)
        hFR.Divide(hden)
    return hFR

def AverageEta2D(h2):
    """
    average eta+ and eta-, assuming eta is on the y axis
    """
    nbinsy = h2.GetNbinsY()
    for ix in range(1, h2.GetNbinsX()+1):
        for iy in range(1, int(nbinsy/2)+1):
            val = (h2.GetBinContent(ix, iy) + h2.GetBinContent(ix, nbinsy-iy+1)) * 0.5
            err = np.sqrt((h2.GetBinError(ix, iy)**2 + h2.GetBinError(ix, nbinsy-iy+1)**2)*0.5)
            h2.SetBinContent(ix, iy, val)
            h2.SetBinError(ix, iy, err)
            h2.SetBinContent(ix, nbinsy-iy+1, val)
            h2.SetBinError(ix, nbinsy-iy+1, err)
    return h2

def StatUnc2SysUnc(h1, prefix = "stat"):
    hs = []
    for ix in range(1, h1.GetNbinsX()+1):
        hup = h1.Clone(f"{h1.GetName()}_{prefix}_bin{ix}Up")
        hdn = h1.Clone(f"{h1.GetName()}_{prefix}_bin{ix}Down")
        val = h1.GetBinContent(ix)
        err = h1.GetBinError(ix)
        hup.SetBinContent(ix, val+err)
        hdn.SetBinContent(ix, val-err)
        hs.append(hup)
        hs.append(hdn)
    return hs

def GetAbsYield(hdata, hmc, suffix=""):
    habs = hdata.Clone(f"{hdata.GetName()}_{suffix}")
    habs.Add(hmc, -1)
    return habs


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--doPtVsEta", action="store_true", help="do pT vs eta")
    parser.add_argument("--doPtVsMet", action="store_true", help="do pT vs met")
    parser.add_argument("--is5TeV", action="store_true", help="is 5TeV")
    parser.add_argument("--doElectron", action="store_true", help="do electron")
    parser.add_argument("--useChgIso", action="store_true", help="use charged isolation")
    parser.add_argument("--useTwoPhis", action="store_true", help="separate phis > pi / 2 and phis < pi / 2")
    parser.add_argument("--useQCDMC", action="store_true", help="use QCD MC")
    parser.add_argument("--useGenMET", action="store_true", help="use gen MET (Works only on 2017 QCD MC)")
    
    args = parser.parse_args()
    
    doPtVsEta = args.doPtVsEta
    doPtVsMet = args.doPtVsMet
    is5TeV = args.is5TeV
    doElectron = args.doElectron
    useChgIso = args.useChgIso
    useTwoPhis = args.useTwoPhis
    useQCDMC = args.useQCDMC
    useGenMET = args.useGenMET

    doPtVsDeltaPhi = not doPtVsEta and not doPtVsMet
    lepname = "mu" if not doElectron else "e"
    sqrtS = "13TeV" if not is5TeV else "5TeV"
    
    if useQCDMC:
        if useGenMET:
            sqrtS = "QCDMC_GenMET"
        else:
            sqrtS = "QCDMC_PFMET"
    
    isovar = "relIso" if not useChgIso else "pfChIso"

    if doPtVsEta:
        iname = f"root/output_qcdLepPtVsEtaMean_{isovar}_{lepname}nu_{sqrtS}.root"
        var1 = f"{lepname}_pt"
        var2 = "eta"
        x1label = "p_{T}^{l} [GeV]"
        x1min = 25
        x1max = 60
        x2label = "#eta^{l}"
        x2min = -2.4
        x2max = 2.4  
    elif doPtVsMet:
        iname = f"root/output_qcdLepPtVsMetMean_{isovar}_{lepname}nu_{sqrtS}.root"
        var1 = f"{lepname}_pt"
        var2 = "met"
        x1label = "p_{T}^{l} [GeV]"
        x1min = 25
        x1max = 60
        x2label = "p_{T}^{miss} [GeV]"
        x2min = 0
        x2max = 60
    else:
        iname = f"root/output_qcdLepPtVsDeltaPhiMean_{isovar}_{lepname}nu_{sqrtS}.root"
        var1 = f"{lepname}_pt"
        var2 = "deltaPhi"
        x1label = "p_{T}^{l} [GeV]"
        x1min = 25
        x1max = 60
        x2label = "#Delta#phi"
        x2min = 0.
        x2max = 3.2

    var2D = f"{var1}_vs_{var2}"
    outdir = f"plots/{sqrtS}/FRAndClosure/{isovar}/{var2D}"
    if useTwoPhis:
        outdir += "_TwoPhis"
    
    print("input root file name: ", iname)
    f = ROOT.TFile.Open(iname) 

    isogroups = OrderedDict()
    if not doElectron:
        #isogroups["SR"] = ["iso0", "iso1", "iso2", "iso3"]
        isogroups["SR"] = ["iso0", "iso1"]
        isogroups["CR0"] = ["iso2", "iso3"]
        isogroups["CR1"] = ["iso4", "iso5", "iso6"]
        isogroups["CRAll"] = isogroups["CR0"] + isogroups["CR1"]
        #isogroups["CR0"] = ["iso4", "iso5", "iso6", "iso7", "iso8"]
        #isogroups["CR1"] = ["iso9", "iso10", "iso11", "iso12"]
        #isogroups["CR2"] = ["iso13", "iso14", "iso15", "iso16", "iso17", "iso18"]
        #isogroups["CRAll"] = isogroups["CR0"] + isogroups["CR1"] + isogroups["CR2"]
    else:
        isogroups["SR"] = ["iso3", "iso4"]
        isogroups["CR0"] = ["iso5", "iso6", "iso7"]
        isogroups["CR1"] = ["iso8", "iso9"]
        isogroups["CRAll"] = isogroups["CR0"] + isogroups["CR1"]

    etabins = ["lepEta_bin0"]

    wptbins = ["WpT_bin0"]

    chgbins = [lepname+"plus", lepname + "minus"]
    
    if not useTwoPhis:
        phibins = ["_"]
    else:
        phibins = ["_deltaPhiM_", "_deltaPhiP_"]

    # mass bins and groups
    # CR1 is the one used to derive FRs
    # CR2 is for the closure test
    # SR is the signal region
    mtbins = mass_bins_forqcd
    bmaxs = OrderedDict()
    if not doElectron:
        bmaxs["MTCR1"] = 20
    else:
        bmaxs["MTCR1"] = 20
    bmaxs["MTCR2"] = 40
    bmaxs["MTSR"] = mtbins[-1]
    
    mtnames = []
    for imt in range(len(mtbins)-1):
        bname = f"mt{imt}"
        mtnames.append(bname)
    
    if not doElectron:
        frbins = ["mt0", "mt1"]
    else:
        frbins = ["mt0", "mt1"]
    
    mass_bins_groups = OrderedDict()
    mass_bins = OrderedDict()
    mass_bins_groups["MTCR1"] = []
    mass_bins["MTCR1"] = []
    mass_bins_groups["MTCR2"] = []
    mass_bins["MTCR2"] = []
    mass_bins_groups["MTSR"] = []
    mass_bins["MTSR"] = []
    for imt in range(len(mtbins)-1):
        bname = f"mt{imt}"
        bmin = mtbins[imt]
        if bmin < bmaxs["MTCR1"]:
            mass_bins_groups["MTCR1"].append(bname)
            mass_bins["MTCR1"].append(bmin)
        if bmin < bmaxs["MTCR2"]:
            mass_bins_groups["MTCR2"].append(bname)
            mass_bins["MTCR2"].append(bmin)
        if bmin >= bmaxs["MTCR2"]:
            mass_bins_groups["MTSR"].append(bname)
            mass_bins["MTSR"].append(bmin)
            
    for k, v in mass_bins.items():
        v.append( bmaxs[k] )
        mass_bins[k] = np.array(v)
    print("mass_bins", mass_bins)
    print("mass_bins_groups", mass_bins_groups)
    
    mc_suffix = "MC"
    data_suffix = "Data"

    # get the histograms
    histos_count_Data = OrderedDict()
    histos_count_MC = OrderedDict()
    histos_count_QCD = OrderedDict()
    histos_FR = OrderedDict()
    histos_ToSave = []
   
    # to save the central values of predictions
    # for systematic variations 
    for wpt in wptbins:
        histos_count_Data[wpt] = OrderedDict()
        histos_count_MC[wpt] = OrderedDict()
        histos_count_QCD[wpt] = OrderedDict()
        histos_FR[wpt] = OrderedDict()
        for lepeta in etabins:
            histos_count_Data[wpt][lepeta] = OrderedDict()
            histos_count_MC[wpt][lepeta] = OrderedDict()
            histos_count_QCD[wpt][lepeta] = OrderedDict()
            histos_FR[wpt][lepeta] = OrderedDict()
            for chg in chgbins:
                histos_count_Data[wpt][lepeta][chg] = OrderedDict()
                histos_count_MC[wpt][lepeta][chg] = OrderedDict()
                histos_count_QCD[wpt][lepeta][chg] = OrderedDict()
                histos_FR[wpt][lepeta][chg] = OrderedDict()
                for mtn in mtnames:
                    histos_count_Data[wpt][lepeta][chg][mtn] = OrderedDict()
                    histos_count_MC[wpt][lepeta][chg][mtn] = OrderedDict()
                    histos_count_QCD[wpt][lepeta][chg][mtn] = OrderedDict()
                    
                    for phin in phibins:
                        
                        histos_count_Data[wpt][lepeta][chg][mtn][phin] = OrderedDict()
                        histos_count_MC[wpt][lepeta][chg][mtn][phin] = OrderedDict()
                        histos_count_QCD[wpt][lepeta][chg][mtn][phin] = OrderedDict()

                        # read and combine all histograms in different iso groups
                        for isog in isogroups:
                            hdata = None
                            hmc = None
                            for iso in isogroups[isog]:
                                strname = "weight_{}_{}_{}_{}_{}{}".format(chg, iso, wpt, lepeta, mtn, phin)
                                hname = f"histo_wjets_{var2D}_{strname}"

                                if not hdata:
                                    print("Getting " + hname + f"{data_suffix}")
                                    hdata = f.Get(hname + f"{data_suffix}").Clone(hname+"PureData")
                                else:
                                    hdata.Add(f.Get(hname + f"{data_suffix}"))
                                if not hmc:
                                    hmc = f.Get(hname + f"{mc_suffix}").Clone(hname+"PureMC")
                                else:
                                    hmc.Add(f.Get(hname + f"{mc_suffix}"))
                                    
                            # include overflow and underflow bins 
                            hdata = IncludeOverflow2D(hdata, True)
                            hmc = IncludeOverflow2D(hmc, True)
                            
                            if doPtVsEta: 
                                #hdata.RebinX(3)
                                #hmc.RebinX(3)
                                hdata.RebinY(3)
                                hmc.RebinY(3)
                                
                            histos_count_Data[wpt][lepeta][chg][mtn][phin][isog] = hdata
                            histos_count_MC[wpt][lepeta][chg][mtn][phin][isog] = hmc       
                            
                # draw the data-MC comparisons in different mt and iso groups
                # in order to validate signal contaminations
                # plot isolation CR only
                # as the SR will be shown later in the Closure test
                for mtn in mtnames:
                    for phin in phibins:
                        for isog in isogroups:
                            if isog == "SR":
                                continue
                            hdata = histos_count_Data[wpt][lepeta][chg][mtn][phin][isog]
                            hmc = histos_count_MC[wpt][lepeta][chg][mtn][phin][isog]  

                            hdata_var2 = hdata.ProjectionY(f"{hdata.GetName()}_Proj_{var2}_in{isog}")
                            hmc_var2 = hmc.ProjectionY(f"{hmc.GetName()}_Proj_{var2}_in{isog}")
                            DrawDataMCStack(hdata_var2, hmc_var2, None, x2min, x2max, x2label, f"{outdir}/ISO_{isog}/DataMC/histos_wjets_{var2}_DataMCComp_{chg}_{wpt}_{lepeta}_{mtn}{phin}stacked", is5TeV = is5TeV)

                            hdata_var1 = hdata.ProjectionX(f"{hdata.GetName()}_Proj_{var1}_for{isog}")
                            hmc_var1 = hmc.ProjectionX(f"{hmc.GetName()}_Proj_{var1}_for{isog}")
                            DrawDataMCStack(hdata_var1, hmc_var1, None, x1min, x1max, x1label, f"{outdir}/ISO_{isog}/DataMC/histos_wjets_{var1}_DataMCComp_{chg}_{wpt}_{lepeta}_{mtn}{phin}stacked", is5TeV = is5TeV)

                # calculate fake/transfer factor in different mt and iso groups, 
                # as a function of var1 and var2, or var1 (projX) or var2 (projY)
                for phin in phibins:
                    histos_FR[wpt][lepeta][chg][phin] = OrderedDict()
                    for isog in isogroups: 
                        histos_FR[wpt][lepeta][chg][phin][isog] = OrderedDict()
                        hnums_data = []
                        hnums_mc = []
                        hdens_data = []
                        hdens_mc = []
                        # accumulate all histograms in different mt groups
                        for mtn in frbins:
                            hnum_data = histos_count_Data[wpt][lepeta][chg][mtn][phin]["SR"]
                            hnum_mc = histos_count_MC[wpt][lepeta][chg][mtn][phin]["SR"]
                            hden_data = histos_count_Data[wpt][lepeta][chg][mtn][phin][isog]
                            hden_mc = histos_count_MC[wpt][lepeta][chg][mtn][phin][isog]
                            
                            hnums_data.append(hnum_data)
                            hnums_mc.append(hnum_mc)
                            hdens_data.append(hden_data)
                            hdens_mc.append(hden_mc)

                        doProjX = False
                        doProjY = False
                        doAverageEta = (True and doPtVsEta)
                        suffix = phin
                        suffix += "ProjY" if doProjY else "ProjX" if doProjX else "2D" 
                        
                        hFR = GetFR(hnums_data, hnums_mc, hdens_data, hdens_mc, projX=doProjX, projY=doProjY, avergeEta=doAverageEta, suffix = f"{wpt}_{lepeta}_{chg}_{iso}_{suffix}")
                        histos_FR[wpt][lepeta][chg][phin][isog]["central"] = hFR
                       
                        # scale the MC contribution in the anti-isolation region
                        # down by 20% 
                        hFR_crUp = GetFR(hnums_data, hnums_mc, hdens_data, hdens_mc, projX=doProjX, projY=doProjY, avergeEta=doAverageEta, suffix = f"{wpt}_{lepeta}_{chg}_{iso}_{suffix}", scaleDenMC=0.50) 
                        histos_FR[wpt][lepeta][chg][phin][isog]["CRUp"] = hFR_crUp
                       
                        # scale the MC contribution in the signal region 
                        # up by 10% 
                        hFR_srUp = GetFR(hnums_data, hnums_mc, hdens_data, hdens_mc, projX=doProjX, projY=doProjY, avergeEta=doAverageEta, suffix = f"{wpt}_{lepeta}_{chg}_{iso}_{suffix}", scaleNumMC=1.10)
                        histos_FR[wpt][lepeta][chg][phin][isog]["SRUp"] = hFR_srUp

                        DrawHistos([hFR], [], x1min, x1max, x1label, x2min, x2max, x2label, f"{outdir}/ISO_{isog}/FR/histo_wjets_{var2D}_FR_{chg}_{wpt}_{lepeta}_{isog}{suffix}", False, False, False, dologz=False, doth2=True, drawoptions="COLZ,texte", is5TeV=is5TeV) 
                        
                        ymax = 0.35 if not doElectron else 6.0

                        # Draw FRs in 1D
                        hFRxs, labelsy = TH2ToTH1s(hFR, False, x2label)
                        DrawHistos(hFRxs, labelsy, x1min, x1max, x1label, 0, ymax, "Fake Factor", f"{outdir}/ISO_{isog}/FR/histo_wjets_{var1}_FR_{chg}_{wpt}_{lepeta}_{isog}{suffix}_ProjX", False, False, False, is5TeV=is5TeV)

                        hFRys, labelsx = TH2ToTH1s(hFR, True, x1label)
                        DrawHistos(hFRys, labelsx, x2min, x2max, x2label, 0, ymax, "Fake Factor", f"{outdir}/ISO_{isog}/FR/histo_wjets_{var2}_FR_{chg}_{wpt}_{lepeta}_{isog}{suffix}_ProjY", False, False, False, is5TeV=is5TeV)
                        
                # get the predictions in the SR, 
                # using the FR/TF in that anti-isolation region
                hpreds_mt_sr = OrderedDict()
                for frsys in ["central", "CRUp", "SRUp"]:
                    hpreds_mt_sr[frsys] = OrderedDict()
                    for isog in isogroups: 
                        hdatas_var1 = OrderedDict()
                        hmcs_var1 = OrderedDict()
                        hpreds_var1 = OrderedDict()
                        hdatas_var2 = OrderedDict()
                        hmcs_var2 = OrderedDict()
                        hpreds_var2 = OrderedDict()
                        if isog == "SR":
                            continue
                        
                        for phin in phibins:
                            hfr = histos_FR[wpt][lepeta][chg][phin][isog][frsys]
                            for mtn in mtnames:
                                # get the abs yield in the anti-isolation region
                                # then yield * FR = QCD prediction in the signal region
                                hdata_cr = histos_count_Data[wpt][lepeta][chg][mtn][phin][isog]
                                hmc_cr = histos_count_MC[wpt][lepeta][chg][mtn][phin][isog]
                                hcounts_cr = GetAbsYield(hdata_cr, hmc_cr, f"CR_for_{isog}_with_FR{frsys}")

                                h = hcounts_cr.Clone(hcounts_cr.GetName() + "_CT")
                                PositiveProtection(h)

                                hpred = MultiplyH2(h, hfr)
                                histos_count_QCD[wpt][lepeta][chg][mtn][phin]['SR'] = hpred
                               
                                # get the abs yield in the signal region 
                                hdata_sr = histos_count_Data[wpt][lepeta][chg][mtn][phin]["SR"]
                                hmc_sr = histos_count_MC[wpt][lepeta][chg][mtn][phin]["SR"]
                                hcounts_sr = GetAbsYield(hdata_sr, hmc_sr, "SR")

                                for ibinx in range(1, h.GetNbinsX()+1):
                                    for ibiny in range(1, h.GetNbinsY()+1):
                                        if hpred.GetBinContent(ibinx, ibiny) < 0:
                                            print("Found negative: ", ibinx, ibiny, hpred.GetBinContent(ibinx, ibiny), h.GetBinContent(ibinx, ibiny), hfr.GetBinContent(ibinx, ibiny))
                                            sys.exit(1)

                                hpred_var1 = hpred.ProjectionX(f"{hpred.GetName()}_Proj_{var1}")
                                hcounts_sr_var1 = hcounts_sr.ProjectionX(f"{hcounts_sr.GetName()}_Proj_{var1}")

                                hpred_var2 = hpred.ProjectionY(f"{hpred.GetName()}_Proj_{var2}")
                                hcounts_sr_var2 = hcounts_sr.ProjectionY(f"{hcounts_sr.GetName()}_Proj_{var2}")

                                hratio = hpred.Clone(hpred.GetName() + "_ratio")
                                hratio.Divide(hcounts_sr)

                                # compare predicted QCD with abs yield in SR (i.e., data-MC)
                                DrawHistos([hcounts_sr_var1, hpred_var1], ["Obs", "Pred"], x1min, x1max, x1label, 0, None, "Events", f"{outdir}/{frsys}/ISO_{isog}/Closure_mT/histos_wjets_{var1}_CT_{chg}_{wpt}_{lepeta}_{mtn}", dology=False, showratio=True, yrmin = 0.8, yrmax = 1.2, is5TeV = is5TeV, mycolors=[ROOT.kBlack, ROOT.kRed])

                                DrawHistos([hcounts_sr_var2, hpred_var2], ["Obs", "Pred"], x2min, x2max, x2label, 0., None, "Events", f"{outdir}/{frsys}/ISO_{isog}/Closure_mT/histos_wjets_{var2}_CT_{chg}_{wpt}_{lepeta}_{mtn}", dology=False, showratio=True, yrmin = 0.8, yrmax = 1.2, is5TeV = is5TeV, mycolors=[ROOT.kBlack, ROOT.kRed])
                                
                                # draw 2D closure ratio
                                DrawHistos([hratio], [], x1min, x1max, x1label, x2min, x2max, x2label, f"{outdir}/{frsys}/ISO_{isog}/Closure_mT_2D/histo_wjets_{var2D}_CT_{chg}_{wpt}_{lepeta}_{mtn}{phin}", False, False, False, dologz=False, doth2=True, drawoptions="COLZ,text", zmin=0.8, zmax=1.2, is5TeV=is5TeV)

                                # make the stacked histograms 
                                # so that the signal contribution can be visualized
                                hdata_var1 = hdata_sr.ProjectionX(f"{hdata_sr.GetName()}_Proj_{var1}_for{isog}")
                                hmc_var1 = hmc_sr.ProjectionX(f"{hmc_sr.GetName()}_Proj_{var1}_for{isog}")
                                DrawDataMCStack(hdata_var1, hmc_var1, hpred_var1, x1min, x1max, x1label, f"{outdir}/{frsys}/ISO_{isog}/Closure_mT_Stacked/histos_wjets_{var1}_CT_{chg}_{wpt}_{lepeta}_{mtn}{phin}_stacked", is5TeV = is5TeV)

                                kname = mtn + phin
                                hdatas_var1[kname] = hdata_var1
                                hmcs_var1[kname] = hmc_var1
                                hpreds_var1[kname] = hpred_var1
                                
                                hdata_var2 = hdata_sr.ProjectionY(f"{hdata_sr.GetName()}_Proj_{var2}_for{isog}")
                                hmc_var2 = hmc_sr.ProjectionY(f"{hmc_sr.GetName()}_Proj_{var2}_for{isog}")
                                DrawDataMCStack(hdata_var2, hmc_var2, hpred_var2, x2min, x2max, x2label, f"{outdir}/{frsys}/ISO_{isog}/Closure_mT_Stacked/histos_wjets_{var2}_CT_{chg}_{wpt}_{lepeta}_{mtn}{phin}_stacked", is5TeV = is5TeV)

                                hdatas_var2[kname] = hdata_var2
                                hmcs_var2[kname] = hmc_var2
                                hpreds_var2[kname] = hpred_var2

                        # check the var1, var2, and mt distributions in different mt groups
                        for mtgn, mtns in mass_bins_groups.items():
                            # check lepton var1 distribution in that mtgroup
                            strname = f"{lepname}_{var1}_{chg}_{wpt}_{lepeta}_{isog}_{mtgn}_{frsys}"
                            hdata_var1 = LHistos2Hist([hdatas_var1[mtn+phin] for mtn in mtns for phin in phibins], f"hdata_{strname}_Data")
                            hmc_var1   = LHistos2Hist([hmcs_var1  [mtn+phin] for mtn in mtns for phin in phibins], f"hmc_{strname}_MC")
                            hpred_var1 = LHistos2Hist([hpreds_var1[mtn+phin] for mtn in mtns for phin in phibins], f"hqcd_{strname}_QCD")
                            DrawDataMCStack(hdata_var1, hmc_var1, hpred_var1, x1min, x1max, x1label, f"{outdir}/{frsys}/ISO_{isog}/MTG_{mtgn}/histos_wjets_{strname}_stacked", is5TeV = is5TeV)

                            # check lepton var2 distribution in that mtgroup
                            strname = f"{lepname}_{var2}_{chg}_{wpt}_{lepeta}_{isog}_{mtgn}_{frsys}"
                            hdata_var2 = LHistos2Hist([hdatas_var2[mtn+phin] for mtn in mtns for phin in phibins], f"hdata_{strname}_Data")
                            hmc_var2   = LHistos2Hist([hmcs_var2  [mtn+phin] for mtn in mtns for phin in phibins], f"hmc_{strname}_MC")
                            hpred_var2 = LHistos2Hist([hpreds_var2[mtn+phin] for mtn in mtns for phin in phibins], f"hqcd_{strname}_QCD")
                            DrawDataMCStack(hdata_var2, hmc_var2, hpred_var2, x2min, x2max, x2label, f"{outdir}/{frsys}/ISO_{isog}/MTG_{mtgn}/histos_wjets_{strname}_stacked", is5TeV = is5TeV)

                            # check the mT distribution in that mtgroup
                            strname = f"{lepname}_mT_{chg}_{wpt}_{lepeta}_{isog}_{mtgn}_{frsys}"
                            hdata_mt = ROOT.TH1D(f"hdata_{strname}_Data", f"hdata_{strname}_Data", len(mass_bins[mtgn])-1, mass_bins[mtgn]) 
                            hdata_mt.Sumw2()
                            hmc_mt = ROOT.TH1D(f"hmc_{strname}_MC", f"hmc_{strname}_MC", len(mass_bins[mtgn])-1, mass_bins[mtgn])
                            hmc_mt.Sumw2()
                            hpred_mt = ROOT.TH1D(f"hqcd_{strname}_QCD", f"hqcd_{strname}_QCD", len(mass_bins[mtgn])-1, mass_bins[mtgn])
                            hpred_mt.Sumw2()
                            for imtn, mtn in enumerate(mtns):
                                val_data, err_data = IntegralAndError2D([histos_count_Data[wpt][lepeta][chg][mtn][phin]["SR"] for phin in phibins])
                                val_mc, err_mc     = IntegralAndError2D([histos_count_MC  [wpt][lepeta][chg][mtn][phin]["SR"] for phin in phibins])
                                val_pred, err_pred = IntegralAndError2D([histos_count_QCD [wpt][lepeta][chg][mtn][phin]["SR"] for phin in phibins])
                                hdata_mt.SetBinContent(imtn + 1, val_data)
                                hdata_mt.SetBinError(imtn + 1, err_data)
                                hmc_mt.SetBinContent(imtn + 1, val_mc)
                                hmc_mt.SetBinError(imtn + 1, err_mc)
                                hpred_mt.SetBinContent(imtn + 1, val_pred)
                                hpred_mt.SetBinError(imtn + 1, err_pred)

                            DrawDataMCStack(hdata_mt, hmc_mt, hpred_mt, 0, 140, "m_{T}", f"{outdir}/{frsys}/ISO_{isog}/MTG_{mtgn}/histos_wjets_{strname}_stacked", is5TeV = is5TeV)

                            if mtgn == "MTSR":
                                # estimated qcd mT distributions in the SR
                                # estimated using different anti-isolation regions
                                # and different fr systematics
                                # save for further processing in order to be used as
                                # combine inputs
                                hpreds_mt_sr[frsys][isog] = hpred_mt
                                
                            #histos_ToSave.append(hdata_mt)
                            #histos_ToSave.append(hmc_mt)
                            #histos_ToSave.append(hpred_mt)
                
                # proccess the estimated qcd mT distributions in the SR
                # (including stat. and sys. variations)
                # to be used as combine inputs                
                for isog in isogroups:
                    if isog == "SR":
                        continue
                    mtgn = "MTSR"
                    strname = f"{lepname}_mT_{chg}_{wpt}_{lepeta}_{isog}_{mtgn}" 
                    hname = f"hqcd_{strname}_QCD"
                    
                    h_mtsr_central = hpreds_mt_sr["central"][isog]
                    h_mtsr_central.SetName(hname)
                    histos_ToSave.append(h_mtsr_central)
                    
                    # save stat uncs    
                    # use prefix to the variations such that they can be uncorrelated
                    prefix = f"{chg}_{lepeta}_{wpt}_{sqrtS}"
                    hmtsys = StatUnc2SysUnc(h_mtsr_central, prefix = prefix)
                    histos_ToSave += hmtsys 
                    
                    # systematic uncs from CRs 
                    h_mtsr_CRUp = hpreds_mt_sr["CRUp"][isog]
                    h_mtsr_CRUp.SetName(f"{hname}_{prefix}_CRMCUp")
                    h_mtsr_CRDown = SymmetrizeHisto(h_mtsr_central, h_mtsr_CRUp, f"{hname}_{prefix}_CRMCDown")
                    histos_ToSave.append(h_mtsr_CRUp)
                    histos_ToSave.append(h_mtsr_CRDown)
                    
                    # systematic uncs from SRs
                    h_mtsr_SRUp = hpreds_mt_sr["SRUp"][isog]
                    h_mtsr_SRUp.SetName(f"{hname}_{prefix}_SRMCUp")    
                    h_mtsr_SRDown = SymmetrizeHisto(h_mtsr_central, h_mtsr_SRUp, f"{hname}_{prefix}_SRMCDown")
                    histos_ToSave.append(h_mtsr_SRUp)
                    histos_ToSave.append(h_mtsr_SRDown)

    ofile = ROOT.TFile.Open(f"{outdir}/histos_qcdFR_{lepname}_{sqrtS}.root", "RECREATE")
    for h in histos_ToSave:
        h.SetDirectory(ofile)
        h.Write()
    ofile.Close()
    
if __name__ == "__main__":
    main()
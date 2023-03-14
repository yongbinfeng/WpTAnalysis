"""
script to measure FR and carry out closure test
"""
import ROOT
import numpy as np
import sys
import argparse
from collections import OrderedDict
from modules.Binnings import mass_bins_test
from CMSPLOTS.myFunction import DrawHistos, MultiplyH2, PositiveProtection, IntegralAndError2D, CombineOneBin2D, GetRatioPanel, LHistos2Hist, IncludeOverflow2D, TH2ToTH1s

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
    
def GetFR(hnums, hdens, projX = False, projY = False):
    hnum = LHistos2Hist(hnums, f"{hnums[0].GetName()}_Clone")
    hden = LHistos2Hist(hdens, f"{hdens[0].GetName()}_Clone")
    
    assert projX + projY < 2, "only one projection can be done at a time"
    hFR = hnum.Clone(f"{hnum.GetName()}_FR")
    
    if projX:
        # 1D histo along x axis
        hnumN = hnum.ProjectionX(f"{hnum.GetName()}_ProjX")
        hdenN = hden.ProjectionX(f"{hden.GetName()}_ProjX")
        hnumN.Divide(hdenN)
        for iy in range(1, hFR.GetNbinsY()+1):
            for ix in range(1, hFR.GetNbinsX()+1):
                hFR.SetBinContent(ix, iy, hnumN.GetBinContent(ix))
                hFR.SetBinError(ix, iy, hnumN.GetBinError(ix))
    elif projY:
        # 1D histo along y axis
        hnumN = hnum.ProjectionY(f"{hnum.GetName()}_ProjY")
        hdenN = hden.ProjectionY(f"{hden.GetName()}_ProjY")
        hnumN.Divide(hdenN)
        for ix in range(1, hFR.GetNbinsX()+1):
            for iy in range(1, hFR.GetNbinsY()+1):
                hFR.SetBinContent(ix, iy, hnumN.GetBinContent(iy))
                hFR.SetBinError(ix, iy, hnumN.GetBinError(iy))
    else:
        hFR.Divide(hden)
    return hFR


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--doPtVsEta", action="store_true", help="do pT vs eta")
    parser.add_argument("--doPtVsMet", action="store_true", help="do pT vs met")
    parser.add_argument("--is5TeV", action="store_true", help="is 5TeV")
    parser.add_argument("--doElectron", action="store_true", help="do electron")
    parser.add_argument("--useChgIso", action="store_true", help="use charged isolation")
    
    args = parser.parse_args()
    
    doPtVsEta = args.doPtVsEta
    doPtVsMet = args.doPtVsMet
    is5TeV = args.is5TeV
    doElectron = args.doElectron
    useChgIso = args.useChgIso

    doPtVsDeltaPhi = not doPtVsEta and not doPtVsMet
    lepname = "mu" if not doElectron else "e"
    sqrtS = "13TeV" if not is5TeV else "5TeV"
    
    isovar = "relIso" if not useChgIso else "pfChIso"

    if doPtVsEta:
        iname = f"root/output_qcdLepPtVsEtaMean_{isovar}_{lepname}nu_{sqrtS}.root"
        var1 = f"{lepname}_pt"
        var2 = "eta"
        x1label = "p_{T}^{l} [GeV]"
        x1min = 25
        x1max = 70
        x2label = "#eta^{l}"
        x2min = -2.4
        x2max = 2.4  
    elif doPtVsMet:
        iname = f"root/output_qcdLepPtVsMetMean_{isovar}_{lepname}nu_{sqrtS}.root"
        var1 = f"{lepname}_pt"
        var2 = "met"
        x1label = "p_{T}^{l} [GeV]"
        x1min = 25
        x1max = 70
        x2label = "p_{T}^{miss} [GeV]"
        x2min = 0
        x2max = 70
    else:
        iname = f"root/output_qcdLepPtVsDeltaPhiMean_{isovar}_{lepname}nu_{sqrtS}.root"
        var1 = f"{lepname}_pt"
        var2 = "deltaPhi"
        x1label = "p_{T}^{l} [GeV]"
        x1min = 25
        x1max = 70
        x2label = "#Delta#phi"
        x2min = 0.
        x2max = ROOT.TMath.Pi()

    var2D = f"{var1}_vs_{var2}"
    outdir = f"plots/{sqrtS}/FRAndClosure/{isovar}/{var2D}/"

    f = ROOT.TFile.Open(iname) 

    isogroups = OrderedDict()
    if not doElectron:
        isogroups["SR"] = ["iso0", "iso1", "iso2", "iso3"]
        isogroups["CR0"] = ["iso4", "iso5", "iso6", "iso7", "iso8"]
        isogroups["CR1"] = ["iso9", "iso10", "iso11", "iso12"]
        isogroups["CR2"] = ["iso13", "iso14", "iso15", "iso16"]
        isogroups["CRAll"] = isogroups["CR0"] + isogroups["CR1"] + isogroups["CR2"]
    else:
        isogroups["SR"] = ["iso3", "iso4"]
        isogroups["CR0"] = ["iso5", "iso6", "iso7", "iso8", "iso9", "iso10"]

    etabins = ["lepEta_bin0"]

    wptbins = ["WpT_bin0"]

    chgbins = [lepname+"plus", lepname + "minus"]

    # mass bins and groups
    # CR1 is the one used to derive FRs
    # CR2 is for the closure test
    # SR is the signal region
    mtbins = mass_bins_test[0]
    bmaxs = OrderedDict()
    if not doElectron:
        bmaxs["MTCR1"] = 10
    else:
        bmaxs["MTCR1"] = 20
    bmaxs["MTCR2"] = 40
    bmaxs["MTSR"] = mtbins[-1]
    
    mtnames = []
    for imt in range(len(mtbins)-1):
        bname = f"mt{imt}"
        mtnames.append(bname)
    
    if doElectron:
        frbins = ["mt0", "mt1"]
    else:
        frbins = ["mt0"]
    
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
    
    mc_suffix = "MC"
    data_suffix = "Data"

    # get the histograms
    histos_count = OrderedDict()
    histos_count_Data = OrderedDict()
    histos_count_MC = OrderedDict()
    histos_count_QCD = OrderedDict()
    histos_FR = OrderedDict()
    histos_ToSave = []

    for wpt in wptbins:
        histos_count[wpt] = OrderedDict()
        histos_count_Data[wpt] = OrderedDict()
        histos_count_MC[wpt] = OrderedDict()
        histos_count_QCD[wpt] = OrderedDict()
        histos_FR[wpt] = OrderedDict()
        for lepeta in etabins:
            histos_count[wpt][lepeta] = OrderedDict()
            histos_count_Data[wpt][lepeta] = OrderedDict()
            histos_count_MC[wpt][lepeta] = OrderedDict()
            histos_count_QCD[wpt][lepeta] = OrderedDict()
            histos_FR[wpt][lepeta] = OrderedDict()
            for chg in chgbins:
                histos_count[wpt][lepeta][chg] = OrderedDict()
                histos_count_Data[wpt][lepeta][chg] = OrderedDict()
                histos_count_MC[wpt][lepeta][chg] = OrderedDict()
                histos_count_QCD[wpt][lepeta][chg] = OrderedDict()
                histos_FR[wpt][lepeta][chg] = OrderedDict()
                for mtn in mtnames:
                    histos_count[wpt][lepeta][chg][mtn] = OrderedDict()
                    histos_count_Data[wpt][lepeta][chg][mtn] = OrderedDict()
                    histos_count_MC[wpt][lepeta][chg][mtn] = OrderedDict()
                    histos_count_QCD[wpt][lepeta][chg][mtn] = OrderedDict()

                    # read and combine all histograms in different iso groups
                    for isog in isogroups:
                        hdata = None
                        hmc = None
                        for iso in isogroups[isog]:
                            strname = "weight_{}_{}_{}_{}_{}".format(chg, iso, wpt, lepeta, mtn)
                            hname = f"histo_wjets_{var2D}_{strname}"

                            if not hdata:
                                print("Getting " + hname + f"_{data_suffix}")
                                hdata = f.Get(hname + f"_{data_suffix}").Clone(hname+"_PureData")
                            else:
                                hdata.Add(f.Get(hname + f"_{data_suffix}"))
                            if not hmc:
                                hmc = f.Get(hname + f"_{mc_suffix}").Clone(hname+"_PureMC")
                            else:
                                hmc.Add(f.Get(hname + f"_{mc_suffix}"))

                        # include overflow and underflow bins 
                        IncludeOverflow2D(hdata, True)
                        IncludeOverflow2D(hmc, True)
                       
                        if 1: 
                            hdata.RebinX(3)
                            hmc.RebinX(3)
                            hdata.RebinY(5)
                            hmc.RebinY(5)

                        # count = count_Data - count_MC
                        h = hdata.Clone(hname + f"_Subtracted")
                        h.Add(hmc, -1)
                        #PositiveProtection(h)

                        if isog not in histos_count[wpt][lepeta][chg][mtn]:
                            histos_count[wpt][lepeta][chg][mtn][isog] = h
                            histos_count_Data[wpt][lepeta][chg][mtn][isog] = hdata
                            histos_count_MC[wpt][lepeta][chg][mtn][isog] = hmc       
                        else:
                            histos_count[wpt][lepeta][chg][mtn][isog].Add(h) 
                            histos_count_Data[wpt][lepeta][chg][mtn][isog].Add(hdata)
                            histos_count_MC[wpt][lepeta][chg][mtn][isog].Add(hmc)
                                
                # draw the data-MC comparisons in different mt and iso groups
                # in order to validate signal contaminations
                # plot isolation CR only
                # as the SR will be shown later in the Closure test
                for mtn in mtnames:
                    for isog in isogroups:
                        if isog == "SR":
                            continue
                        hdata = histos_count_Data[wpt][lepeta][chg][mtn][isog]
                        hmc = histos_count_MC[wpt][lepeta][chg][mtn][isog]  

                        hdata_var2 = hdata.ProjectionY(f"{hdata.GetName()}_Proj_{var2}_in{isog}")
                        hmc_var2 = hmc.ProjectionY(f"{hmc.GetName()}_Proj_{var2}_in{isog}")
                        DrawDataMCStack(hdata_var2, hmc_var2, None, x2min, x2max, x2label, f"{outdir}/ISO_{isog}/DataMC/histos_wjets_{var2}_DataMCComp_{chg}_{wpt}_{lepeta}_{mtn}_stacked", is5TeV = is5TeV)

                        hdata_var1 = hdata.ProjectionX(f"{hdata.GetName()}_Proj_{var1}_for{isog}")
                        hmc_var1 = hmc.ProjectionX(f"{hmc.GetName()}_Proj_{var1}_for{isog}")
                        DrawDataMCStack(hdata_var1, hmc_var1, None, x1min, x1max, x1label, f"{outdir}/ISO_{isog}/DataMC/histos_wjets_{var1}_DataMCComp_{chg}_{wpt}_{lepeta}_{mtn}_stacked", is5TeV = is5TeV)

                # calculate fake/transfer factor in different mt and iso groups, 
                # as a function of var1 and var2, or var1 (projX) or var2 (projY)
                for isog in isogroups: 
                    hnums = []
                    hdens = []
                    for mtn in frbins:
                        # PositiveProtection(histos_count[wpt][lepeta][chg][mtn][isog])
                        hnums.append( histos_count[wpt][lepeta][chg][mtn]["SR"] )
                        hdens.append( histos_count[wpt][lepeta][chg][mtn][isog] )
                        #hnum.Divide(hden)

                    doProjX = False
                    doProjY = False
                    suffix = "ProjY" if doProjY else "ProjX" if doProjX else "2D" 
                    hFR = GetFR(hnums, hdens, projX=doProjX, projY=doProjY)
                    PositiveProtection(hFR)
                    histos_FR[wpt][lepeta][chg][isog] = hFR

                    DrawHistos([hFR], [], x1min, x1max, x1label, x2min, x2max, x2label, f"{outdir}/ISO_{isog}/FR/histo_wjets_{var2D}_FR_{chg}_{wpt}_{lepeta}_{isog}_{suffix}", False, False, False, dologz=False, doth2=True, drawoptions="COLZ,texte", is5TeV=is5TeV) 
                    
                    # Draw FRs in 1D
                    hFRxs, labelsy = TH2ToTH1s(hFR, False, x2label)
                    DrawHistos(hFRxs, labelsy, x1min, x1max, x1label, 0, None, "Fake Factor", f"{outdir}/ISO_{isog}/FR/histo_wjets_{var1}_FR_{chg}_{wpt}_{lepeta}_{isog}_{suffix}_ProjX", False, False, False, is5TeV=is5TeV)
                    
                    hFRys, labelsx = TH2ToTH1s(hFR, True, x1label)
                    DrawHistos(hFRys, labelsx, x2min, x2max, x2label, 0, None, "Fake Factor", f"{outdir}/ISO_{isog}/FR/histo_wjets_{var2}_FR_{chg}_{wpt}_{lepeta}_{isog}_{suffix}_ProjY", False, False, False, is5TeV=is5TeV)
                        
                # get the predictions in the SR, 
                # using the FR/TF in that anti-isolation region
                for isog in isogroups: 
                    hdatas_var1 = OrderedDict()
                    hmcs_var1 = OrderedDict()
                    hpreds_var1 = OrderedDict()
                    hdatas_var2 = OrderedDict()
                    hmcs_var2 = OrderedDict()
                    hpreds_var2 = OrderedDict()
                    if isog == "SR":
                        continue
                    
                    hfr = histos_FR[wpt][lepeta][chg][isog]
                    for mtn in mtnames:
                        hcounts_cr = histos_count[wpt][lepeta][chg][mtn][isog]
                        hcounts_sr = histos_count[wpt][lepeta][chg][mtn]["SR"]

                        h = hcounts_cr.Clone(hcounts_cr.GetName() + "_CT")
                        PositiveProtection(h)

                        hpred = MultiplyH2(h, hfr)
                        histos_count_QCD[wpt][lepeta][chg][mtn]['SR'] = hpred
                        
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

                        # compare predicted QCD with data-MC
                        DrawHistos([hcounts_sr_var1, hpred_var1], ["Obs", "Pred"], x1min, x1max, x1label, 0, None, "Events", f"{outdir}/ISO_{isog}/Closure_mT/histos_wjets_{var1}_CT_{chg}_{wpt}_{lepeta}_{mtn}", dology=False, showratio=True, yrmin = 0.8, yrmax = 1.2, is5TeV = is5TeV, mycolors=[ROOT.kBlack, ROOT.kRed])

                        DrawHistos([hcounts_sr_var2, hpred_var2], ["Obs", "Pred"], x2min, x2max, x2label, 0., None, "Events", f"{outdir}/ISO_{isog}/Closure_mT/histos_wjets_{var2}_CT_{chg}_{wpt}_{lepeta}_{mtn}", dology=False, showratio=True, yrmin = 0.8, yrmax = 1.2, is5TeV = is5TeV, mycolors=[ROOT.kBlack, ROOT.kRed])

                        # make the stacked histograms 
                        # so that the signal contribution can be visualized
                        hdata = histos_count_Data[wpt][lepeta][chg][mtn]["SR"]
                        hmc = histos_count_MC[wpt][lepeta][chg][mtn]["SR"]  

                        hdata_var2 = hdata.ProjectionY(f"{hdata.GetName()}_Proj_{var2}_for{isog}")
                        hmc_var2 = hmc.ProjectionY(f"{hmc.GetName()}_Proj_{var2}_for{isog}")
                        DrawDataMCStack(hdata_var2, hmc_var2, hpred_var2, x2min, x2max, x2label, f"{outdir}/ISO_{isog}/Closure_mT_Stacked/histos_wjets_{var2}_CT_{chg}_{wpt}_{lepeta}_{mtn}_stacked", is5TeV = is5TeV)
                        
                        hdatas_var2[mtn] = hdata_var2
                        hmcs_var2[mtn] = hmc_var2
                        hpreds_var2[mtn] = hpred_var2
                        
                        hdata_var1 = hdata.ProjectionX(f"{hdata.GetName()}_Proj_{var1}_for{isog}")
                        hmc_var1 = hmc.ProjectionX(f"{hmc.GetName()}_Proj_{var1}_for{isog}")
                        DrawDataMCStack(hdata_var1, hmc_var1, hpred_var1, x1min, x1max, x1label, f"{outdir}/ISO_{isog}/Closure_mT_Stacked/histos_wjets_{var1}_CT_{chg}_{wpt}_{lepeta}_{mtn}_stacked", is5TeV = is5TeV)
                        
                        hdatas_var1[mtn] = hdata_var1
                        hmcs_var1[mtn] = hmc_var1
                        hpreds_var1[mtn] = hpred_var1

                        # draw 2D closure ratio
                        DrawHistos([hratio], [], x1min, x1max, x1label, x2min, x2max, x2label, f"{outdir}/ISO_{isog}/Closure_mT_2D/histo_wjets_{var2D}_CT_{chg}_{wpt}_{lepeta}_{mtn}", False, False, False, dologz=False, doth2=True, drawoptions="COLZ,text", zmin=0.8, zmax=1.2, is5TeV=is5TeV)

                    # check the var1, var2, and mt distributions in different mt groups
                    for mtgn, mtns in mass_bins_groups.items():
                        # check lepton var1 distribution in that mtgroup
                        strname = f"{lepname}_{var1}_{chg}_{wpt}_{lepeta}_{isog}_{mtgn}"
                        hdata_var1 = LHistos2Hist([hdatas_var1[mtn] for mtn in mtns], f"hdata_{strname}_Data")
                        hmc_var1   = LHistos2Hist([hmcs_var1  [mtn] for mtn in mtns], f"hmc_{strname}_MC")
                        hpred_var1 = LHistos2Hist([hpreds_var1[mtn] for mtn in mtns], f"hqcd_{strname}_QCD")
                        DrawDataMCStack(hdata_var1, hmc_var1, hpred_var1, x1min, x1max, x1label, f"{outdir}/ISO_{isog}/MTG_{mtgn}/histos_wjets_{strname}_stacked", is5TeV = is5TeV)

                        # check lepton var2 distribution in that mtgroup
                        strname = f"{lepname}_{var2}_{chg}_{wpt}_{lepeta}_{isog}_{mtgn}"
                        hdata_var2 = LHistos2Hist([hdatas_var2[mtn] for mtn in mtns], f"hdata_{strname}_Data")
                        hmc_var2   = LHistos2Hist([hmcs_var2  [mtn] for mtn in mtns], f"hmc_{strname}_MC")
                        hpred_var2 = LHistos2Hist([hpreds_var2[mtn] for mtn in mtns], f"hqcd_{strname}_QCD")
                        DrawDataMCStack(hdata_var2, hmc_var2, hpred_var2, x2min, x2max, x2label, f"{outdir}/ISO_{isog}/MTG_{mtgn}/histos_wjets_{strname}_stacked", is5TeV = is5TeV)

                        # check the mT distribution in that mtgroup
                        strname = f"{lepname}_mT_{chg}_{wpt}_{lepeta}_{isog}_{mtgn}"
                        hdata_mt = ROOT.TH1D(f"hdata_{strname}_Data", f"hdata_{strname}_Data", len(mass_bins[mtgn])-1, mass_bins[mtgn]) 
                        hdata_mt.Sumw2()
                        hmc_mt = ROOT.TH1D(f"hmc_{strname}_MC", f"hmc_{strname}_MC", len(mass_bins[mtgn])-1, mass_bins[mtgn])
                        hmc_mt.Sumw2()
                        hpred_mt = ROOT.TH1D(f"hqcd_{strname}_QCD", f"hqcd_{strname}_QCD", len(mass_bins[mtgn])-1, mass_bins[mtgn])
                        hpred_mt.Sumw2()
                        for imtn, mtn in enumerate(mtns):
                            val_data, err_data = IntegralAndError2D(histos_count_Data[wpt][lepeta][chg][mtn]["SR"])
                            val_mc, err_mc = IntegralAndError2D(histos_count_MC[wpt][lepeta][chg][mtn]["SR"])
                            val_pred, err_pred = IntegralAndError2D(histos_count_QCD[wpt][lepeta][chg][mtn]["SR"])
                            hdata_mt.SetBinContent(imtn + 1, val_data)
                            hdata_mt.SetBinError(imtn + 1, err_data)
                            hmc_mt.SetBinContent(imtn + 1, val_mc)
                            hmc_mt.SetBinError(imtn + 1, err_mc)
                            hpred_mt.SetBinContent(imtn + 1, val_pred)
                            hpred_mt.SetBinError(imtn + 1, err_pred)

                        DrawDataMCStack(hdata_mt, hmc_mt, hpred_mt, 0, 140, "m_{T}", f"{outdir}/ISO_{isog}/MTG_{mtgn}/histos_wjets_{strname}_stacked", is5TeV = is5TeV)

                        histos_ToSave.append(hdata_mt)
                        histos_ToSave.append(hmc_mt)
                        histos_ToSave.append(hpred_mt)

    ofile = ROOT.TFile.Open(f"{outdir}/histos_qcdFR_{lepname}_{sqrtS}.root", "RECREATE")
    for h in histos_ToSave:
        h.SetDirectory(ofile)
        h.Write()
    ofile.Write()
    ofile.Close()
    
if __name__ == "__main__":
    main()
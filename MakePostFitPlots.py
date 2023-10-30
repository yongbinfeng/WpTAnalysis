"""
script to make postfit comparisons
"""
from modules.postFitScripts import MakeDataMCPlot, result2json, MakeWpTPostFitPlots, GetPOIValue, ComparePOIs, DumpGroupImpacts, WriteOutputToText
from modules.CombineHarvester.plotImpacts import plotImpacts
from modules.Binnings import mass_bins_w, mass_bins_z, mass_bins_test
from modules.Utils import FormatTable
import ROOT
import numpy as np
import argparse
from collections import OrderedDict
ROOT.gROOT.SetBatch(True)

parser = argparse.ArgumentParser(description="Make postfit plots and tables")
parser.add_argument("--doInclusive", action="store_true", dest="doInclusive",help="Run on inclusive results; false runs on fiducial results.")
parser.add_argument("--do13TeV", action="store_true", dest="do13TeV",help="Run on 13TeV results.")
parser.add_argument("--do5TeV", action="store_true", dest="do5TeV",help="Run on 5TeV results.")
parser.add_argument("--combineSqrtS", action="store_true", dest="combineSqrtS",help="Combine the results of 5TeV and 13TeV; false runs on each year separately.")
parser.add_argument("--doElectron", action="store_true", dest="doElectron",help="Run on electron channel;")
parser.add_argument("--doMuon", action="store_true", dest="doMuon",help="Run on muon channel;")
parser.add_argument("--reweightZpt", action="store_true", dest="reweightZpt", help="use the Z pt reweighted histograms")

args = parser.parse_args()

doDifferential = False
# boolean flag to config if the pull
# distribution should be included in the plots
showPULL = False

doInclusive = args.doInclusive
do13TeV = args.do13TeV
do5TeV = args.do5TeV
combineSqrtS = args.combineSqrtS

doZptReweight = args.reweightZpt
cmb_suffix = "default" if doZptReweight else "noZptReweight"

doElectron = args.doElectron
doMuon = args.doMuon
doCombineChannel = (doElectron and doMuon)

assert doElectron + doMuon >= 1, "At least one of doElectron, doMuon must be True"
assert do13TeV + do5TeV >= 1, "At least one of do13TeV, do5TeV must be True"

sqrtSs = []
if do13TeV:
    sqrtSs.append("13TeV")
if do5TeV:
    sqrtSs.append("5TeV")

if not doDifferential:
    idir = "Inclusive" if doInclusive else "Fiducial"
    for fits in ["asimov", "data"]:
    #for fits in ["data"]:
        ntests = len(mass_bins_test)

        dvals_lep_pos = {}
        dvals_lep_pos["5TeV"] = []
        dvals_lep_pos["13TeV"] = []
        dvals_lep_neg = {}
        dvals_lep_neg["5TeV"] = []
        dvals_lep_neg["13TeV"] = []
        dvals_leplep = {}
        dvals_leplep["5TeV"] = []
        dvals_leplep["13TeV"] = []
        derrs_lep_pos = {}
        derrs_lep_pos["5TeV"] = []
        derrs_lep_pos["13TeV"] = []
        derrs_lep_neg = {}
        derrs_lep_neg["5TeV"] = []
        derrs_lep_neg["13TeV"] = []
        derrs_leplep = {}
        derrs_leplep["5TeV"] = []
        derrs_leplep["13TeV"] = []

        #sqrtSs = ["13TeV", "5TeV"]
        #sqrtSs = ["5TeV"]

        for idx in range(ntests):
            if idx != 0:
                continue
            mass_bins = mass_bins_test[idx]

            starts_mT = []
            starts_mT.append( mass_bins[0])

            binBase = 1

            for sqrtS in sqrtSs:
                #binBase = 1
                is5TeV = (sqrtS == "5TeV")
                lumi_unc = 0.017 if not is5TeV else 0.019

                workdir = f"forCombine_{cmb_suffix}/{idir}/test{idx}/commands/scripts/"
                suffix = "combined" if doCombineChannel else "mu" if doMuon else "e"
                if combineSqrtS:
                    filename = workdir + f"card_{suffix}"
                else:
                    filename = workdir + f"card_{suffix}_{sqrtS}"
                
                if fits == "asimov":
                    filename += "_asimov"
                filename += ".root"

                outdir = f"forCombine_{cmb_suffix}/{idir}/test{idx}/results_{suffix}"
                if fits == "asimov":
                    outdir += "_asimov"

                val_lep_pos, err_lep_pos = GetPOIValue(filename, f"lepplus_{sqrtS}_sig_mu")
                val_lep_neg, err_lep_neg = GetPOIValue(filename, f"lepminus_{sqrtS}_sig_mu")
                val_leplep,  err_leplep  = GetPOIValue(filename, f"leplep_{sqrtS}_sig_mu")
                dvals_lep_pos[sqrtS].append( val_lep_pos )
                derrs_lep_pos[sqrtS].append( err_lep_pos )
                dvals_lep_neg[sqrtS].append( val_lep_neg )
                derrs_lep_neg[sqrtS].append( err_lep_neg )
                dvals_leplep[sqrtS] .append( val_leplep  )
                derrs_leplep[sqrtS] .append( err_leplep  )

                # todo: binBase needs to be recomputed if not combine all years
                #binBase =  len(mass_bins)*4 + len(mass_bins_z)*2 - 5 if is5TeV and combineSqrtS else 1

                nevts_postfit = OrderedDict()
                nevts_withCut_postfit = OrderedDict()
                nevts_prefit = OrderedDict()
                nevts_withCut_prefit = OrderedDict()
                if doMuon:
                    # prefit
                    nevts_prefit['muplus'], nevts_withCut_prefit['muplus']   = MakeDataMCPlot(filename, "muplus",  mass_bins, f"_mT{idx}_{sqrtS}", False, startbin = binBase, is5TeV=is5TeV, outdir = f"{outdir}/prefits", doPostfit=False)
                    nevts_prefit['muminus'], nevts_withCut_prefit['muminus'] = MakeDataMCPlot(filename, "muminus", mass_bins, f"_mT{idx}_{sqrtS}", False, startbin = binBase + len(mass_bins) - 1, is5TeV=is5TeV, outdir = f"{outdir}/prefits", doPostfit=False)
                    nevts_prefit['mumu'], nevts_withCut_prefit['mumu']       = MakeDataMCPlot(filename, "mumu", mass_bins_z,  f"_mT{idx}_{sqrtS}", False, startbin = binBase + len(mass_bins)*2 - 2, is5TeV=is5TeV, outdir = f"{outdir}/prefits", doPostfit=False)
                    # postfit
                    nevts_postfit['muplus'], nevts_withCut_postfit["muplus"]   = MakeDataMCPlot(filename, "muplus",  mass_bins, f"_mT{idx}_{sqrtS}", showPULL, startbin = binBase, is5TeV=is5TeV, outdir = f"{outdir}/postfits")
                    nevts_postfit['muminus'], nevts_withCut_postfit["muminus"] = MakeDataMCPlot(filename, "muminus", mass_bins, f"_mT{idx}_{sqrtS}", showPULL, startbin = binBase + len(mass_bins) - 1, is5TeV=is5TeV, outdir = f"{outdir}/postfits")
                    nevts_postfit['mumu'], nevts_withCut_postfit["mumu"]       = MakeDataMCPlot(filename, "mumu", mass_bins_z,  f"_mT{idx}_{sqrtS}", showPULL, startbin = binBase + len(mass_bins)*2 - 2, is5TeV=is5TeV, outdir = f"{outdir}/postfits", yrmin = 0.86, yrmax = 1.14)

                    # postfit log
                    MakeDataMCPlot(filename, "muplus",  mass_bins, f"_mT{idx}_{sqrtS}", showPULL, startbin = binBase, is5TeV=is5TeV, outdir = f"{outdir}/postfits_log", dology=True)
                    MakeDataMCPlot(filename, "muminus", mass_bins, f"_mT{idx}_{sqrtS}", showPULL, startbin = binBase + len(mass_bins) - 1, is5TeV=is5TeV, outdir = f"{outdir}/postfits_log", dology=True)
                    MakeDataMCPlot(filename, "mumu", mass_bins_z,  f"_mT{idx}_{sqrtS}", showPULL, startbin = binBase + len(mass_bins)*2 - 2, is5TeV=is5TeV, outdir = f"{outdir}/postfits_log", dology=True, yrmin = 0.86, yrmax = 1.14)

                    binBase = binBase + len(mass_bins)*2 + len(mass_bins_z) - 3

                if doElectron:
                    # prefit
                    nevts_prefit['eplus'], nevts_withCut_prefit['eplus']   = MakeDataMCPlot(filename, "eplus",  mass_bins,  f"_mT{idx}_{sqrtS}", False, startbin = binBase, is5TeV=is5TeV, outdir = f"{outdir}/prefits", doPostfit=False)
                    nevts_prefit['eminus'], nevts_withCut_prefit['eminus']  = MakeDataMCPlot(filename, "eminus", mass_bins,  f"_mT{idx}_{sqrtS}", False, startbin = binBase + len(mass_bins) - 1, is5TeV=is5TeV, outdir = f"{outdir}/prefits", doPostfit=False)
                    nevts_prefit['ee'], nevts_withCut_prefit['ee']          = MakeDataMCPlot(filename, "ee",   mass_bins_z,  f"_mT{idx}_{sqrtS}", False, startbin = binBase + len(mass_bins)*2 - 2, is5TeV=is5TeV, outdir = f"{outdir}/prefits", doPostfit=False)
                    # postfit
                    nevts_postfit['eplus'], nevts_withCut_postfit['eplus']   = MakeDataMCPlot(filename, "eplus",  mass_bins,  f"_mT{idx}_{sqrtS}", showPULL, startbin = binBase, is5TeV=is5TeV, outdir = f"{outdir}/postfits")
                    nevts_postfit['eminus'], nevts_withCut_postfit['eminus'] = MakeDataMCPlot(filename, "eminus", mass_bins,  f"_mT{idx}_{sqrtS}", showPULL, startbin = binBase + len(mass_bins) - 1, is5TeV=is5TeV, outdir = f"{outdir}/postfits")
                    nevts_postfit['ee'], nevts_withCut_postfit['ee']         = MakeDataMCPlot(filename, "ee",   mass_bins_z,  f"_mT{idx}_{sqrtS}", showPULL, startbin = binBase + len(mass_bins)*2 - 2, is5TeV=is5TeV, outdir = f"{outdir}/postfits", yrmin = 0.86, yrmax = 1.14)

                    MakeDataMCPlot(filename, "eplus",  mass_bins,  f"_mT{idx}_{sqrtS}", showPULL, startbin = binBase, is5TeV=is5TeV, outdir = f"{outdir}/postfits_log", dology=True)
                    MakeDataMCPlot(filename, "eminus", mass_bins,  f"_mT{idx}_{sqrtS}", showPULL, startbin = binBase + len(mass_bins) - 1, is5TeV=is5TeV, outdir = f"{outdir}/postfits_log", dology=True)
                    MakeDataMCPlot(filename, "ee",   mass_bins_z,  f"_mT{idx}_{sqrtS}", showPULL, startbin = binBase + len(mass_bins)*2 - 2, is5TeV=is5TeV, outdir = f"{outdir}/postfits_log", yrmin = 0.86, yrmax = 1.14, dology=True)

                    binBase = binBase + len(mass_bins)*2 + len(mass_bins_z) - 3

                result2json(filename, f"lepplus_{sqrtS}_sig_mu",  f"{outdir}/json/impacts_lepplus.json")
                result2json(filename, f"lepminus_{sqrtS}_sig_mu", f"{outdir}/json/impacts_lepminus.json")
                result2json(filename, f"leplep_{sqrtS}_sig_mu",   f"{outdir}/json/impacts_leplep.json")
                result2json(filename, f"Winc_{sqrtS}_sig_sumxsec", f"{outdir}/json/impacts_Winc.json", "nuisance_impact_sumpois")
                result2json(filename, f"WZRatio_{sqrtS}_ratiometaratio",   f"{outdir}/json/impacts_WZRatio.json",   "nuisance_impact_ratiometapois")
                result2json(filename, f"WchgRatio_{sqrtS}_ratiometaratio", f"{outdir}/json/impacts_WchgRatio.json", "nuisance_impact_ratiometapois")
                #result2json(filename, "WchgAsym_chargemetaasym",  f"{outdir}/impacts_WchgAsym.json",  "nuisance_impact_chargemetapois")

                plotImpacts(f"{outdir}/json/impacts_lepplus.json",  f"{outdir}/impacts/impacts_lepplus_mT{idx}_{sqrtS}", poiname="#mu (W^{+})")
                plotImpacts(f"{outdir}/json/impacts_lepminus.json", f"{outdir}/impacts/impacts_lepminus_mT{idx}_{sqrtS}", poiname="#mu (W^{-})")
                plotImpacts(f"{outdir}/json/impacts_leplep.json",   f"{outdir}/impacts/impacts_leplep_mll_mT{idx}_{sqrtS}", poiname="#mu (Z)")
                plotImpacts(f"{outdir}/json/impacts_Winc.json",     f"{outdir}/impacts/impacts_Winc_mT{idx}_{sqrtS}", poiname="#mu (W^{#pm})")
                plotImpacts(f"{outdir}/json/impacts_WZRatio.json",   f"{outdir}/impacts/impacts_WZRatio_mT{idx}_{sqrtS}", poiname="#mu (W^{#pm}/Z)")
                plotImpacts(f"{outdir}/json/impacts_WchgRatio.json", f"{outdir}/impacts/impacts_WchgRatio_mT{idx}_{sqrtS}", poiname="#mu (W^{+}/W^{-})")

                impacts = OrderedDict()
                impacts['lepplus']  = DumpGroupImpacts(filename, f"lepplus_{sqrtS}_sig_mu")
                impacts['lepminus'] = DumpGroupImpacts(filename, f"lepminus_{sqrtS}_sig_mu")
                impacts['leplep']   = DumpGroupImpacts(filename, f"leplep_{sqrtS}_sig_mu")
                impacts['Winc']     = DumpGroupImpacts(filename, f"Winc_{sqrtS}_sig_sumxsec",         "nuisance_group_impact_sumpois")
                impacts['WOverZ']   = DumpGroupImpacts(filename, f"WZRatio_{sqrtS}_ratiometaratio",   "nuisance_group_impact_ratiometapois")
                impacts['WpOverWm'] = DumpGroupImpacts(filename, f"WchgRatio_{sqrtS}_ratiometaratio", "nuisance_group_impact_ratiometapois")
                #impacts['WAsym']    = DumpGroupImpacts(filename, "WchgAsym_chargemetaasym",  "nuisance_group_impact_chargemetapois")

                ## print out the nevts pre and post fit information
                # prefit
                outputs = FormatTable(nevts_prefit, isYield=True)
                print(outputs)
                WriteOutputToText(outputs, f"{outdir}/tables/nevts_prefit_{sqrtS}.tex")
                print("\n\n\n\n")
                outputs = FormatTable(nevts_withCut_prefit, isYield=True)
                print(outputs)
                WriteOutputToText(outputs, f"{outdir}/tables/nevts_withCut_prefit_{sqrtS}.tex")
                print("\n\n\n\n")
                # postfit 
                outputs = FormatTable(nevts_postfit, isYield=True)
                print(outputs)
                WriteOutputToText(outputs, f"{outdir}/tables/nevts_postfit_{sqrtS}.tex")
                print("\n\n\n\n")
                outputs = FormatTable(nevts_withCut_postfit, isYield=True)
                print(outputs)
                WriteOutputToText(outputs, f"{outdir}/tables/nevts_withCut_postfit_{sqrtS}.tex")
                print("\n\n\n\n")
                # impacts
                outputs = FormatTable(impacts, precision=2)
                print(outputs)
                WriteOutputToText(outputs, f"{outdir}/tables/impacts_{sqrtS}.tex")

            if combineSqrtS:
                result2json(filename, "sqrtS_Wplus_ratio_ratiometaratio",   f"{outdir}/json/impacts_Wplus_ratio_sqrtS.json", "nuisance_impact_ratiometapois")
                result2json(filename, "sqrtS_Wminus_ratio_ratiometaratio",  f"{outdir}/json/impacts_Wminus_ratio_sqrtS.json",   "nuisance_impact_ratiometapois")
                result2json(filename, "sqrtS_Winc_ratio_ratiometaratio",    f"{outdir}/json/impacts_Winc_ratio_sqrtS.json", "nuisance_impact_ratiometapois")
                result2json(filename, "sqrtS_Zinc_ratio_ratiometaratio",    f"{outdir}/json/impacts_Zinc_ratio_sqrtS.json",  "nuisance_impact_ratiometapois")
                # double ratios
                result2json(filename, "sqrtS_WchgRatio_ratio_doubleratiometaratio",     f"{outdir}/json/impacts_WchgRatio_sqrtS.json", "nuisance_impact_doubleratiometapois")
                result2json(filename, "sqrtS_WZRatio_ratio_doubleratiometaratio",       f"{outdir}/json/impacts_WZRatio_sqrtS.json", "nuisance_impact_doubleratiometapois")

                plotImpacts(f"{outdir}/json/impacts_Wplus_ratio_sqrtS.json",  f"{outdir}/impacts/impacts_Wplus_ratio_mT{idx}_sqrtS", poiname="#mu (W^{+})")
                plotImpacts(f"{outdir}/json/impacts_Wminus_ratio_sqrtS.json", f"{outdir}/impacts/impacts_Wminus_ratio_mT{idx}_sqrtS", poiname="#mu (W^{-})")
                plotImpacts(f"{outdir}/json/impacts_Winc_ratio_sqrtS.json",   f"{outdir}/impacts/impacts_Winc_ratio_mT{idx}_sqrtS", poiname="#mu (W)")
                plotImpacts(f"{outdir}/json/impacts_Zinc_ratio_sqrtS.json",   f"{outdir}/impacts/impacts_Zinc_ratio_mT{idx}_sqrtS", poiname="#mu (Z)")
                plotImpacts(f"{outdir}/json/impacts_WchgRatio_sqrtS.json",    f"{outdir}/impacts/impacts_WchgRatio_mT{idx}_sqrtS", poiname="#mu (W^{+}/W^{-})")
                plotImpacts(f"{outdir}/json/impacts_WZRatio_sqrtS.json",      f"{outdir}/impacts/impacts_WZRatio_mT{idx}_sqrtS", poiname="#mu (W/Z)")

                impacts = OrderedDict()
                impacts['Wplus']  = DumpGroupImpacts(filename, "sqrtS_Wplus_ratio_ratiometaratio",   "nuisance_group_impact_ratiometapois")
                impacts['Wminus'] = DumpGroupImpacts(filename, "sqrtS_Wminus_ratio_ratiometaratio",  "nuisance_group_impact_ratiometapois")
                impacts['Winc']   = DumpGroupImpacts(filename, "sqrtS_Winc_ratio_ratiometaratio",    "nuisance_group_impact_ratiometapois")
                impacts['Zinc']   = DumpGroupImpacts(filename, "sqrtS_Zinc_ratio_ratiometaratio",    "nuisance_group_impact_ratiometapois")
                impacts['WOverZ']   = DumpGroupImpacts(filename, "sqrtS_WZRatio_ratio_doubleratiometaratio",     "nuisance_group_impact_doubleratiometapois")
                impacts['WpOverWm'] = DumpGroupImpacts(filename, "sqrtS_WchgRatio_ratio_doubleratiometaratio",   "nuisance_group_impact_doubleratiometapois")

                outputs = FormatTable(impacts, precision=2)
                print(outputs)
                WriteOutputToText(outputs, f"{outdir}/tables/impacts_sqrtS_mT{idx}.tex")

    
        for sqrtS in sqrtSs:
            vals_mT = np.array(starts_mT)
            vals_lep_pos = np.array(dvals_lep_pos[sqrtS])
            errs_lep_pos = np.array(derrs_lep_pos[sqrtS])
            vals_lep_neg = np.array(dvals_lep_neg[sqrtS])
            errs_lep_neg = np.array(derrs_lep_neg[sqrtS])
            vals_leplep  = np.array(dvals_leplep[sqrtS])
            errs_leplep  = np.array(derrs_leplep[sqrtS])

            errs_lep_pos = np.sqrt(errs_lep_pos**2 - (lumi_unc * vals_lep_pos)**2)
            errs_lep_neg = np.sqrt(errs_lep_neg**2 - (lumi_unc * vals_lep_neg)**2)
            errs_leplep  = np.sqrt(errs_leplep**2  - (lumi_unc * vals_leplep)**2)

            # scale back to 1
            vals_lep_pos = vals_lep_pos / vals_lep_pos[0]
            vals_lep_neg = vals_lep_neg / vals_lep_neg[0]
            vals_leplep  = vals_leplep  / vals_leplep [0]

            labels = [
                "W^{+} #rightarrow #mu^{+}#nu",
                "W^{-} #rightarrow #mu^{-}#bar{#nu}",
                "W^{+} #rightarrow e^{+}#nu",
                "W^{-} #rightarrow e^{-}#bar{#nu}",
                "W^{+} #rightarrow l^{+}#nu",
                "W^{-} #rightarrow l^{-}#bar{#nu}",
                "Z #rightarrow l^{+}l^{-}",
            ]
            markers = [22, 23, 34]
            colors = [2, 4, 12]
            #ComparePOIs(vals_mT, [vals_mu_pos, vals_mu_neg], [errs_mu_pos, errs_mu_neg], labels[:2], colors[:2], markers[:2], "test_mu")
            #ComparePOIs(vals_mT, [vals_e_pos, vals_e_neg], [errs_e_pos, errs_e_neg], labels[2:], colors[2:], markers[2:], "test_e")
            #ComparePOIs(vals_mT, [vals_mu_pos, vals_mu_neg, vals_e_pos, vals_e_neg], [errs_mu_pos, errs_mu_neg, errs_e_pos, errs_e_neg], labels, colors, markers, "test")
            ComparePOIs(vals_mT, [vals_lep_pos, vals_lep_neg, vals_leplep], [errs_lep_pos, errs_lep_neg, errs_leplep], labels[4:], colors, markers, f"poi_vs_mT_{sqrtS}", is5TeV = is5TeV)

            # plot the impacts on QCD
            #result2json("cards/datacard_muplus_lepEta_bin0_WpT_bin0", "QCD_muplus_lepEta_bin0_WpT_bin0_mu", "cards/impacts_muplus_lepEta_bin0_WpT_bin0_qcd.json")
            #result2json("cards/datacard_eplus",       "QCD_eplus_lepEta_bin1_WpT_bin0_mu",  "cards/impacts_eplus_qcd.json")

            #plotImpacts("cards/impacts_muplus_lepEta_bin0_WpT_bin0_qcd.json", "impacts_muplus_lepEta_bin0_qcd")
            #plotImpacts("cards/impacts_eplus_qcd.json",       "impacts_eplus_lepEta_qcd")


if doDifferential:
    MakeDataMCPlot("cards/datacard_muplus_lepEta_bin0_WpT.root", "muplus", "WpT", showPULL)

    wpttruthbins  = ["WpT_truth_bin1", "WpT_truth_bin2", "WpT_truth_bin3", "WpT_truth_bin4", "WpT_truth_bin5", "WpT_truth_bin6", "WpT_truth_bin7", "WpT_truth_bin8", "WpT_truth_bin9"]
   
    jsonNames = []
    for wpttruthbin in wpttruthbins:
        result2json("cards/datacard_muplus_lepEta_bin0_WpT.root", "w_muplus_{}_sig_mu".format(wpttruthbin), "cards/impacts_muplus_{}_sig_mu.json".format(wpttruthbin))
        jsonNames.append("cards/impacts_muplus_{}_sig_mu.json".format(wpttruthbin))
        plotImpacts("cards/impacts_muplus_{}_sig_mu.json".format(wpttruthbin), "impacts_muplus_{}".format(wpttruthbin))

    MakeWpTPostFitPlots(jsonNames)


"""
script to make postfit comparisons
"""
from modules.postFitScripts import MakePostPlot, result2json, MakeWpTPostFitPlots, GetPOIValue, ComparePOIs, DumpGroupImpacts, WriteOutputToText
from modules.CombineHarvester.plotImpacts import plotImpacts
from modules.Binnings import mass_bins_w, mass_bins_z, mass_bins_test
from modules.Utils import FormatTable
import ROOT
import numpy as np
from collections import OrderedDict
ROOT.gROOT.SetBatch(True)

doInclusive = True
doWpT = False
doCombineYear = True

# boolean flag to config if the pull
# distribution should be included in the plots
showPULL = True
if doInclusive:
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

    sqrtSs = ["5TeV", "13TeV"]

    for idx in range(ntests):
        if idx != 4:
            continue
        mass_bins = mass_bins_test[idx]

        starts_mT = []
        starts_mT.append( mass_bins[0])

        for sqrtS in sqrtSs:
            is5TeV = (sqrtS == "5TeV")
            lumi_unc = 0.017 if not is5TeV else 0.019

            workdir = f"forCombine/test/test{idx}/commands/scripts/"
            if doCombineYear:
                filename = workdir + "card_combined.root"
            else:
                filename = workdir + f"card_combined_{sqrtS}.root"

            outdir = f"forCombine/test/test{idx}/results/"

            val_lep_pos, err_lep_pos = GetPOIValue(filename, f"lepplus_{sqrtS}_sig_mu")
            val_lep_neg, err_lep_neg = GetPOIValue(filename, f"lepminus_{sqrtS}_sig_mu")
            val_leplep,  err_leplep  = GetPOIValue(filename, f"leplep_{sqrtS}_sig_mu")
            dvals_lep_pos[sqrtS].append( val_lep_pos )
            derrs_lep_pos[sqrtS].append( err_lep_pos )
            dvals_lep_neg[sqrtS].append( val_lep_neg )
            derrs_lep_neg[sqrtS].append( err_lep_neg )
            dvals_leplep[sqrtS] .append( val_leplep  )
            derrs_leplep[sqrtS] .append( err_leplep  )

            binBase =  len(mass_bins)*4 + len(mass_bins_z)*2 - 5 if is5TeV and doCombineYear else 1

            nevts = OrderedDict()
            nevts['muplus']  = MakePostPlot(filename, "muplus",  mass_bins, f"_mT{idx}_{sqrtS}", showPULL, startbin = binBase, is5TeV=is5TeV, outdir = f"{outdir}/postfits")
            nevts['muminus'] = MakePostPlot(filename, "muminus", mass_bins, f"_mT{idx}_{sqrtS}", showPULL, startbin = binBase + len(mass_bins) - 1, is5TeV=is5TeV, outdir = f"{outdir}/postfits")
            nevts['mumu']    = MakePostPlot(filename, "mumu", mass_bins_z,  f"_mT{idx}_{sqrtS}", showPULL, startbin = binBase + len(mass_bins)*2 - 2, is5TeV=is5TeV, outdir = f"{outdir}/postfits")
            nevts['eplus']   = MakePostPlot(filename, "eplus",  mass_bins,  f"_mT{idx}_{sqrtS}", showPULL, startbin = binBase + len(mass_bins)*2 + len(mass_bins_z) - 3, is5TeV=is5TeV, outdir = f"{outdir}/postfits")
            nevts['eminus']  = MakePostPlot(filename, "eminus", mass_bins,  f"_mT{idx}_{sqrtS}", showPULL, startbin = binBase + len(mass_bins)*3 + len(mass_bins_z) - 4, is5TeV=is5TeV, outdir = f"{outdir}/postfits")
            nevts['ee']      = MakePostPlot(filename, "ee",   mass_bins_z,  f"_mT{idx}_{sqrtS}", showPULL, startbin = binBase + len(mass_bins)*4 + len(mass_bins_z) - 5, is5TeV=is5TeV, outdir = f"{outdir}/postfits")

            result2json(filename, f"lepplus_{sqrtS}_sig_mu",  f"{outdir}/json/impacts_lepplus.json")
            result2json(filename, f"lepminus_{sqrtS}_sig_mu", f"{outdir}/json/impacts_lepminus.json")
            result2json(filename, f"leplep_{sqrtS}_sig_mu",   f"{outdir}/json/impacts_leplep.json")
            result2json(filename, f"Winc_{sqrtS}_sig_sumxsec", f"{outdir}/json/impacts_Winc.json", "nuisance_impact_sumpois")
            result2json(filename, f"WZRatio_{sqrtS}_ratiometaratio",   f"{outdir}/json/impacts_WZRatio.json",   "nuisance_impact_ratiometapois")
            result2json(filename, f"WchgRatio_{sqrtS}_ratiometaratio", f"{outdir}/json/impacts_WchgRatio.json", "nuisance_impact_ratiometapois")
            #result2json(filename, "WchgAsym_chargemetaasym",  f"{outdir}/impacts_WchgAsym.json",  "nuisance_impact_chargemetapois")

            plotImpacts(f"{outdir}/json/impacts_lepplus.json",  f"{outdir}/impacts/impacts_lepplus_mT{idx}_{sqrtS}")
            plotImpacts(f"{outdir}/json/impacts_lepminus.json", f"{outdir}/impacts/impacts_lepminus_mT{idx}_{sqrtS}")
            plotImpacts(f"{outdir}/json/impacts_leplep.json",   f"{outdir}/impacts/impacts_leplep_mll_mT{idx}_{sqrtS}")
            plotImpacts(f"{outdir}/json/impacts_Winc.json",     f"{outdir}/impacts/impacts_Winc_mT{idx}_{sqrtS}")
            plotImpacts(f"{outdir}/json/impacts_WZRatio.json",   f"{outdir}/impacts/impacts_WZRatio_mT{idx}_{sqrtS}")
            plotImpacts(f"{outdir}/json/impacts_WchgRatio.json", f"{outdir}/impacts/impacts_WchgRatio_mT{idx}_{sqrtS}")

            impacts = OrderedDict()
            impacts['lepplus']  = DumpGroupImpacts(filename, f"lepplus_{sqrtS}_sig_mu")
            impacts['lepminus'] = DumpGroupImpacts(filename, f"lepminus_{sqrtS}_sig_mu")
            impacts['leplep']   = DumpGroupImpacts(filename, f"leplep_{sqrtS}_sig_mu")
            impacts['Winc']     = DumpGroupImpacts(filename, f"Winc_{sqrtS}_sig_sumxsec",         "nuisance_group_impact_sumpois")
            impacts['WOverZ']   = DumpGroupImpacts(filename, f"WZRatio_{sqrtS}_ratiometaratio",   "nuisance_group_impact_ratiometapois")
            impacts['WpOverWm'] = DumpGroupImpacts(filename, f"WchgRatio_{sqrtS}_ratiometaratio", "nuisance_group_impact_ratiometapois")
            #impacts['WAsym']    = DumpGroupImpacts(filename, "WchgAsym_chargemetaasym",  "nuisance_group_impact_chargemetapois")

            ## print out the nevts information
            outputs = FormatTable(nevts, caption="Event yield at 13 TeV", label = f"tab:nevts_{sqrtS}")
            print(outputs)
            WriteOutputToText(outputs, f"{outdir}/tables/nevts_{sqrtS}.tex")
            print("\n\n\n\n")
            outputs = FormatTable(impacts, caption=f"Systematic uncertainties in percentage at {sqrtS}", label = f"tab:impacts_{sqrtS}", precision=2)
            print(outputs)
            WriteOutputToText(outputs, f"{outdir}/tables/impacts_{sqrtS}.tex")

        if doCombineYear:
            result2json(filename, "sqrtS_Wplus_ratio_ratiometaratio",   f"{outdir}/json/impacts_Wplus_ratio_sqrtS.json", "nuisance_impact_ratiometapois")
            result2json(filename, "sqrtS_Wminus_ratio_ratiometaratio",  f"{outdir}/json/impacts_Wminus_ratio_sqrtS.json",   "nuisance_impact_ratiometapois")
            result2json(filename, "sqrtS_Winc_ratio_ratiometaratio",    f"{outdir}/json/impacts_Winc_ratio_sqrtS.json", "nuisance_impact_ratiometapois")
            result2json(filename, "sqrtS_Zinc_ratio_ratiometaratio",    f"{outdir}/json/impacts_Zinc_ratio_sqrtS.json",  "nuisance_impact_ratiometapois")
            # double ratios
            result2json(filename, "sqrtS_WchgRatio_ratio_doubleratiometaratio",     f"{outdir}/json/impacts_WchgRatio_sqrtS.json", "nuisance_impact_doubleratiometapois")
            result2json(filename, "sqrtS_WZRatio_ratio_doubleratiometaratio",       f"{outdir}/json/impacts_WZRatio_sqrtS.json", "nuisance_impact_doubleratiometapois")

            plotImpacts(f"{outdir}/json/impacts_Wplus_ratio_sqrtS.json",  f"{outdir}/impacts/impacts_Wplus_ratio_mT{idx}_sqrtS")
            plotImpacts(f"{outdir}/json/impacts_Wminus_ratio_sqrtS.json", f"{outdir}/impacts/impacts_Wminus_ratio_mT{idx}_sqrtS")
            plotImpacts(f"{outdir}/json/impacts_Winc_ratio_sqrtS.json",   f"{outdir}/impacts/impacts_Winc_ratio_mT{idx}_sqrtS")
            plotImpacts(f"{outdir}/json/impacts_Zinc_ratio_sqrtS.json",   f"{outdir}/impacts/impacts_Zinc_ratio_mT{idx}_sqrtS")
            plotImpacts(f"{outdir}/json/impacts_WchgRatio_sqrtS.json",    f"{outdir}/impacts/impacts_WchgRatio_mT{idx}_sqrtS")
            plotImpacts(f"{outdir}/json/impacts_WZRatio_sqrtS.json",      f"{outdir}/impacts/impacts_WZRatio_mT{idx}_sqrtS")

            impacts = OrderedDict()
            impacts['Wplus']  = DumpGroupImpacts(filename, "sqrtS_Wplus_ratio_ratiometaratio",   "nuisance_group_impact_ratiometapois")
            impacts['Wminus'] = DumpGroupImpacts(filename, "sqrtS_Wminus_ratio_ratiometaratio",  "nuisance_group_impact_ratiometapois")
            impacts['Winc']   = DumpGroupImpacts(filename, "sqrtS_Winc_ratio_ratiometaratio",    "nuisance_group_impact_ratiometapois")
            impacts['Zinc']   = DumpGroupImpacts(filename, "sqrtS_Zinc_ratio_ratiometaratio",    "nuisance_group_impact_ratiometapois")
            impacts['WOverZ']   = DumpGroupImpacts(filename, "sqrtS_WZRatio_ratio_doubleratiometaratio",     "nuisance_group_impact_doubleratiometapois")
            impacts['WpOverWm'] = DumpGroupImpacts(filename, "sqrtS_WchgRatio_ratio_doubleratiometaratio",   "nuisance_group_impact_doubleratiometapois")

            outputs = FormatTable(impacts, caption=f"Systematic uncertainties in percentage for the cross section ratios between 13TeV and 5TeV", label = f"tab:impacts_sqrtS_mT{idx}", precision=2)
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


if doWpT:
    MakePostPlot("cards/datacard_muplus_lepEta_bin0_WpT.root", "muplus", "WpT", showPULL)

    wpttruthbins  = ["WpT_truth_bin1", "WpT_truth_bin2", "WpT_truth_bin3", "WpT_truth_bin4", "WpT_truth_bin5", "WpT_truth_bin6", "WpT_truth_bin7", "WpT_truth_bin8", "WpT_truth_bin9"]
   
    jsonNames = []
    for wpttruthbin in wpttruthbins:
        result2json("cards/datacard_muplus_lepEta_bin0_WpT.root", "w_muplus_{}_sig_mu".format(wpttruthbin), "cards/impacts_muplus_{}_sig_mu.json".format(wpttruthbin))
        jsonNames.append("cards/impacts_muplus_{}_sig_mu.json".format(wpttruthbin))
        plotImpacts("cards/impacts_muplus_{}_sig_mu.json".format(wpttruthbin), "impacts_muplus_{}".format(wpttruthbin))

    MakeWpTPostFitPlots(jsonNames)


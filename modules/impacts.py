"""
save some pre-computed impacts of some systematics
"""

# from eff calculations on W
unc_effstat_w = {}
unc_effstat_w['effstat_muplus_13TeV'] = 1.0022
unc_effstat_w['effstat_muminus_13TeV'] = 1.0021
unc_effstat_w['effstat_eplus_13TeV'] = 1.0059
unc_effstat_w['effstat_eminus_13TeV'] = 1.0053
unc_effstat_w['effstat_muplus_5TeV'] = 1.0023
unc_effstat_w['effstat_muminus_5TeV'] = 1.0021
unc_effstat_w['effstat_eplus_5TeV'] = 1.0080
unc_effstat_w['effstat_eminus_5TeV'] = 1.0077

# eff stat on Z
unc_effstat_z = {}
unc_effstat_z['effstat_muplus_13TeV'] = 1.0020
unc_effstat_z['effstat_muminus_13TeV'] = 1.0020
unc_effstat_z['effstat_eplus_13TeV'] = 1.0041
unc_effstat_z['effstat_eminus_13TeV'] = 1.0041
unc_effstat_z['effstat_muplus_5TeV'] = 1.0019
unc_effstat_z['effstat_muminus_5TeV'] = 1.0019
unc_effstat_z['effstat_eplus_5TeV'] = 1.0061
unc_effstat_z['effstat_eminus_5TeV'] = 1.0061

# resummation uncertainty
# only affect acceptance, therefore only inclusive xsec
unc_resumm_w = {}
unc_resumm_w["resumm_muplus_13TeV"] = 1.004
unc_resumm_w["resumm_muminus_13TeV"] = 1.004
unc_resumm_w["resumm_eplus_13TeV"] = 1.004
unc_resumm_w["resumm_eminus_13TeV"] = 1.004
unc_resumm_w["resumm_muplus_5TeV"] = 1.003
unc_resumm_w["resumm_muminus_5TeV"] = 1.004
unc_resumm_w["resumm_eplus_5TeV"] = 1.004
unc_resumm_w["resumm_eminus_5TeV"] = 1.004

unc_resumm_z = {}
unc_resumm_z["resumm_mumu_13TeV"] = 1.001
unc_resumm_z["resumm_ee_13TeV"] = 1.001
unc_resumm_z["resumm_mumu_5TeV"] = 1.004
unc_resumm_z["resumm_ee_5TeV"] = 1.004

# FSR uncertainty on W
unc_fsr_w = {}
unc_fsr_w["fsr_muplus_13TeV"] = 1.0014
unc_fsr_w["fsr_muminus_13TeV"] = 1.0006
unc_fsr_w["fsr_eplus_13TeV"] = 1.0033
unc_fsr_w["fsr_eminus_13TeV"] = 1.0025
unc_fsr_w["fsr_muplus_5TeV"] = 1.0014
unc_fsr_w["fsr_muminus_5TeV"] = 1.0006
unc_fsr_w["fsr_eplus_5TeV"] = 1.0033
unc_fsr_w["fsr_eminus_5TeV"] = 1.0025

unc_fsr_z = {}
unc_fsr_z["fsr_mumu_13TeV"] = 1.0005
unc_fsr_z["fsr_ee_13TeV"] = 1.0008
unc_fsr_z["fsr_mumu_5TeV"] = 1.0005
unc_fsr_z["fsr_ee_5TeV"] = 1.0008
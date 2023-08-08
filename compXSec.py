import numpy as np
from collections import OrderedDict
"""
theory results

from Giulio Lenzi

ppbar:

Energy              Z                                 W-                                W+
500        45868.10 + 1647.86 - 1138.33       235396.01 + 8221.80 - 6494.89
235547.50 + 7810.82 - 7008.10
700        77325.39 + 2489.40 - 1832.86       398429.20 + 12360.64 - 9678.81
398399.97 + 12183.04 - 9649.49
1000       119953.55 + 3351.49 - 2753.88      626677.43 + 16986.85 - 12963.71
627224.15 + 15250.08 - 14680.24
1400       173615.31 + 4141.43 - 3621.40      919368.68 + 20593.77 - 17527.57
919645.06 + 19653.62 - 18071.74
3000       392940.24 + 8029.79 - 5916.10      2104471.95 + 40384.46 - 33948.49
2103337.73 + 44050.51 - 29759.76

pp:
Energy              Z                                 W-                                 W+
3000       338131.34 + 7298.24 - 5763.35     1443365.36 + 32699.45 - 28515.93
2322984.02 + 52428.03 - 37592.27
5000       640346.62 + 13007.13 - 9767.54    2791283.06 + 62491.67 - 49378.95
4138688.76 + 89912.67 - 66157.58
7000       949277.02 + 20003.11 - 13375.32   4180328.21 + 84106.48 - 70776.93
5935443.82 + 125116.95 - 97107.21
10000      1416441.13 + 28962.92 - 20799.72  6280739.78 + 115056.96 - 111296.06
8579164.86 + 165284.75 - 149423.39
14000      2037515.51 + 39964.02 - 33317.92  9062406.16 + 183981.50 - 152528.91
11986835.40 + 254025.44 - 194849

"""
xsecs_Z_ppbar = OrderedDict()
xsecs_Z_ppbar[500] = (49.0, 0., 0.)
xsecs_Z_ppbar[700] = (81.3, 0., 0.)
xsecs_Z_ppbar[1000] = (123.0, 0., 0.)
xsecs_Z_ppbar[1400] = (176.0, 0., 0.)
xsecs_Z_ppbar[3000] = (393.0, 0., 0.)

xsecs_Winc_ppbar = OrderedDict()
xsecs_Winc_ppbar[500] = (509.0, 0., 0.)
xsecs_Winc_ppbar[700] = (839.0, 0., 0.)
xsecs_Winc_ppbar[1000] = (1290.0, 0., 0.)
xsecs_Winc_ppbar[1400] = (1870.0, 0., 0.)
xsecs_Winc_ppbar[3000] = (4270.0, 0., 0.)

xsecs_Z_pp = OrderedDict()
xsecs_Z_pp[3000] = (338.0, 0., 0.)
xsecs_Z_pp[5000] = (639.0, 0., 0.)
xsecs_Z_pp[7000] = (946.0, 0., 0.)
xsecs_Z_pp[8000] = (1100.0, 0., 0.)
xsecs_Z_pp[10000] = (1410.0, 0., 0.)
xsecs_Z_pp[13000] = (1870.0, 0., 0.)
xsecs_Z_pp[14000] = (2020.0, 0., 0.)
xsecs_Z_pp[20000] = (2920.0, 0., 0.)

xsecs_Wp_pp = OrderedDict()
xsecs_Wp_pp[3000] = (2390.0, 0., 0.)
xsecs_Wp_pp[5000] = (4540.0, 0., 0.)
xsecs_Wp_pp[7000] = (6070.0, 0., 0.)
xsecs_Wp_pp[8000] = (6960.0, 0., 0.)
xsecs_Wp_pp[10000] = (8730.0, 0., 0.)
xsecs_Wp_pp[13000] = (11330.0, 0., 0.)
xsecs_Wp_pp[14000] = (12170.0, 0., 0.)
xsecs_Wp_pp[20000] = (17180.0, 0., 0.)

xsecs_Wm_pp = OrderedDict()
xsecs_Wm_pp[3000] = (1450.0, 0., 0.)
xsecs_Wm_pp[5000] = (2810.0, 0., 0.)
xsecs_Wm_pp[7000] = (4190.0, 0., 0.)
xsecs_Wm_pp[8000] = (4890.0, 0., 0.)
xsecs_Wm_pp[10000] = (6290.0, 0., 0.)
xsecs_Wm_pp[13000] = (8370.0, 0., 0.)
xsecs_Wm_pp[14000] = (9060.0, 0., 0.)
xsecs_Wm_pp[20000] = (13150.0, 0., 0.)

xsecs_Winc_pp = OrderedDict()
for key in xsecs_Wp_pp:
    xsecs_Winc_pp[key] = (xsecs_Wp_pp[key][0] + xsecs_Wm_pp[key][0], 0., 0.)
    
def ThXSec2TGraph(xsecs, markerstyle = 20, color = 4):
    n = len(xsecs)
    sqrtSs = np.array([key/1.0e3 for key in xsecs.keys()], dtype='d')
    values = np.array([xsecs[key][0] for key in xsecs.keys()], dtype='d')
    print("sqrtSs = ", sqrtSs)
    print("values = ", values)
    g = ROOT.TGraph(n, sqrtSs, values)
    g.SetLineColor(color)
    g.SetLineWidth(4)
    return g

def ExpXsec2TGraph(xsecs, markerstyle = 20, color = 1):
    n = len(xsecs)
    sqrts = xsecs.sqrts
    vals = []
    errs = []
    for xsec in xsecs.xsecs:
        val, err = xsec.GetXsec()
        vals.append(val)
        errs.append(err)
    sqrtSs = np.array([sqrts for i in range(n)], dtype='d')
    vals = np.array(vals, dtype='d')
    errs = np.array(errs, dtype='d')
    g = ROOT.TGraphErrors(n, sqrtSs, vals, np.zeros(n, dtype='d'), errs)
    g.SetMarkerColor(color)
    g.SetLineWidth(2)
    g.SetLineColor(color)
    g.SetMarkerStyle(markerstyle)
    g.SetMarkerSize(1.5)
    return g
        
    
import ROOT
ROOT.gROOT.SetBatch(True)
from CMSPLOTS.myFunction import DrawHistos
gZ_ppbar = ThXSec2TGraph(xsecs_Z_ppbar)
gW_ppbar = ThXSec2TGraph(xsecs_Winc_ppbar)

gZ_pp = ThXSec2TGraph(xsecs_Z_pp)
gWp_pp = ThXSec2TGraph(xsecs_Wp_pp)
gWm_pp = ThXSec2TGraph(xsecs_Wm_pp)
gW_pp = ThXSec2TGraph(xsecs_Winc_pp)

toDraws = [gZ_ppbar, gW_ppbar, gZ_pp, gWp_pp, gWm_pp, gW_pp]
drawoptions = ["PL", "PL", "PL", "PL", "PL", "PL"]
legends = ["" for i in range(len(toDraws))]
legendoptions = ["L", "L", "L", "L", "L", "L"]

from expXsecResults import *
g_CMS_13TeV = ExpXsec2TGraph(xsecs_CMS_13TeV, color = 2, markerstyle = 43)
g_CMS_5TeV = ExpXsec2TGraph(xsecs_CMS_5TeV, color = 2, markerstyle = 45)

toDraws += [g_CMS_13TeV]
drawoptions += ["PE"]
legends += [f"CMS, {xsecs_CMS_13TeV.GetLegend()}"]
legendoptions += ["P"]

g_CMS_8TeV = ExpXsec2TGraph(xsecs_CMS_8TeV)
g_CMS_7TeV = ExpXsec2TGraph(xsecs_CMS_7TeV, color = 1, markerstyle = 24)

toDraws += [g_CMS_8TeV, g_CMS_7TeV]
drawoptions += ["PE", "PE"]
legendoptions += ["P", "P"]
legends += [f"CMS, {xsecs_CMS_8TeV.GetLegend()}", f"CMS, {xsecs_CMS_7TeV.GetLegend()}"]

xsecs_ATLAS_7TeV.sqrts = 7.0-0.2
g_ATLAS_7TeV = ExpXsec2TGraph(xsecs_ATLAS_7TeV, color = 1, markerstyle = 41)

g_CMS_5TeV = ExpXsec2TGraph(xsecs_CMS_5TeV, color = 2, markerstyle = 45)
toDraws += [g_CMS_5TeV]
drawoptions += ["PE"]
legends += [f"CMS, {xsecs_CMS_5TeV.GetLegend()}"]
legendoptions += ["P"]

g_CDF_Run2 = ExpXsec2TGraph(xsecs_CDF_Run2, color = 1, markerstyle = 21)
g_D0_Run1 = ExpXsec2TGraph(xsecs_D0_Run1, color = 1, markerstyle = 25)

toDraws += [g_CDF_Run2, g_D0_Run1]
drawoptions += ["PE", "PE"]
legendoptions += ["P", "P"]
legends += [f"CDF, {xsecs_CDF_Run2.GetLegend()}", f"D0, {xsecs_D0_Run1.GetLegend()}"]

g_UA2 = ExpXsec2TGraph(xsecs_UA2, color = 1, markerstyle = 22)
# shift the x value so that they do not overlap
xsecs_UA1.sqrts = xsecs_UA1.sqrts - 0.035
g_UA1 = ExpXsec2TGraph(xsecs_UA1, color = 1, markerstyle = 32)
# for the label, change it back
xsecs_UA1.sqrts = xsecs_UA1.sqrts + 0.035

toDraws += [g_UA2, g_UA1]
drawoptions += ["PE", "PE"]
legendoptions += ["P", "P"]
legends += [f"UA2, {xsecs_UA2.GetLegend()}", f"UA1, {xsecs_UA1.GetLegend()}"]

Theory = ROOT.TPaveText(2., 40, 24., 80)
Theory.AddText("Theory: NNLO, DYTURBO and NNPDF 3.1 PDFs")
Theory.SetTextAlign(12)
Theory.SetFillColor(4000)
Theory.SetTextColor(4)
Theory.SetBorderSize(0)

ppbar = ROOT.TPaveText(1., 280, 1.4, 500)
ppbar.AddText("p#bar{p}")
ppbar.SetTextAlign(12)
ppbar.SetFillColor(4000)
ppbar.SetTextColor(4)
ppbar.SetBorderSize(0)

pp = ROOT.TPaveText(4., 900, 6., 1700)
pp.AddText("pp")
pp.SetTextAlign(12)
pp.SetFillColor(4000)
pp.SetTextColor(4)
pp.SetBorderSize(0)

scale1_ = 1.0 / 800.0
txts = []
for energy in [0.5, 1, 2, 5, 7, 10, 20]:
    txt_ = ROOT.TLatex()
    xlab_ = 40.
    ytxt_ = 16.
    txt_.SetTextAngle(0)
    txt_.SetTextAlign(23)
    txt_.SetTextFont(42)
    txt_.SetTextSize(xlab_ * scale1_)
    txt_.SetText(energy, ytxt_, str(energy))
    txts.append(txt_)
    
ylabels = {40000: "W^{#pm}", 24000: "W^{+}", 13000: "W^{-}", 3000: "Z"}
for pos, ch in ylabels.items():
    txt_ = ROOT.TLatex()
    xlab_ = 40.
    ytxt_ = pos
    txt_.SetTextAngle(0)
    txt_.SetTextAlign(23)
    txt_.SetTextFont(42)
    txt_.SetTextSize(xlab_ * scale1_)
    txt_.SetText(26, ytxt_, ch)
    txts.append(txt_)


DrawHistos(toDraws, legends, 0.3, 48, "Center-of-mass energy [TeV]", 20, 2e5, "#sigma x B [pb]", "xsec", dologx=True, W_ref = 800, H_ref = 600, noCMS=True, noLumi=True, additionalToDraw = [Theory, ppbar, pp] + txts, drawoptions = drawoptions, legendoptions = legendoptions, legendPos = [0.18, 0.59, 0.5, 0.91], legendTextSize = 0.028, nolabel = True)





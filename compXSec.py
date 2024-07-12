import numpy as np
from collections import OrderedDict
from data.FLAG import doPAS
from modules.Utils import GrepXsecs

xsecs = GrepXsecs()

xsecs['pp']['Z'][2000] = 207.4
xsecs['pp']['Z'][3000] = 361.6
xsecs['pp']['Z'][4000] = 520.5
xsecs['pp']['Z'][5000] = 681.6
xsecs['pp']['Z'][5020] = 684.8
xsecs['pp']['Z'][6000] = 843.6
xsecs['pp']['Z'][7000] = 1006.0
xsecs['pp']['Z'][8000] = 1168.0
xsecs['pp']['Z'][10000] = 1492.0
xsecs['pp']['Z'][13000] = 1975.0
xsecs['pp']['Z'][14000] = 2133.0
xsecs['pp']['Z'][20000] = 3073.0

xsecs['pp']['Winc'] = OrderedDict()
for key in xsecs['pp']['Wp'].keys():
    xsecs['pp']['Winc'][key] = xsecs['pp']['Wp'][key] + xsecs['pp']['Wm'][key]

xsecs['ppbar']['Winc'] = OrderedDict()
for key in xsecs['ppbar']['Wp'].keys():
    xsecs['ppbar']['Winc'][key] = xsecs['ppbar']['Wp'][key] + xsecs['ppbar']['Wm'][key]
    
print("xsecs: ", xsecs)

def ThXSec2TGraph(xsecs, color = 4, name = None, title = None):
    n = len(xsecs)
    sqrtSs = np.array([key/1.0e3 for key in xsecs.keys()], dtype='d')
    values = np.array([xsecs[key] for key in xsecs.keys()], dtype='d')
    print("sqrtSs = ", sqrtSs)
    print("values = ", values)
    g = ROOT.TGraph(n, sqrtSs, values)
    g.SetLineColor(color)
    g.SetLineWidth(4)
    if name:
        g.SetName(name)
    if title:
        g.SetTitle(title)
    return g

def ExpXsec2TGraph(xsecs, markerstyle = 20, color = 1, name = None, title = None):
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
    g.SetLineWidth(1)
    g.SetLineColor(color)
    g.SetMarkerStyle(markerstyle)
    g.SetMarkerSize(1.5)
    if name:
        g.SetName(name)
    if title:
        g.SetTitle(title)
    return g
        
    
import ROOT
ROOT.gROOT.SetBatch(True)
from CMSPLOTS.myFunction import DrawHistos
gZ_ppbar = ThXSec2TGraph(xsecs['ppbar']['Z'], name = "gZ_ppbar")
gW_ppbar = ThXSec2TGraph(xsecs['ppbar']['Winc'], name = "gW_ppbar")

gZ_pp = ThXSec2TGraph(xsecs['pp']['Z'], name = "gZ_pp")
gWp_pp = ThXSec2TGraph(xsecs['pp']['Wp'], name = "gWp_pp")
gWm_pp = ThXSec2TGraph(xsecs['pp']['Wm'], name = "gWm_pp")
gW_pp = ThXSec2TGraph(xsecs['pp']['Winc'], name = "gW_pp")

toDraws = [gZ_ppbar, gW_ppbar, gZ_pp, gWp_pp, gWm_pp, gW_pp]
drawoptions = ["PL", "PL", "PL", "PL", "PL", "PL"]
legends = ["" for i in range(len(toDraws))]
legendoptions = ["L", "L", "L", "L", "L", "L"]

from data.expXsecResults import *
g_CMS_13TeV = ExpXsec2TGraph(xsecs_CMS_13TeV, color = 2, markerstyle = 43, name = "g_CMS_13TeV", title = f"CMS, {xsecs_CMS_13TeV.GetLegend()}")
#g_CMS_5TeV = ExpXsec2TGraph(xsecs_CMS_5TeV, color = 2, markerstyle = 45, name = "g_CMS_5TeV")

toDraws += [g_CMS_13TeV]
drawoptions += ["PZ"]
legends += [f"CMS, {xsecs_CMS_13TeV.GetLegend()}"]
legendoptions += ["P"]

g_CMS_8TeV = ExpXsec2TGraph(xsecs_CMS_8TeV, name = "g_CMS_8TeV", title = f"CMS, {xsecs_CMS_8TeV.GetLegend()}")
g_CMS_7TeV = ExpXsec2TGraph(xsecs_CMS_7TeV, color = 1, markerstyle = 24, name = "g_CMS_7TeV", title = f"CMS, {xsecs_CMS_7TeV.GetLegend()}")

toDraws += [g_CMS_8TeV, g_CMS_7TeV]
drawoptions += ["PZ", "PZ"]
legendoptions += ["P", "P"]
legends += [f"CMS, {xsecs_CMS_8TeV.GetLegend()}", f"CMS, {xsecs_CMS_7TeV.GetLegend()}"]

xsecs_ATLAS_7TeV.sqrts = 7.0-0.2
g_ATLAS_7TeV = ExpXsec2TGraph(xsecs_ATLAS_7TeV, color = 1, markerstyle = 41, name = "g_ATLAS_7TeV", title = f"ATLAS, {xsecs_ATLAS_7TeV.GetLegend()}")

g_CMS_5TeV = ExpXsec2TGraph(xsecs_CMS_5TeV, color = 2, markerstyle = 45, name = "g_CMS_5TeV", title = f"CMS, {xsecs_CMS_5TeV.GetLegend()}")
toDraws += [g_CMS_5TeV]
drawoptions += ["PZ"]
legends += [f"CMS, {xsecs_CMS_5TeV.GetLegend()}"]
legendoptions += ["P"]

g_CMS_2760GeV_2 = ExpXsec2TGraph(xsecs_CMS_2760GeV_2, color = 1, markerstyle = 47, name = "g_CMS_2760GeV_2", title = f"CMS, {xsecs_CMS_2760GeV_2.GetLegend()}")
toDraws += [g_CMS_2760GeV_2]
drawoptions += ["PZ"]
legends += [f"CMS, {xsecs_CMS_2760GeV_2.GetLegend()}"]
legendoptions += ["P"]

g_CMS_2760GeV = ExpXsec2TGraph(xsecs_CMS_2760GeV, color = 1, markerstyle = 46, name = "g_CMS_2760GeV", title = f"CMS, {xsecs_CMS_2760GeV.GetLegend()}")
toDraws += [g_CMS_2760GeV]
drawoptions += ["PZ"]
legends += [f"CMS, {xsecs_CMS_2760GeV.GetLegend()}"]
legendoptions += ["P"]

g_CDF_Run2 = ExpXsec2TGraph(xsecs_CDF_Run2, color = 1, markerstyle = 21, name = "g_CDF_Run2", title = f"CDF, {xsecs_CDF_Run2.GetLegend()}")
g_D0_Run1 = ExpXsec2TGraph(xsecs_D0_Run1, color = 1, markerstyle = 25, name = "g_D0_Run1", title = f"D0, {xsecs_D0_Run1.GetLegend()}")

toDraws += [g_CDF_Run2, g_D0_Run1]
drawoptions += ["PZ", "PZ"]
legendoptions += ["P", "P"]
legends += [f"CDF, {xsecs_CDF_Run2.GetLegend()}", f"D0, {xsecs_D0_Run1.GetLegend()}"]

g_UA2 = ExpXsec2TGraph(xsecs_UA2, color = 1, markerstyle = 22, name = "g_UA2", title = f"UA2, {xsecs_UA2.GetLegend()}")
# shift the x value so that they do not overlap
xsecs_UA1.sqrts = xsecs_UA1.sqrts - 0.035
g_UA1 = ExpXsec2TGraph(xsecs_UA1, color = 1, markerstyle = 32, name = "g_UA1", title = f"UA1, {xsecs_UA1.GetLegend()}")
# for the label, change it back
xsecs_UA1.sqrts = xsecs_UA1.sqrts + 0.035

toDraws += [g_UA2, g_UA1]
drawoptions += ["PZ", "PZ"]
legendoptions += ["P", "P"]
legends += [f"UA2, {xsecs_UA2.GetLegend()}", f"UA1, {xsecs_UA1.GetLegend()}"]

Theory = ROOT.TPaveText(2., 40, 30, 200)
Theory.SetTextSize(0.030)
Theory.SetTextFont(42)
Theory.AddText("Theory: NNLO, DYTURBO and NNPDF 4.0 PDFs")
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
    
ylabels = {40000: "W", 24000: "W^{+}", 13000: "W^{-}", 3000: "Z"}
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


DrawHistos(toDraws, legends, 0.3, 48, "#sqrt{s} [TeV]", 20, 3e5, "#sigma x B [pb]", "xsec", dologx=True, W_ref = 800, H_ref = 600, noLumi=True, additionalToDraw = [Theory, ppbar, pp] + txts, drawoptions = drawoptions, legendoptions = legendoptions, legendPos = [0.18, 0.55, 0.5, 0.91], legendTextSize = 0.028, doPAS = doPAS, nolabel = True, noSqrtS=True)


# CMS only
toDraws_CMSOnly = toDraws[2:-4]
drawoptions = drawoptions[2:-4]
legends = legends[2:-4]
legends = [legend.replace("CMS, ", "") for legend in legends]
legendoptions = legendoptions[2:-4]

txts2 = []
for energy in [2, 5, 7, 10, 20]:
    txt_ = ROOT.TLatex()
    xlab_ = 40.
    ytxt_ = 88.
    txt_.SetTextAngle(0)
    txt_.SetTextAlign(23)
    txt_.SetTextFont(42)
    txt_.SetTextSize(xlab_ * scale1_)
    txt_.SetText(energy, ytxt_, str(energy))
    txts2.append(txt_)
    
ylabels = {40000: "W", 23000: "W^{+}", 14000: "W^{-}", 3500: "Z"}
for pos, ch in ylabels.items():
    txt_ = ROOT.TLatex()
    xlab_ = 40.
    ytxt_ = pos
    txt_.SetTextAngle(0)
    txt_.SetTextAlign(23)
    txt_.SetTextFont(42)
    txt_.SetTextSize(xlab_ * scale1_)
    txt_.SetText(23, ytxt_, ch)
    txts2.append(txt_)
    
Theory2 = ROOT.TPaveText(4.6, 5e4, 27., 3e5)
Theory2.AddText("Theory: NNLO, DYTURBO and NNPDF 4.0 PDFs")
Theory2.SetTextAlign(12)
Theory2.SetFillColor(4000)
Theory2.SetTextColor(4)
Theory2.SetFillColor(0)
Theory2.SetBorderSize(0)

DrawHistos(toDraws_CMSOnly, legends, 1.5, 29, "#sqrt{s} [TeV]", 1.1e2, 3e5, "#sigma x B [pb]", "xsec_CMS_Only_PAS", dologx=True, W_ref = 800, H_ref = 600, noLumi=True, additionalToDraw = [Theory2] + txts2, drawoptions = drawoptions, legendoptions = legendoptions, legendPos = [0.18, 0.64, 0.5, 0.91], legendTextSize = 0.028, doPAS = doPAS, nolabel = True, noSqrtS=True)


# save the graphs for HEPData
fout = ROOT.TFile("xsecs.root", "RECREATE")
for g in toDraws:
    g.Write()
fout.Close()
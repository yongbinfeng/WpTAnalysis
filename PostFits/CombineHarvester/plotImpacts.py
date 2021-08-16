"""
script modified from https://raw.githubusercontent.com/cms-analysis/CombineHarvester/master/CombineTools/scripts/plotImpacts.py
used to plot the impacts from a json file
"""
import ROOT
import math
import json
import argparse
import plotting as plot
import rounding as rounding

perPage = 30
maxPages = None
cmsLabel = "Internal"
checkBoxes = False
colorGroups = None

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.TH1.AddDirectory(0)

def GetRounded(nom, e_hi, e_lo):
    if e_hi < 0.0:
        e_hi = 0.0
    if e_lo < 0.0:
        e_lo = 0.0
    rounded = rounding.PDGRoundAsym(nom, e_hi if e_hi != 0.0 else 1.0, e_lo if e_lo != 0.0 else 1.0)
    s_nom = rounding.downgradePrec(rounded[0],rounded[2])
    s_hi = rounding.downgradePrec(rounded[1][0][0],rounded[2]) if e_hi != 0.0 else '0'
    s_lo = rounding.downgradePrec(rounded[1][0][1],rounded[2]) if e_lo != 0.0 else '0'
    return (s_nom, s_hi, s_lo)

def plotImpacts(ifilename, ofilename, blind=True):
    ROOT.gStyle.Reset()
    # Load the json output of combineTool.py -M Impacts
    data = {}
    with open(ifilename) as jsonfile:
        data = json.load(jsonfile)
    
    # Set the global plotting style
    plot.ModTDRStyle(l=0.4, b=0.10, width=(900 if checkBoxes else 700), height=600)
    
    # We will assume the first POI is the one to plot
    POIs = [ele['name'] for ele in data['POIs']]
    POI = POIs[0]
    
    for ele in data['POIs']:
        if ele['name'] == POI:
            POI_info = ele
            break
    
    POI_fit = POI_info['fit']
    
    # Sort parameters by largest absolute impact on this POI
    data['params'].sort(key=lambda x: abs(x['impact_%s' % POI]), reverse=True)
    
    if checkBoxes:
        cboxes = data['checkboxes']
    
    # Set the number of parameters per page (show) and the number of pages (n)
    show = perPage
    n = int(math.ceil(float(len(data['params'])) / float(show)))
    if maxPages is not None and maxPages > 0:
        n = maxPages
    
    colors = {
        'Gaussian': 1,
        'Poisson': 8,
        'AsymmetricGaussian': 9,
        'Unconstrained': 39,
        'Unrecognised': 2
    }
    color_hists = {}
    color_group_hists = {}
    
    if colorGroups is not None:
        color_groups = {
            x.split('=')[0]: int(x.split('=')[1]) for x in colorGroups.split(',')
        }
    
    seen_types = set()
    
    for name, col in colors.iteritems():
        color_hists[name] = ROOT.TH1F()
        plot.Set(color_hists[name], FillColor=col, Title=name)
    
    if colorGroups is not None:
        for name, col in colorGroups.iteritems():
            color_group_hists[name] = ROOT.TH1F()
            plot.Set(color_group_hists[name], FillColor=col, Title=name)
    
    for page in xrange(n):
        canv = ROOT.TCanvas("plots/"+ofilename, "plots/"+ofilename)
        n_params = len(data['params'][show * page:show * (page + 1)])
        pdata = data['params'][show * page:show * (page + 1)]
        print '>> Doing page %i, have %i parameters' % (page, n_params)
    
        boxes = []
        for i in xrange(n_params):
            y1 = ROOT.gStyle.GetPadBottomMargin()
            y2 = 1. - ROOT.gStyle.GetPadTopMargin()
            h = (y2 - y1) / float(n_params)
            y1 = y1 + float(i) * h
            y2 = y1 + h
            box = ROOT.TPaveText(0, y1, 1, y2, 'NDC')
            plot.Set(box, TextSize=0.02, BorderSize=0, FillColor=0, TextAlign=12, Margin=0.005)
            if i % 2 == 0:
                box.SetFillColor(18)
            box.AddText('%i' % (n_params - i + page * show))
            box.Draw()
            boxes.append(box)
    
        # Crate and style the pads
        if checkBoxes:
            pads = plot.MultiRatioSplitColumns([0.54, 0.24], [0., 0.], [0., 0.])
            pads[2].SetGrid(1, 0)
        else:
            pads = plot.MultiRatioSplitColumns([0.7], [0.], [0.])
        pads[0].SetGrid(1, 0)
        pads[0].SetTickx(1)
        pads[1].SetGrid(1, 0)
        pads[1].SetTickx(1)
    
        h_pulls = ROOT.TH2F("pulls", "pulls", 6, -2.9, 2.9, n_params, 0, n_params)
        g_pulls = ROOT.TGraphAsymmErrors(n_params)
        g_impacts_hi = ROOT.TGraphAsymmErrors(n_params)
        g_impacts_lo = ROOT.TGraphAsymmErrors(n_params)
        g_check = ROOT.TGraphAsymmErrors()
        g_check_i = 0
    
        max_impact = 0.
    
        text_entries = []
        redo_boxes = []
        for p in xrange(n_params):
            i = n_params - (p + 1)
            pre = pdata[p]['prefit']
            fit = pdata[p]['fit']
            tp = pdata[p]['type']
            seen_types.add(tp)
            if pdata[p]['type'] != 'Unconstrained':
                pre_err_hi = (pre[2] - pre[1])
                pre_err_lo = (pre[1] - pre[0])
    
                pull = fit[1] - pre[1]
                pull = (pull/pre_err_hi) if pull >= 0 else (pull/pre_err_lo)
                pull_hi = fit[2] - pre[1]
                pull_hi = (pull_hi/pre_err_hi) if pull_hi >= 0 else (pull_hi/pre_err_lo)
                pull_hi = pull_hi - pull
                pull_lo = fit[0] - pre[1]
                pull_lo = (pull_lo/pre_err_hi) if pull_lo >= 0 else (pull_lo/pre_err_lo)
                pull_lo =  pull - pull_lo
    
                g_pulls.SetPoint(i, pull, float(i) + 0.5)
                g_pulls.SetPointError(
                    i, pull_lo, pull_hi, 0., 0.)
            else:
                # Hide this point
                g_pulls.SetPoint(i, 0., 9999.)
                y1 = ROOT.gStyle.GetPadBottomMargin()
                y2 = 1. - ROOT.gStyle.GetPadTopMargin()
                x1 = ROOT.gStyle.GetPadLeftMargin()
                h = (y2 - y1) / float(n_params)
                y1 = y1 + ((float(i)+0.5) * h)
                x1 = x1 + (1 - pads[0].GetRightMargin() -x1)/2.
                s_nom, s_hi, s_lo = GetRounded(fit[1], fit[2] - fit[1], fit[1] - fit[0])
                text_entries.append((x1, y1, '%s^{#plus%s}_{#minus%s}' % (s_nom, s_hi, s_lo)))
                redo_boxes.append(i)
            g_impacts_hi.SetPoint(i, 0, float(i) + 0.5)
            g_impacts_lo.SetPoint(i, 0, float(i) + 0.5)
            if checkBoxes:
                pboxes = pdata[p]['checkboxes']
                for pbox in pboxes:
                    cboxes.index(pbox)
                    g_check.SetPoint(g_check_i, cboxes.index(pbox) + 0.5, float(i) + 0.5)
                    g_check_i += 1
            imp = pdata[p][POI]
            g_impacts_hi.SetPointError(i, 0, imp[2] - imp[1], 0.5, 0.5)
            g_impacts_lo.SetPointError(i, imp[1] - imp[0], 0, 0.5, 0.5)
            max_impact = max(
                max_impact, abs(imp[1] - imp[0]), abs(imp[2] - imp[1]))
            col = colors.get(tp, 2)
            if colorGroups is not None and len(pdata[p]['groups']) >= 1:
                for p_grp in pdata[p]['groups']:
                    if p_grp in color_groups:
                        col = color_groups[p_grp]
                        break
            h_pulls.GetYaxis().SetBinLabel(
                i + 1, ('#color[%i]{%s}'% (col, pdata[p]['name'])))
    
        # Style and draw the pulls histo
        plot.Set(h_pulls.GetXaxis(), TitleSize=0.04, LabelSize=0.03, Title='(#hat{#theta}-#theta_{0})/#Delta#theta')
    
        plot.Set(h_pulls.GetYaxis(), LabelSize=0.021, TickLength=0.0)
        h_pulls.GetYaxis().LabelsOption('v')
        h_pulls.Draw()
    
        for i in redo_boxes:
            newbox = boxes[i].Clone()
            newbox.Clear()
            newbox.SetY1(newbox.GetY1()+0.005)
            newbox.SetY2(newbox.GetY2()-0.005)
            newbox.SetX1(ROOT.gStyle.GetPadLeftMargin()+0.001)
            newbox.SetX2(0.7-0.001)
            newbox.Draw()
            boxes.append(newbox)
        latex = ROOT.TLatex()
        latex.SetNDC()
        latex.SetTextFont(42)
        latex.SetTextSize(0.02)
        latex.SetTextAlign(22)
        for entry in text_entries:
            latex.DrawLatex(*entry)
    
        # Go to the other pad and draw the impacts histo
        pads[1].cd()
        if max_impact == 0.: max_impact = 1E-6  # otherwise the plotting gets screwed up
        h_impacts = ROOT.TH2F(
            "impacts", "impacts", 6, -max_impact * 1.1, max_impact * 1.1, n_params, 0, n_params)
        plot.Set(h_impacts.GetXaxis(), LabelSize=0.03, TitleSize=0.04, Ndivisions=505, Title=
            '#Delta#hat{%s}' % (POI))
        plot.Set(h_impacts.GetYaxis(), LabelSize=0, TickLength=0.0)
        h_impacts.Draw()
    
        if checkBoxes:
            pads[2].cd()
            h_checkboxes = ROOT.TH2F(
                "checkboxes", "checkboxes", len(cboxes), 0, len(cboxes), n_params, 0, n_params)
            for i, cbox in enumerate(cboxes):
                h_checkboxes.GetXaxis().SetBinLabel(i+1, cbox)
            plot.Set(h_checkboxes.GetXaxis(), LabelSize=0.03, LabelOffset=0.002)
            h_checkboxes.GetXaxis().LabelsOption('v')
            plot.Set(h_checkboxes.GetYaxis(), LabelSize=0, TickLength=0.0)
            h_checkboxes.Draw()
            # g_check.SetFillColor(ROOT.kGreen)
            g_check.Draw('PSAME')
    
        # Back to the first pad to draw the pulls graph
        pads[0].cd()
        plot.Set(g_pulls, MarkerSize=0.8, LineWidth=2)
        g_pulls.Draw('PSAME')
    
        # And back to the second pad to draw the impacts graphs
        pads[1].cd()
        alpha = 0.7
    
        lo_color = {
            'default': 38,
            'hesse': ROOT.kOrange - 3,
            'robust': ROOT.kGreen + 1
        }
        hi_color = {
            'default': 46,
            'hesse': ROOT.kBlue,
            'robust': ROOT.kAzure - 5
        }
        method = 'default'
        if 'method' in data and data['method'] in lo_color:
            method = data['method']
        g_impacts_hi.SetFillColor(plot.CreateTransparentColor(hi_color[method], alpha))
        g_impacts_hi.Draw('2SAME')
        g_impacts_lo.SetFillColor(plot.CreateTransparentColor(lo_color[method], alpha))
        g_impacts_lo.Draw('2SAME')
        pads[1].RedrawAxis()
    
        legend = ROOT.TLegend(0.02, 0.02, 0.40, 0.06, '', 'NBNDC')
        legend.SetNColumns(3)
        legend.AddEntry(g_pulls, 'Pull', 'LP')
        legend.AddEntry(g_impacts_hi, '+1#sigma Impact', 'F')
        legend.AddEntry(g_impacts_lo, '-1#sigma Impact', 'F')
        legend.Draw()
    
        leg_width = pads[0].GetLeftMargin() - 0.01
        if colorGroups is not None:
            legend2 = ROOT.TLegend(0.01, 0.94, leg_width, 0.99, '', 'NBNDC')
            legend2.SetNColumns(2)
            for name, h in color_group_hists.iteritems():
                legend2.AddEntry(h, name, 'F')
            legend2.Draw()
        elif len(seen_types) > 1:
            legend2 = ROOT.TLegend(0.01, 0.94, leg_width, 0.99, '', 'NBNDC')
            legend2.SetNColumns(2)
            for name, h in color_hists.iteritems():
                if name == 'Unrecognised': continue
                legend2.AddEntry(h, name, 'F')
            legend2.Draw()
    
        plot.DrawCMSLogo(pads[0], 'CMS', cmsLabel, 0, 0.25, 0.00, 0.00)
        s_nom, s_hi, s_lo = GetRounded(POI_fit[1], POI_fit[2] - POI_fit[1], POI_fit[1] - POI_fit[0])
        if not blind:
            plot.DrawTitle(pads[1], '#hat{%s} = %s^{#plus%s}_{#minus%s}%s' % (
                POI, s_nom, s_hi, s_lo,
                ''), 3, 0.27)
        extra = ''
        if page == 0:
            extra = '('
        if page == n - 1:
            extra = ')'
        canv.Print('.pdf%s' % extra)

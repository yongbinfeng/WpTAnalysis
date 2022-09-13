import ROOT
from ROOT import TH1F, TCanvas, TPad, Math, TF1, TLegend
from . import CMS_lumi
from . import tdrstyle


class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    RED = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    RED = '\033[1;31m'
    UNDERLINE = '\033[4m'


def DumpHist(hist):
    # dump the number of events in histogram
    print("histogram %s with %d bins" % (hist.GetName, hist.GetNbinsX()))
    for ibin in range(0, hist.GetNbinsX()+2):
        print(" %d bin with %.2f events" % (ibin, hist.GetBinContent(ibin)))


def PrepareHisto(myfile, listhistos, printlist=True):
    if printlist:
        print(listhistos)
    myhisto = myfile.Get("%s" % listhistos[0])
    for idx in range(1, len(listhistos)):
        myhisto.Add(myfile.Get("%s" % listhistos[idx]))
    return myhisto


def RebinHisto(histo, *args):
    if len(args) == 0:
        return histo
    elif len(args) == 1:
        return histo.Rebin(int(args[0]))
    else:
        return histo.Rebin(len(args[0])-1, "Rebinned_%s" % args[1], args[0])


def GetHisto(myfile, lhistos, *args):
    htemp = PrepareHisto(myfile, lhistos[0].split("+"))
    htemp_rebin = RebinHisto(htemp, *args)

    if len(lhistos) > 1:
        htemp_den = PrepareHisto(myfile, lhistos[1].split("+"))
        htemp_den_rebin = RebinHisto(htemp_den, *args)
        # deal with the case with numerator is 0 but denominator is not
        # temporary fix, use the upper error as the error
        # https://twiki.cern.ch/twiki/bin/viewauth/CMS/PoissonErrorBars
        for ibin in range(1, htemp_rebin.GetNbinsX()+1):
            if htemp_rebin.GetBinContent(ibin) == 0 and htemp_den_rebin.GetBinContent(ibin) != 0:
                htemp_rebin.SetBinError(
                    ibin, Math.gamma_quantile_c((1 - 0.6827)/2.0, 1, 1))
                htemp_den_rebin.SetBinError(ibin, Math.gamma_quantile_c(
                    (1 - 0.6827)/2.0, htemp_den_rebin.GetBinContent(ibin)+1, 1))
        # divide
        # set the uncertainty to binomial, see line 2800 of https://root.cern.ch/doc/master/TH1_8cxx_source.html
        htemp_rebin.Divide(htemp_rebin, htemp_den_rebin, 1.0, 1.0, "B")
        ##
    return htemp_rebin


def Normalize(th1, option=0):
    if option == 0:
        th1.Scale(1.0/(th1.Integral(0, th1.GetNbinsX()+1)+1e-6))
    elif option == 1:
        th1.Scale(1.0/(th1.Integral()+1e-6))


def ScaleWithWidth(th1):
    th1.Scale(1.0, "width")


def myDivide(num, den):
    # set the ratio to 1 if denominator is 0
    if den == 0:
        #print bcolors.RED, "denominator is zero. Set the ratio to 1", bcolors.ENDC
        return 0.
    else:
        return num/den


def CalculateChi(hobsclone, hexp, doNewman=False, doPearson=False, ignoreHistError=False, printchi2=True):

    if doNewman and doPearson:
        print(bcolors.RED, "can NOT do Newman and Pearson chi2 at the same time!!! Use Pearson instead", bcolors.ENDC)

    for ibin in range(1, hexp.GetNbinsX()+1):
        num = hobsclone.GetBinContent(ibin) - hexp.GetBinContent(ibin)
        if doPearson and ignoreHistError:
            den = ROOT.TMath.Sqrt(hexp.GetBinContent(ibin))
        elif doNewman and ignoreHistError:
            den = ROOT.TMath.Sqrt(hobsclone.GetBinContent(ibin))
        elif doPearson and not ignoreHistError:
            den = hexp.GetBinError(ibin)
        elif doNewman and not ignoreHistError:
            den = hobsclone.GetBinError(ibin)

        #print ibin, hobsclone.GetBinContent(ibin), hexp.GetBinContent(ibin), num, den

        hobsclone.SetBinContent(ibin, myDivide(num, den))
        hobsclone.SetBinError(ibin, 0.)

    if printchi2:
        chi2 = 0
        for ibin in range(1, hobsclone.GetNbinsX()+1):
            chi2 += ROOT.TMath.Power(hobsclone.GetBinContent(ibin), 2.0)
        print("Chi2 for %s is %.2f" % (hobsclone.GetName(), chi2))


'''
Read file, color, label and list of histograms information from the input file
ignore lines starting with #
'''


def myRead(inputfile):

    list_input = open("%s" % inputfile)

    files = []
    colors = []
    labels = []
    lhistos = []

    for line in list_input:
        tline = []
        tline = line.rstrip().split()

        if tline[0][0] == '#':  # ignore the line starting with #
            continue

        if len(tline) < 4 or len(tline) > 5:
            print("Input Error in line\n   %s of file %s" % (line, inputfile))
            sys.exit(0)

        files.append(tline[0])
        colors.append(int(tline[1]))
        labels.append(tline[2].replace("_", " "))

        if len(tline) == 4:
            lhistos.append([tline[3]])
        else:
            lhistos.append([tline[3], tline[4]])

    return files, colors, labels, lhistos


def getResolution(th2d):
    # return the resolution distribution of th2d, as a format of th1d
    hqup = th2d.QuantilesX(0.84, th2d.GetName()+"_q84")
    hqdn = th2d.QuantilesX(0.16, th2d.GetName()+"_q16")

    hresol = hqup.Clone(th2d.GetName() + "_resol")
    hresol.Add(hqdn, -1.0)
    hresol.Scale(0.5)

    return hresol


def getMedian(th2d):
    # return the median of the th2d
    hmedian = th2d.QuantilesX(0.5, th2d.GetName()+"_q50")
    return hmedian


def getMean(th2d):
    hpmean = th2d.ProfileX("_mean")
    hmean = hpmean.ProjectionX()
    return hmean


def getErrors(hprof, verbose=False):
    herr = hprof.ProjectionX("%s_Error" % hprof.GetName())
    # to clone the binning of input hprof. If there is a better way,
    # this can be changed
    if verbose:
        print("input histogram %s has %d bins" % (
            hprof.GetName(), hprof.GetNbinsX()))
    for ibin in range(1, hprof.GetNbinsX()+1):
        if verbose:
            print("%d bin, with RMS %f" % (ibin, hprof.GetBinError(ibin)))
        herr.SetBinContent(ibin, hprof.GetBinError(ibin))
        herr.SetBinError(ibin, 0)
    for ibin in range(1, hprof.GetNbinsX()+1):
        if verbose:
            print("%d bin, with Content %f" % (ibin, herr.GetBinContent(ibin)))
    return herr


def THStack2TH1(hstack, postfix=""):
    assert isinstance(hstack, ROOT.THStack), "input must be a ROOT.THStack"
    hlist = list(hstack.GetHists())
    h1 = hlist[0].Clone("{}_histogram_{}".format(hstack.GetName(), postfix))
    list(map(h1.Add, hlist[1:]))
    return h1

def AddOverflowsTH1(h1, dolastbin=True):
    assert isinstance(h1, ROOT.TH1), "input must be a ROOT.TH1"
    nbins = h1.GetNbinsX()
    if dolastbin:
        # merge overflow bin to the last bin
        h1.SetBinContent(nbins, h1.GetBinContent(nbins)+h1.GetBinContent(nbins+1))
        # treat the uncertainties as uncorrelated
        h1.SetBinError(nbins, ROOT.TMath.Sqrt((h1.GetBinError(nbins))**2 + (h1.GetBinError(nbins+1))**2))
        # clean the old bins
        h1.SetBinContent(nbins + 1, 0)
        h1.SetBinError(nbins + 1, 0)
    else:
        h1.SetBinContent(1, h1.GetBinContent(0)+h1.GetBinContent(1))
        # treat the uncertainties as uncorrelated
        h1.SetBinError(1, ROOT.TMath.Sqrt((h1.GetBinError(1))**2 + (h1.GetBinError(0))**2))
        # clean the old bins
        h1.SetBinContent(0, 0)
        h1.SetBinError(0, 0)

def AddOverflows(hinput, dolastbin=True):
    '''
    move the over/under flow bin to the last/first bin
    '''
    if isinstance(hinput, ROOT.TH1):
        AddOverflowsTH1(hinput, dolastbin)

    if isinstance(hinput, ROOT.THStack):
        # do the AddOverflowsTH1 for all the histograms in THStack
        hlist = list(hinput.GetHists())
        list(map(AddOverflowsTH1, hlist, [dolastbin]*len(hlist)))

def Ratio2Diff(hratio, inpercent=True):
    '''
    convert the ratio plot to diff, by substracting 1.0,
    and may show it in percentile
    '''
    for ibin in range(0,hratio.GetNbinsX()+1):
        hratio.SetBinContent(ibin, hratio.GetBinContent(ibin)-1.0)
    if inpercent:
        hratio.Scale(100.0)


'''
draw histograms with the CMS tdr style
'''


def DrawHistos(myhistos, mylabels, xmin, xmax, xlabel, ymin, ymax, ylabel, outputname, dology=True, showratio=False, dologx=False, lheader=None, donormalize=False, binomialratio=False, yrmax=2.0, yrmin=0.0, yrlabel=None, MCOnly=False, leftlegend=False, mycolors=[], legendPos=[], legendNCols=1, linestyles=[], markerstyles=[], showpull=False, doNewman=False, doPearson=False, ignoreHistError=False, ypullmin=-2.99, ypullmax=2.99, drawashist=False, padsize=(2, 0.9, 1.1), setGridx=False, setGridy=False, drawoptions=[], legendoptions=[], ratiooptions=[], dologz=False, doth2=False, ratiobase=0, redrawihist=-1, extraText=None, noCMS=False, noLumi=False, nMaxDigits=None, addOverflow=False, addUnderflow=False, plotdiff=False, hratiopanel = None, doratios=None, hpulls=None, W_ref = 600, is5TeV = False):
    # set the tdr style
    tdrstyle.setTDRStyle()
    # not sure why need this...
    ROOT.gStyle.SetErrorX(0.5)

    # change the CMS_lumi variables (see CMS_lumi.py)
    # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
    CMS_lumi.lumi_sqrtS = "(13 TeV)"
    CMS_lumi.relPosX = 0.12
    #CMS_lumi.extraText = "Internal"
    #if MCOnly:
    #    CMS_lumi.extraText = "Simulation"
    if isinstance(extraText, str):
        CMS_lumi.extraText = extraText
    if extraText == None:
        CMS_lumi.extraText = ""

    if nMaxDigits:
        ROOT.TGaxis.SetMaxDigits(nMaxDigits)
    else:
        # default val
        ROOT.TGaxis.SetMaxDigits(5)

    H_ref = 500
    W = W_ref
    H = H_ref

    iPeriod = 0

    npads = 1 + showratio + showpull

    if npads == 2:
        canvas = TCanvas("c2"+outputname, "c2", 50, 50, W, 600)
        padsize1 = float(padsize[0])/(padsize[0]+padsize[1])
        padsize2 = float(padsize[1])/(padsize[0]+padsize[1])
        padsize3 = 0.
    elif npads == 1:
        canvas = TCanvas("c2"+outputname, "c2", 50, 50, W, 600)
        canvas.SetGrid(setGridx, setGridy)
        canvas.SetTicks(1,1)
        padsize1 = 1.0
        if doth2:
            canvas.SetRightMargin(0.10)
        padsize2 = 0.
        padsize3 = 0.
        canvas.cd()
    else:
        canvas = TCanvas("c2"+outputname, "c2", 50, 50, W, 800)
        padsize1 = float(padsize[0])/(padsize[0]+padsize[1]+padsize[2])
        padsize2 = float(padsize[1])/(padsize[0]+padsize[1]+padsize[2])
        padsize3 = float(padsize[2])/(padsize[0]+padsize[1]+padsize[2])

    if dology:
        #print "set y axis to log scale"
        canvas.SetLogy()
    if dologx:
        #print "set x axis to log scale"
        canvas.SetLogx()
    if dologz:
        #print "set z axis to log scale"
        canvas.SetLogz()

    if npads == 2:
        pad1 = TPad("pad1" + outputname, "pad1", 0, padsize2, 1, 1)
        pad2 = TPad("pad2" + outputname, "pad1", 0, 0, 1, padsize2)
        pad1.SetTopMargin(0.06/padsize1)
        pad1.SetBottomMargin(0.012/padsize1)
        pad1.SetLeftMargin(0.15 * (600.0)/W)
        pad2.SetTopMargin(0.010/padsize2)
        pad2.SetBottomMargin(0.13/padsize2)
        pad2.SetLeftMargin(0.15 * (600.0)/W)
        pad2.SetGridy(1)
        pad2.SetTicks(1,1)
        pad1.Draw()
        pad2.Draw()

    if npads == 3:
        pad1 = TPad("pad1" + outputname, "pad1", 0, 1-padsize1, 1, 1)
        pad2 = TPad("pad2" + outputname, "pad2", 0, padsize3, 1, 1-padsize1)
        pad3 = TPad("pad3" + outputname, "pad3", 0, 0.,   1, padsize3)
        pad1.SetTopMargin(0.06/(padsize1+padsize3))
        pad1.SetBottomMargin(0.012/padsize1)
        pad1.SetLeftMargin(0.15 * (600.0)/W)
        pad2.SetTopMargin(0.010/padsize2)
        pad2.SetBottomMargin(0.02/padsize2)
        pad2.SetLeftMargin(0.15 * (600.0)/W)
        pad2.SetGridy(1)
        pad2.SetTicks(1,1)
        pad3.SetTopMargin(0.01/padsize3)
        pad3.SetBottomMargin(0.13/padsize3)
        pad3.SetLeftMargin(0.15 * (600.0)/W)
        pad3.SetGridy(1)
        pad3.SetTicks(1,1)
        pad1.Draw()
        pad2.Draw()
        pad3.Draw()

    if npads > 1:
        pad1.cd()
        pad1.SetGrid(setGridx, setGridy)
        pad1.SetTicks(1,1)
        if dology:
            pad1.SetLogy()
        if dologx:
            pad1.SetLogx()

    h1 = TH1F("h1" + outputname, "h1", 80, xmin, xmax)

    h1.SetMinimum(ymin)
    h1.SetMaximum(ymax)

    #print "xmin : %f xmax : %f"%(xmin, xmax)

    h1.GetXaxis().SetNdivisions(6, 5, 0)

    h1.GetYaxis().SetNdivisions(6, 5, 0)
    h1.GetYaxis().SetTitle("%s" % ylabel)
    h1.GetYaxis().SetTitleSize(0.050/(padsize1+padsize3))
    h1.GetYaxis().SetLabelSize(0.045/(padsize1+padsize3))
    h1.GetXaxis().SetTitleSize(0.050/(padsize1+padsize3))
    h1.GetXaxis().SetLabelSize(0.045/(padsize1+padsize3))
    h1.GetYaxis().SetTitleOffset(1.35*(padsize1+padsize3)*(600.0/W))
    h1.GetXaxis().SetTitleOffset(1.1*(padsize1+padsize3))

    if showratio or showpull:
        h1.GetXaxis().SetLabelSize(0)
    else:
        h1.GetXaxis().SetTitle("%s" % xlabel)

    h1.Draw()

    x1_l = 0.92
    y1_l = 0.88

    #dx_l = 0.20
    dy_l = 0.22*0.66/padsize1
    dx_l = 0.35*0.66/padsize1
    x0_l = x1_l-dx_l
    y0_l = y1_l-dy_l

    if len(legendPos) == 4:
        x0_l = legendPos[0]
        y0_l = legendPos[1]
        x1_l = legendPos[2]
        y1_l = legendPos[3]

    if npads > 2:
        y1_l -= 0.03
        y0_l -= 0.03

    if not leftlegend:
        legend = TLegend(x0_l, y0_l, x1_l, y1_l)
    else:
        legend = TLegend(x0_l-0.40, y0_l, x1_l-0.40, y1_l)
    if lheader and lheader != "":
        legend.SetHeader(lheader)
    legend.SetNColumns(legendNCols)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.04)
    legend.SetTextFont(42)
    legend.SetFillColor(0)

    myf = []

    myhistos_clone = []
    for ihisto in myhistos:
        ihcolone = ihisto.Clone("%s_Clone" % ihisto.GetName())
        # ihcolone.SetDirectory(0)
        myhistos_clone.append(ihcolone)

    if drawashist:
        drawoptions = ["HIST" for i in range(len(myhistos_clone))]
        legendoptions = ['L' for i in range(len(myhistos_clone))]

    if not isinstance(drawoptions, list):
        # copy the option n times
        tmp = [drawoptions for i in range(len(myhistos_clone))]
        drawoptions = tmp
    if not isinstance(legendoptions, list):
        tmp = [legendoptions for i in range(len(myhistos_clone))]
        legendoptions = tmp

    #print(("legend options ", legendoptions))

    ileg = 0
    for idx in range(0, len(myhistos_clone)):
        if addOverflow:
            # include overflow bin
            AddOverflows(myhistos_clone[idx])
        if addUnderflow:
            # include underflow bin
            AddOverflows(myhistos_clone[idx], dolastbin=False)
        if len(mycolors) == len(myhistos_clone):
            myhistos_clone[idx].SetLineColor(mycolors[idx])
            myhistos_clone[idx].SetMarkerColor(mycolors[idx])
        if len(linestyles) == len(myhistos_clone):
            myhistos_clone[idx].SetLineStyle(linestyles[idx])
        if len(markerstyles) == len(myhistos_clone):
            myhistos_clone[idx].SetMarkerStyle(markerstyles[idx])
        if isinstance(myhistos_clone[idx], ROOT.TH1):
            myhistos_clone[idx].SetLineWidth(3)
        #myhistos_clone[idx].GetXaxis().SetRangeUser(xmin, xmax)
        if donormalize:
            Normalize(myhistos_clone[idx])
        if idx >= len(drawoptions):
            drawoptions.append("EL")
        if idx >= len(legendoptions):
            if isinstance(myhistos_clone[idx], ROOT.TH1) and myhistos_clone[idx].GetMarkerStyle() != 1:
                legendoptions.append("LP")
            else:
                legendoptions.append("LE")
        if isinstance(myhistos_clone[idx], ROOT.THStack):
            drawoptions[idx] = 'HIST'
            legendoptions[idx] = "F"
        myhistos_clone[idx].Draw(" ".join(filter(None,[drawoptions[idx], "same"])))

        if ileg < len(mylabels):
            # number of labels might be different from number of histograms.
            # do the match logically
            if isinstance(myhistos_clone[idx], ROOT.THStack):
                hlist = list(myhistos_clone[idx].GetHists())
                for hist in reversed(hlist):
                    legend.AddEntry(
                        hist, str(mylabels[ileg]), legendoptions[idx])
                    ileg += 1
            else:
                legend.AddEntry(myhistos_clone[idx], str(
                    mylabels[ileg]), legendoptions[idx])
                ileg += 1

    #print("draw options ", drawoptions)
    #print("legend options ", legendoptions)

    if redrawihist >= 0:
        myhistos_clone[redrawihist].Draw(" ".join(filter(None,[drawoptions[redrawihist], "same"])))

    iPosX = 0
    plotCMS = True
    if noCMS:
        # do not draw CMS, only lumi info and extraText
        plotCMS = False
        CMS_lumi.cmsText = ""
        #CMS_lumi.extraText = ""
        #iPosX = 11
    if not is5TeV:
        iPeriod = 4
    else:
        iPeriod = 5
    if MCOnly:
        iPeriod = 0
    if showratio or showpull:
        if not noLumi:
            CMS_lumi.CMS_lumi(pad1, iPeriod, iPosX, plotCMS = plotCMS)
        pad1.Update()
        pad1.RedrawAxis()
        #frame = pad1.GetFrame()
        #frame.Draw()
        if len(mylabels):
            legend.Draw()
        pad1.Update()

    else:
        if not noLumi:
            CMS_lumi.CMS_lumi(canvas, iPeriod, iPosX, plotCMS = plotCMS)
        canvas.cd()
        canvas.Update()
        canvas.RedrawAxis()
        canvas.RedrawAxis("G")
        #frame = canvas.GetFrame()
        #frame.Draw()
        if len(mylabels):
            #print("draw legend", len(mylabels))
            legend.Draw()
        canvas.Update()

    ##
    # ration plot
    ##
    if showratio or showpull:
        if isinstance(myhistos_clone[ratiobase], ROOT.THStack):
            hden = THStack2TH1(myhistos_clone[ratiobase])
        else:
            hden = myhistos_clone[ratiobase].Clone("hden_%s" % idx)

    if showratio:
        hratios = {}
        for idx in range(0, len(myhistos_clone)):
            if doratios != None and doratios[idx]==False:
                continue
            if idx == ratiobase:
                continue
            if isinstance(myhistos_clone[idx], ROOT.THStack):
                hratios[idx] = THStack2TH1(myhistos_clone[idx])
            else:
                hratios[idx] = myhistos_clone[idx].Clone("hratios_%s" % idx)
            if binomialratio:
                hratios[idx].Divide(hratios[idx], hden, 1.0, 1.0, "b")
            else:
                hratios[idx].Divide(hden)
            if plotdiff:
                Ratio2Diff(hratios[idx])

    if showpull:
        if hpulls == None:
            hpulls = []
            for idx in range(0, len(myhistos_clone)):
                if idx == ratiobase:
                    continue
                hpull = myhistos_clone[idx].Clone("hpullsnew_%s" % idx)
                hpulls.append( hpull )
                CalculateChi(hpull, hden, doNewman=doNewman,
                             doPearson=doPearson, ignoreHistError=ignoreHistError)

    if npads > 1:
        pad2.cd()
        if dologx:
            pad2.SetLogx()

        h2 = TH1F("h2", "h2", 80, xmin, xmax)
        h2.GetXaxis().SetTitleSize(0.050/(padsize2+0.3*padsize3))
        h2.GetXaxis().SetLabelSize(0.045/(padsize2+0.3*padsize3))
        h2.GetYaxis().SetTitleSize(0.050/(padsize2+0.3*padsize3))
        h2.GetYaxis().SetLabelSize(0.045/(padsize2+0.3*padsize3))
        h2.GetYaxis().SetTitleOffset(1.35*(padsize2+0.35*padsize3)*(600.0/W))

        h2.GetYaxis().SetNdivisions(8)
        h2.GetYaxis().CenterTitle()
        if yrlabel:
            ytitle = yrlabel
        elif showratio:
            ytitle = "Ratio"
        else:
            ytitle = "Pull"
        h2.GetYaxis().SetTitle(ytitle)
        h2.SetMaximum(yrmax if showratio else ypullmax)
        h2.SetMinimum(yrmin if showratio else ypullmin)

        if npads == 2:
            h2.GetXaxis().SetTitle("%s" % xlabel)
            h2.GetXaxis().SetTitleOffset(1.1)
        else:
            h2.GetXaxis().SetTitleSize(0.0)
            h2.GetXaxis().SetLabelSize(0)

        h2.Draw()

        if showratio:
            if hratiopanel:
                print("draw hratiopanel")
                #hratiopanel.SetFillColor(15)
                hratiopanel.SetFillColorAlpha(15, 0.5)
                hratiopanel.SetLineColor(2)
                hratiopanel.SetMarkerSize(0)
                hratiopanel.Draw("E2 same")
                #hratiopanel.Draw("HIST same")
                #for ibin in xrange(1, hratiopanel.GetNbinsX()+1):
                #    print("bin content ", hratiopanel.GetBinContent(ibin), " error ", hratiopanel.GetBinError(ibin))
                h2.Draw("same")
            idx = 0
            for hratio in list(hratios.values()):
                if idx >= len(ratiooptions):
                    if drawashist:
                        ratiooptions.append("HIST")
                    else:
                        ratiooptions.append("")
                hratio.SetFillColor(0)
                hratio.Draw(ratiooptions[idx] + " same")
                idx += 1
        elif showpull:
            for hpull in hpulls:
                hpull.Draw("HIST same")

        pad2.RedrawAxis("G")
        pad2.Update()

    if npads > 2:
        pad3.cd()
        if dologx:
            pad3.SetLogx()

        h3 = TH1F("h3", "h3", 80, xmin, xmax)
        h3.GetXaxis().SetTitle("%s" % xlabel)
        h3.GetXaxis().SetTitleOffset(1.1)
        h3.GetXaxis().SetTitleSize(0.050/(padsize3+0.3*padsize2))
        h3.GetXaxis().SetLabelSize(0.045/(padsize3+0.3*padsize2))
        h3.GetYaxis().SetTitleSize(0.050/(padsize3+0.3*padsize2))
        h3.GetYaxis().SetLabelSize(0.045/(padsize3+0.3*padsize2))
        h3.GetYaxis().SetTitleOffset(1.35*(padsize3+0.35*padsize2)*(600.0/W))
        h3.GetYaxis().SetNdivisions(8)
        h3.GetYaxis().CenterTitle()
        #h3.GetYaxis().SetTitle("#frac{Obs. - Exp.}{Err}")
        h3.GetYaxis().SetTitle("Pull")

        h3.SetMaximum(ypullmax)
        h3.SetMinimum(ypullmin)
        h3.Draw()

        for hpull in hpulls:
            hpull.SetFillColorAlpha(2, 0.6)
            hpull.SetLineWidth(1)
            hpull.SetLineColor(1)
            hpull.Draw("HIST same")

        pad3.Update()

    canvas.Update()
    print("save plot to plots/%s.pdf" % outputname)
    canvas.Print("plots/%s.C"%outputname)
    canvas.Print("plots/%s.pdf" % outputname)
    #canvas.Print("plots/%s.png" % outputname)
    #canvas.Print("plots/%s.root" % outputname)
    #print("function runs fine")
    return hratios if showratio else None

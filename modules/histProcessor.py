import ROOT
import re
import numpy as np
from CMSPLOTS.myFunction import AddOverflowsTH1, RebinHisto

ROOT.gROOT.SetBatch(True)

def DoRebin(hist, mass_bins):
    """
    rebin the histogram to be used for the mass fit
    """
    hist = RebinHisto(hist, mass_bins, hist.GetName())

    return hist

def SaveHistToFile(hist, ofile):
    """
    post processing for histograms
    """
    hist.SetDirectory(ofile)
    hist.Write()

def ProcessHists(ifile, ofile, mass_bins, includeUnderflow = False, includeOverflow = False):
    """
    process all the hists in ifile (rebinning, include over/underflow, etc),
    and save them into ofile
    """
    finput = ROOT.TFile.Open(ifile)
    hnames = finput.GetListOfKeys()
    hnames = [hname.GetName() for hname in hnames]

    foutput = ROOT.TFile(ofile, "RECREATE")

    for hname in hnames:
        h = finput.Get(hname)
        if includeOverflow:
            AddOverflowsTH1(h)
        if includeUnderflow:
            AddOverflowsTH1(h, False)
        h = DoRebin(h, mass_bins)
        SaveHistToFile(h, foutput)

    foutput.Close()

def CopyandMergeTau(iname, oname):
    """
    copy all the histograms from iname;
    merge the tau processes into the signal process,
    and assign relative systematics;
    save into new file
    """
    finput = ROOT.TFile(iname)
    hnames = finput.GetListOfKeys()
    hnames = [hname.GetName() for hname in hnames]

    # get the substring for wx (tau) and wlnu
    wxstr = None
    wlnustr = None
    for hname in hnames:
        if wxstr is None:
            searched = re.search('wx\d*', hname)
            if searched:
                wxstr = searched.group(0)
        if wlnustr is None:
            searched = re.search('wlnu\d*', hname)
            if searched:
                wlnustr = searched.group(0)
    print("wx string is ", wxstr)
    print("wlnu string is ", wlnustr)

    foutput = ROOT.TFile(oname, "RECREATE")
    for hname in hnames:
        h = finput.Get(hname)

        if wlnustr not in hname and wxstr not in hname:
            SaveHistToFile(h, foutput)

        elif wlnustr in hname:
            print(hname)
            # sig represents both wlnu and wx
            hsig = h.Clone(hname.replace(wlnustr, "wsig"))
            htau = finput.Get(hname.replace(wlnustr, wxstr))
            hsig.Add(htau)
            SaveHistToFile(hsig, foutput)

            if re.search("wlnu\d*$", hname):
                # central histograms
                # assign tau fraction systematic variations on it
                # and save to output
                hsig_up = h.Clone(hsig.GetName() + "_SysTauFracUp")
                hsig_dn = h.Clone(hsig.GetName() + "_SysTauFracDown")

                htau_up = htau.Clone(htau.GetName() + "_SysTauFracUp")
                htau_dn = htau.Clone(htau.GetName() + "_SysTauFracDown")
                varfraction = 0.2
                htau_up.Scale(1.0 + varfraction)
                htau_dn.Scale(1.0 - varfraction)

                hsig_up.Add(htau_up)
                hsig_dn.Add(htau_dn)
                SaveHistToFile(hsig_up, foutput)
                SaveHistToFile(hsig_dn, foutput)

    foutput.Close()


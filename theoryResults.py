from collections import OrderedDict


class XSecResults(object):
    def __init__(self):
        self.xsec_wp_13 = (0.,0., 0.)
        self.xsec_wm_13 = (0.,0.,0.)
        self.xsec_z0_13 = (0.,0.,0.)
        self.xsec_winc_13 = (0.,0.,0.)
        self.xsec_wp_5 = (0.,0.,0.)
        self.xsec_wm_5 = (0.,0.,0.)
        self.xsec_z0_5 = (0.,0.,0.)
        self.xsec_winc_5 = (0.,0.,0.)
        
        self.xsec_WchgRatio_13 = (0.,0.,0.)
        self.xsec_WZRatio_13 = (0.,0.,0.)
        self.xsec_WchgRatio_5 = (0.,0.,0.)
        self.xsec_WZRatio_5 = (0.,0.,0.)
        
        self.xsec_wp_13o5 = (0.,0.,0.)
        self.xsec_wm_13o5 = (0.,0.,0.)
        self.xsec_z0_13o5 = (0.,0.,0.)
        self.xsec_winc_13o5 = (0.,0.,0.)
        
        self.xsec_WchgRatio_13o5 = (0.,0.,0.)
        self.xsec_WZRatio_13o5 = (0.,0.,0.)
        
    def GetXsec(self, ch, sqrtS):
        era = "13" if sqrtS == "13TeV" else "5"
        channelMaps = {"lepplus": "wp", "lepminus": "wm", "leplep": "z0", "lepinc": "winc"}
        assert ch in ["wp", "wm", "winc","z0", "lepplus", "lepminus", "winc", "leplep"], "channel must be wp, wm, winc, z0, or lepplus, lepminus, winc, leplep"
        if ch in channelMaps:
            ch = channelMaps[ch]
        return getattr(self, f"xsec_{ch}_{era}")
        
    def GetWchgRatio(self, is13TeV=True):
        if is13TeV:
            if self.xsec_WchgRatio_13[0] == 0.:
                return (self.xsec_wp_13[0] / self.xsec_wm_13[0], 0., 0.)
            return self.xsec_WchgRatio_13
        else:
            if self.xsec_WchgRatio_5[0] == 0.:
                return (self.xsec_wp_5[0] / self.xsec_wm_5[0], 0., 0.)
            return self.xsec_WchgRatio_5
    
    def GetWZRatio(self, is13TeV=True):
        if is13TeV:
            if self.xsec_WZRatio_13[0] == 0.:
                return ((self.xsec_wp_13[0] + self.xsec_wm_13[0])/ self.xsec_z0_13[0], 0., 0.)
            return self.xsec_WZRatio_13
        else:
            if self.xsec_WZRatio_5[0] == 0.:
                return ((self.xsec_wp_5[0] + self.xsec_wm_5[0])/ self.xsec_z0_5[0], 0., 0.)
            return self.xsec_WZRatio_5
    
    def GetSqrtSRatio(self, channel = "wp"):
        channelMaps = {"lepplus": "wp", "lepminus": "wm", "leplep": "z0", "lepinc": "winc"}
        if channel in channelMaps:
            channel = channelMaps[channel]
        assert channel in ["wp", "wm", "winc","z0"], "channel must be wp, wm, winc, z0, or lepplus, lepminus, winc, leplep"
        #print("sqrts ratio ", getattr(self, f"xsec_{channel}_13o5"))
        if getattr(self, f"xsec_{channel}_13o5")[0] == 0:
            num = getattr(self, f"xsec_{channel}_13")
            den = getattr(self, f"xsec_{channel}_5")
            return (num[0] / den[0], 0., 0.)
        return getattr(self, f"xsec_{channel}_13o5")
    
    def GetSqrtSDoubleRatio(self, channel = "Wchg"):
        channelMaps = {"WchgRatio": "Wchg", "WZRatio": "WZ"}
        if channel in channelMaps:
            channel = channelMaps[channel]
        assert channel in ["Wchg", "WZ"], "channel must be Wchg, WZ, or WchgRatio, WZRatio"
        if getattr(self, f"xsec_{channel}Ratio_13o5")[0] == 0:
            num = getattr(self, f"Get{channel}Ratio")(is13TeV=True)
            den = getattr(self, f"Get{channel}Ratio")(is13TeV=False)
            return (num[0] / den[0], 0., 0.)
        return getattr(self, f"xsec_{channel}Ratio_13o5")
        
    
xsec_fid = OrderedDict()

xsec_fid['nnpdf4.0'] = XSecResults()
xsec_fid['nnpdf3.1'] = XSecResults()
xsec_fid['ct18'] = XSecResults()
xsec_fid['msht20'] = XSecResults()

#
# toDo: change these hard-coded values into a function, 
# with uncertainties as well
#

xsec_fid['nnpdf4.0'].xsec_wp_13 =  (5135.0,0.,0.)
xsec_fid['nnpdf4.0'].xsec_wm_13 =  (3913.0,0.,0.)
xsec_fid['nnpdf4.0'].xsec_z0_13 =  ( 750.0,0.,0.)
xsec_fid['nnpdf4.0'].xsec_wp_5  =  (2496.0,0.,0.)
xsec_fid['nnpdf4.0'].xsec_wm_5  =  (1537.4,0.,0.)
xsec_fid['nnpdf4.0'].xsec_z0_5  =  ( 324.0,0.,0.)

xsec_fid['nnpdf3.1'].xsec_wp_13 = (5075.0,0.,0.)
xsec_fid['nnpdf3.1'].xsec_wm_13 = (3847.0,0.,0.)
xsec_fid['nnpdf3.1'].xsec_z0_13 =  (734.3,0.,0.)
xsec_fid['nnpdf3.1'].xsec_wp_5 =  (2468.2,0.,0.)
xsec_fid['nnpdf3.1'].xsec_wm_5 =  (1521.7,0.,0.)
xsec_fid['nnpdf3.1'].xsec_z0_5 =   (319.0,0.,0.)

xsec_fid['ct18'].xsec_wp_13 = (  4977.0,0,0.)
xsec_fid['ct18'].xsec_wm_13 = (  3825.0,0,0.)
xsec_fid['ct18'].xsec_z0_13 = (   717.1,0,0.)
xsec_fid['ct18'].xsec_wp_5  = ( 2419.0,0,0.)
xsec_fid['ct18'].xsec_wm_5  = ( 1509.0,0,0.)
xsec_fid['ct18'].xsec_z0_5  = (  309.0,0,0.)

xsec_fid['msht20'].xsec_wp_13 = (  4967.0,0,0.)
xsec_fid['msht20'].xsec_wm_13 = (  3840.0,0,0.)
xsec_fid['msht20'].xsec_z0_13 = (   736.6,0,0.)
xsec_fid['msht20'].xsec_wp_5  = (  2417.5,0,0.)
xsec_fid['msht20'].xsec_wm_5  = (  1492.4,0,0.)
xsec_fid['msht20'].xsec_z0_5  = (   313.2,0,0.)

#
# inc
#
xsec_inc = OrderedDict()

xsec_inc['nnpdf4.0'] = XSecResults()
xsec_inc['nnpdf3.1'] = XSecResults()
xsec_inc['ct18'] = XSecResults()
xsec_inc['msht20'] = XSecResults()

xsec_inc['nnpdf4.0'].xsec_wp_13 = ( 11625.6,0, 0.)
xsec_inc['nnpdf4.0'].xsec_wm_13 = ( 8602.6, 0, 0.)
xsec_inc['nnpdf4.0'].xsec_z0_13 = (  1963.5,0, 0.)
xsec_inc['nnpdf4.0'].xsec_wp_5 = ( 4432.8,  0, 0.)
xsec_inc['nnpdf4.0'].xsec_wm_5 = ( 2906.7,  0, 0.)
xsec_inc['nnpdf4.0'].xsec_z0_5 = (  683.0,  0, 0.)

xsec_inc['nnpdf3.1'].xsec_wp_13 = (11500.4,0,0.)
xsec_inc['nnpdf3.1'].xsec_wm_13 = (8492.2, 0,0.)
xsec_inc['nnpdf3.1'].xsec_z0_13 = ( 1933.8,0,0.)
xsec_inc['nnpdf3.1'].xsec_wp_5 = ( 4381.2, 0,0.)
xsec_inc['nnpdf3.1'].xsec_wm_5 = ( 2872.2, 0,0.)
xsec_inc['nnpdf3.1'].xsec_z0_5 = (  673.1, 0,0.)

xsec_inc['ct18'].xsec_wp_13 = ( 11489.8,0,0.)
xsec_inc['ct18'].xsec_wm_13 = ( 8462.8, 0,0.)
xsec_inc['ct18'].xsec_z0_13 = ( 1914.3, 0,0.)
xsec_inc['ct18'].xsec_wp_5 = (4308.8,   0,0.)
xsec_inc['ct18'].xsec_wm_5 = (2844.1,   0,0.)
xsec_inc['ct18'].xsec_z0_5 = ( 658.6,   0,0.)

xsec_inc['msht20'].xsec_wp_13 = ( 11389.0,0,0.)
xsec_inc['msht20'].xsec_wm_13 = ( 8407.9, 0,0.)
xsec_inc['msht20'].xsec_z0_13 = ( 1928.3, 0,0.)
xsec_inc['msht20'].xsec_wp_5 = ( 4284.7,  0,0.)
xsec_inc['msht20'].xsec_wm_5 = ( 2818.3,  0,0.)
xsec_inc['msht20'].xsec_z0_5 = (  661.1,  0,0.)

import json
# load json file
pdfs = ["NNPDF40_nnlo_hessian_pdfas", "NNPDF31_nnlo_hessian_pdfas", "CT18NNLO", "MSHT20nnlo_as118"]
pdfMaps = {
    "NNPDF40_nnlo_hessian_pdfas": "nnpdf4.0",
    "NNPDF31_nnlo_hessian_pdfas": "nnpdf3.1",
    "CT18NNLO": "ct18",
    "MSHT20nnlo_as118": "msht20"
}
for poi in ["fid", "inc"]:
    xsec_results = xsec_fid if poi == "fid" else xsec_inc
    suffix = "_fid" if poi == "fid" else "_inc"
    for pdf in pdfs:
        with open(f"data/results_{pdf}{suffix}.json") as json_file:
           results = json.load(json_file) 
        for sqrtS in ["13", "5"]:
            if sqrtS+"tev" in results[0]:
                for key, val in results[0][sqrtS+"tev"].items():
                    #print("key", key, "val", val, "poi", poi)
                    if poi == "fid":
                        setattr(xsec_results[pdfMaps[pdf]], f"xsec_{key}_{sqrtS}", (val[0] / 1e3, 0., 0.))
                    else:
                        setattr(xsec_results[pdfMaps[pdf]], f"xsec_{key}_{sqrtS}", (val[0] / 1e3, val[1] / 1e3, val[2] / 1e3))
            if sqrtS+"tev" in results[1]:
                for key, val in results[1][sqrtS+"tev"].items():
                    #print("key", key, "val", val)
                    if poi == "fid":
                        setattr(xsec_results[pdfMaps[pdf]], f"xsec_{key}_{sqrtS}", (val[0], 0., 0.))
                    else:
                        setattr(xsec_results[pdfMaps[pdf]], f"xsec_{key}_{sqrtS}", (val[0], val[1], val[2]))
    
        # ratio
        #print("results 2", results[2])
        for key, val in results[2].items():
            #print("key ", ch, " val", val)
            if poi == "fid":
                setattr(xsec_results[pdfMaps[pdf]], f"xsec_{key}_13o5", (val[0], 0., 0.))
            else:
                setattr(xsec_results[pdfMaps[pdf]], f"xsec_{key}_13o5", (val[0], val[1], val[2]))
    
        # double ratio
        #print("results 3", results[3])
        for key, val in results[3].items():
            #print("key ", key, " val ", val)
            if poi == "fid":
                setattr(xsec_results[pdfMaps[pdf]], f"xsec_{key}_13o5", (val[0], 0., 0.))
            else:
                setattr(xsec_results[pdfMaps[pdf]], f"xsec_{key}_13o5", (val[0], val[1], val[2]))

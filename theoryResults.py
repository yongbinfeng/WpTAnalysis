from collections import OrderedDict


class XSecResults(object):
    def __init__(self):
        self.xsec_wp_13 = 0.
        self.xsec_wm_13 = 0.
        self.xsec_z0_13 = 0.
        self.xsec_wp_5 = 0.
        self.xsec_wm_5 = 0.
        self.xsec_z0_5 = 0.
        
    def GetXsec(self, ch, sqrtS):
        era = "13" if sqrtS == "13TeV" else "5"
        ch = "wp" if ch == "lepplus" else "wm" if ch == "lepminus" else "z0"
        return getattr(self, f"xsec_{ch}_{era}")
        
    def GetWchgRatio(self, is13TeV=True):
        if is13TeV:
            return self.xsec_wp_13 / self.xsec_wm_13
        else:
            return self.xsec_wp_5 / self.xsec_wm_5
    
    def GetWZRatio(self, is13TeV=True):
        if is13TeV:
            w_tot = self.xsec_wp_13 + self.xsec_wm_13
            z_tot = self.xsec_z0_13
        else:
            w_tot = self.xsec_wp_5 + self.xsec_wm_5
            z_tot = self.xsec_z0_5
        return w_tot / z_tot
    
    def GetSqrtSRatio(self, channel = "wp"):
        channelMaps = {"lepplus": "wp", "lepminus": "wm", "leplep": "z0"}
        if channel in channelMaps:
            channel = channelMaps[channel]
        assert channel in ["wp", "wm", "z0"], "channel must be wp, wm, z0, or lepplus, lepminus, leplep"
        num = getattr(self, f"xsec_{channel}_13")
        den = getattr(self, f"xsec_{channel}_5")
        return num / den
    
    def GetSqrtSDoubleRatio(self, channel = "Wchg"):
        channelMaps = {"WchgRatio": "Wchg", "WZRatio": "WZ"}
        if channel in channelMaps:
            channel = channelMaps[channel]
        assert channel in ["Wchg", "WZ"], "channel must be Wchg, WZ, or WchgRatio, WZRatio"
        num = getattr(self, f"Get{channel}Ratio")(is13TeV=True)
        den = getattr(self, f"Get{channel}Ratio")(is13TeV=False)
        return num / den
        
    
xsec_fid = OrderedDict()

xsec_fid['nnpdf4.0'] = XSecResults()
xsec_fid['nnpdf3.1'] = XSecResults()
xsec_fid['ct18'] = XSecResults()
xsec_fid['msht20'] = XSecResults()

#
# toDo: change these hard-coded values into a function, 
# with uncertainties as well
#

xsec_fid['nnpdf4.0'].xsec_wp_13 =  5135.0
xsec_fid['nnpdf4.0'].xsec_wm_13 =  3913.0
xsec_fid['nnpdf4.0'].xsec_z0_13 =   750.0
xsec_fid['nnpdf4.0'].xsec_wp_5 =  2496.0
xsec_fid['nnpdf4.0'].xsec_wm_5 =  1537.4
xsec_fid['nnpdf4.0'].xsec_z0_5 =   324.0

xsec_fid['nnpdf3.1'].xsec_wp_13 = 5075.0
xsec_fid['nnpdf3.1'].xsec_wm_13 = 3847.0
xsec_fid['nnpdf3.1'].xsec_z0_13 =  734.3
xsec_fid['nnpdf3.1'].xsec_wp_5 =  2468.2
xsec_fid['nnpdf3.1'].xsec_wm_5 =  1521.7
xsec_fid['nnpdf3.1'].xsec_z0_5 =   319.0

xsec_fid['ct18'].xsec_wp_13 =  4977.0
xsec_fid['ct18'].xsec_wm_13 =  3825.0
xsec_fid['ct18'].xsec_z0_13 =   717.1
xsec_fid['ct18'].xsec_wp_5 = 2419.0
xsec_fid['ct18'].xsec_wm_5 = 1509.0
xsec_fid['ct18'].xsec_z0_5 =  309.0

xsec_fid['msht20'].xsec_wp_13 =  4967.0
xsec_fid['msht20'].xsec_wm_13 =  3840.0
xsec_fid['msht20'].xsec_z0_13 =   736.6
xsec_fid['msht20'].xsec_wp_5 =  2417.5
xsec_fid['msht20'].xsec_wm_5 =  1492.4
xsec_fid['msht20'].xsec_z0_5 =   313.2

#
# inc
#
xsec_inc = OrderedDict()

xsec_inc['nnpdf4.0'] = XSecResults()
xsec_inc['nnpdf3.1'] = XSecResults()
xsec_inc['ct18'] = XSecResults()
xsec_inc['msht20'] = XSecResults()

xsec_inc['nnpdf4.0'].xsec_wp_13 =  11625.6
xsec_inc['nnpdf4.0'].xsec_wm_13 =  8602.6
xsec_inc['nnpdf4.0'].xsec_z0_13 =   1963.5
xsec_inc['nnpdf4.0'].xsec_wp_5 =  4432.8
xsec_inc['nnpdf4.0'].xsec_wm_5 =  2906.7
xsec_inc['nnpdf4.0'].xsec_z0_5 =   683.0

xsec_inc['nnpdf3.1'].xsec_wp_13 = 11500.4
xsec_inc['nnpdf3.1'].xsec_wm_13 = 8492.2
xsec_inc['nnpdf3.1'].xsec_z0_13 =  1933.8
xsec_inc['nnpdf3.1'].xsec_wp_5 =  4381.2
xsec_inc['nnpdf3.1'].xsec_wm_5 =  2872.2
xsec_inc['nnpdf3.1'].xsec_z0_5 =   673.1

xsec_inc['ct18'].xsec_wp_13 =  11489.8
xsec_inc['ct18'].xsec_wm_13 =  8462.8
xsec_inc['ct18'].xsec_z0_13 =  1914.3
xsec_inc['ct18'].xsec_wp_5 = 4308.8
xsec_inc['ct18'].xsec_wm_5 = 2844.1
xsec_inc['ct18'].xsec_z0_5 =  658.6

xsec_inc['msht20'].xsec_wp_13 =  11389.0
xsec_inc['msht20'].xsec_wm_13 =  8407.9
xsec_inc['msht20'].xsec_z0_13 =  1928.3
xsec_inc['msht20'].xsec_wp_5 =  4284.7
xsec_inc['msht20'].xsec_wm_5 =  2818.3
xsec_inc['msht20'].xsec_z0_5 =   661.1

import math
from enum import Enum

class Process(Enum):
    wp = 0
    wm = 1
    winc = 2
    z = 3

class XsecValues(object):
    def __init__(self, sqrts, xsec = 0., stat = 0., syst = 0., theo = 0., lumi = 0., totL = 0., proc = Process.wp):
        self.xsec = xsec
        self.stat = stat
        self.syst = syst
        self.theo = theo
        self.lumi = lumi
        
        self.sqrts = sqrts
        self.totL = totL
        
        if isinstance(proc, str):
            proc = Process[proc]
        self.proc = proc   
        
    def SetXsec(self, xsec, stat = 0., syst = 0., theo = 0., lumi = 0.):
        self.xsec = xsec
        self.stat = stat
        self.syst = syst
        self.theo = theo
        self.lumi = lumi
        
    def GetUnc(self):
        return math.sqrt(self.stat**2 + self.syst**2 + self.theo**2 + self.lumi**2)
    
    def GetXsec(self):
        return self.xsec, self.GetUnc()

    def __str__(self):
        return f"xsec: {self.xsec} +/- {self.GetUnc()} at sqrtS = {self.sqrts} with L = {self.totL}"
    
    
class XsecMeasurements(object):
    def __init__(self, sqrts = 0., totL = "0. pb^{{-1}}"):
        self.sqrts = sqrts
        self.totL = totL
        self.xsecs = [
            XsecValues(sqrts, totL = totL, proc = Process.wp), 
            XsecValues(sqrts, totL = totL, proc = Process.wm), 
            XsecValues(sqrts, totL = totL, proc = Process.winc), 
            XsecValues(sqrts, totL = totL, proc = Process.z)
            ]
    
    def __getitem__(self, key):
        if isinstance(key, str):
            key = Process[key].value
        return self.xsecs[key]
    
    def __setitem__(self, key, value):
        if isinstance(key, str):
            key = Process[key].value
        self.xsecs[key] = value
        
    def __len__(self):
        return len(self.xsecs)
    
    def GetLegend(self):
        return f"{self.totL} ({self.sqrts} TeV)"
    
# CMS 2.76 TeV results
xsecs_CMS_2760GeV = XsecMeasurements(2.76, "5.4 pb^{-1}")
xsecs_CMS_2760GeV['z'].SetXsec(298.0, 11.0)

xsecs_CMS_2760GeV_2 = XsecMeasurements(2.76, "7.3 #mub^{-1}")
xsecs_CMS_2760GeV_2['wp']   .SetXsec(  2380.0,  140.0, 180., 0., 0.)
xsecs_CMS_2760GeV_2['wm']   .SetXsec(  1450.0,  110.0, 110., 0., 0.)
xsecs_CMS_2760GeV_2['winc'] .SetXsec(  3830.0,  180.0, 290., 0., 0.)
    
# CMS 8 TeV result
# from https://arxiv.org/pdf/1402.0923.pdf
# Table 4
xsecs_CMS_8TeV = XsecMeasurements(8, "18.2 pb^{-1}")
xsecs_CMS_8TeV['wp']   .SetXsec(7110,  30, 140, 0., 180.)
xsecs_CMS_8TeV['wm']   .SetXsec(5090,  20, 110, 0., 130.)
xsecs_CMS_8TeV['winc'] .SetXsec(12210, 30, 240, 0., 320.)
xsecs_CMS_8TeV['z']    .SetXsec(1150,  10, 20,  0., 30.)

# CMS 7 TeV result
# from https://arxiv.org/pdf/1107.4789.pdf
xsecs_CMS_7TeV = XsecMeasurements(7, "36 pb^{-1}")
xsecs_CMS_7TeV['wp']   .SetXsec(6040,  20, 60,  80,  240.)
xsecs_CMS_7TeV['wm']   .SetXsec(4260,  10, 40,  70,  170.)
xsecs_CMS_7TeV['winc'] .SetXsec(10300, 20, 100, 100, 410.)
xsecs_CMS_7TeV['z']    .SetXsec(974,    7,  7,  18,  39)

# ATLAS 7 TeV result
# from https://arxiv.org/pdf/1612.03016.pdf
# table 9
xsecs_ATLAS_7TeV = XsecMeasurements(7, "4.6 fb^{-1}")
xsecs_ATLAS_7TeV['wp']   .SetXsec(6350,  2, 30,  100,  110.)
xsecs_ATLAS_7TeV['wm']   .SetXsec(4376,  2, 25,  90,   79.)
xsecs_ATLAS_7TeV['winc'] .SetXsec(10720, 3, 60,  130,  190.)
xsecs_ATLAS_7TeV['z']    .SetXsec(990,   1, 3,   15,   18.)

# CDF Run2
# from https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.94.091803
xsecs_CDF_Run2 = XsecMeasurements(1.96, "72.0 pb^{-1}")
xsecs_CDF_Run2['winc'].SetXsec(2775,  10,  53,  0.,  167)
xsecs_CDF_Run2['z']   .SetXsec(254.9, 3.3, 4.6, 0.,  15.2)

# D0 Run1
# from https://journals.aps.org/prd/pdf/10.1103/PhysRevD.61.072001
xsecs_D0_Run1 = XsecMeasurements(1.80, "84.5 pb^{-1}")
xsecs_D0_Run1['winc'].SetXsec(2310,  10,  50,  0.,  100)
xsecs_D0_Run1['z']   .SetXsec(221,    3,   4,  0.,   10)

# UA2 
# from https://link.springer.com/article/10.1007/BF01551906
xsecs_UA2 = XsecMeasurements(0.63, "7.8 pb^{-1}")
xsecs_UA2['winc'].SetXsec(660,   15, 37,  0., 0.)
xsecs_UA2['z']   .SetXsec(70.4, 5.5, 4.0, 0., 0.)

# UA1
# from https://www.sciencedirect.com/science/article/pii/0370269387915103?via%3Dihub
xsecs_UA1 = XsecMeasurements(0.63, "0.4 pb^{-1}")
xsecs_UA1['winc'].SetXsec(630,   40, 100,  0., 0.)
xsecs_UA1['z']   .SetXsec(71.,   11,  11,  0., 0.)

xsecs_CMS_13TeV = XsecMeasurements(13, "206 pb^{-1}")
xsecs_CMS_13TeV['wp']   .SetXsec( 11800.0,  10.0,  100.0,  0.0,  270.0)
xsecs_CMS_13TeV['wm']   .SetXsec(  8670.0,  10.0,   80.0,  0.0,  200.0)
xsecs_CMS_13TeV['winc'] .SetXsec( 20480.0,  10.0,  170.0,  0.0,  470.0)
xsecs_CMS_13TeV['z']    .SetXsec( 1952.0,    4.0,   18.0,  0.0,  45.0)

xsecs_CMS_5TeV = XsecMeasurements(5.02, "298 pb^{-1}")
xsecs_CMS_5TeV['wp']   .SetXsec(  4401.0,  4.0,  36.0,  0.0,  84.0)
xsecs_CMS_5TeV['wm']   .SetXsec(  2897.0,  4.0,  26.0,  0.0,  55.0)
xsecs_CMS_5TeV['winc'] .SetXsec(  7300.0,  10.0, 60.0,  0.0,  140.0)
xsecs_CMS_5TeV['z']    .SetXsec(   669.0,  2.0,  6.0,  0.0,  13.0)

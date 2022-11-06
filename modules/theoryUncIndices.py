
class TheoryUncMap(object):
    def __init__(self, MuRUp, MuRDown, MuFUp, MuFDown, MuFMuRUp, MuFMuRDown, PDFStart, PDFEnd, alphaSUp, alphaSDown):
        self.MuRUp = MuRUp
        self.MuRDown = MuRDown
        self.MuFUp = MuFUp
        self.MuFDown = MuFDown
        self.MuFMuRUp = MuFMuRUp
        self.MuFMuRDown = MuFMuRDown
        self.PDFStart = PDFStart
        self.PDFEnd = PDFEnd
        self.alphaSUp = alphaSUp
        self.alphaSDown = alphaSDown
    
    def getIndex(self, var):
        return getattr(self, var)

theoryUnc_13TeV_wjets = TheoryUncMap(1, 2, 3, 6, 4, 8, 10, 109, 110, None)
theoryUnc_5TeV_wjets  = TheoryUncMap(0, 1, 2, 5, 3, 7,  9, 108, 109, 110)

theoryUnc_13TeV_zjets = TheoryUncMap(0, 1, 2, 5, 3, 7,  9, 108, 109, 110)
theoryUnc_5TeV_zjets  = TheoryUncMap(0, 1, 2, 5, 3, 7,  9, 108, 109, 110)
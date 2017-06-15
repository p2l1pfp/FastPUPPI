import ROOT, copy, string, os

class AbsCollectionDrawer(object):
    def __init__(self,label="",filter=None):
        self.filter = filter
        self.label = label
    def drawOne(self,index,object):
        pass
    def draw(self,objs):
        for i,o in enumerate(objs):
            if self.filter and not self.filter(o): 
                continue
            self.drawOne(i+1,o)
        return True
    def tobjForLegend(self):
        return None
    def clone(self, name):
        raise RuntimeError, "Not implemented"

class MultiGraphMarkerDrawerVsPt(AbsCollectionDrawer):
    def __init__(self,markerStyle=20,markerColor=1,markerSizes=[2],ptEdges=[0,9e9],maxItems=99999,filter=None,label="",getcoord=lambda o : (o.eta(),o.phi())):
        AbsCollectionDrawer.__init__(self,filter=filter,label=label)
        self._markers = [ ]
        self._getcoord = getcoord
        for i,markerSize in enumerate(markerSizes):
            marker = ROOT.TGraph()
            marker.SetMarkerStyle(markerStyle)
            marker.SetMarkerColor(markerColor)
            marker.SetMarkerSize(markerSize)
            marker.ptEdges = (ptEdges[i],ptEdges[i+1])
            self._markers.append(marker)
        self._maxItems = maxItems
    def draw(self,objs):
        todraw = [[] for m in self._markers]
        for o in objs:
            if self.filter and not self.filter(o): 
                continue
            for im,m in enumerate(self._markers):
                if o.pt() < m.ptEdges[0]: break
                if o.pt() < m.ptEdges[1]:
                    todraw[im].append(o)
                    break
        for im,m in enumerate(self._markers):
            if len(todraw[im]) == 0: continue
            m.Set(len(todraw[im]))
            for i,o in enumerate(todraw[im]):
                eta, phi = self._getcoord(o)
                m.SetPoint(i, eta, phi) 
            m.Draw("P")
        return True
    def tobjForLegend(self):
        return self._markers[-1];
    def clone(self, name):
        ret = copy.copy(self)
        ret._markers = [ g.Clone() for g in self._markers ]
        for src,dst in zip(self._markers, ret._markers): dst.ptEdges = src.ptEdges[:]
        return ret

class CircleDrawer(AbsCollectionDrawer):
    def __init__(self,radius,lineColor=1,lineWidth=2,lineStyle=1,maxItems=99,filter=None):
        AbsCollectionDrawer.__init__(self,filter)
        self._radius = radius
        self._circle = ROOT.TEllipse(0,0,radius)
        self._circle.SetLineColor(lineColor)
        self._circle.SetLineWidth(lineWidth)
        self._circle.SetLineStyle(lineStyle)
        self._circle.SetFillStyle(0)
        self._maxItems = maxItems
    def drawOne(self,index,object):
        if (index >= self._maxItems): return
        self._circle.DrawEllipse(object.eta(),object.phi(),self._radius,self._radius,0.,360.,0.)
    def clone(self, name):
        return copy.copy(self)

class ExtrapToCaloDrawer(AbsCollectionDrawer):
    def __init__(self,lineColor=1,lineWidth=2,lineStyle=1,maxItems=99,filter=None):
        AbsCollectionDrawer.__init__(self,filter)
        self._line = ROOT.TLine(0,0,1,1)
        self._line.SetLineColor(lineColor)
        self._line.SetLineWidth(lineWidth)
        self._line.SetLineStyle(lineStyle)
        self._maxItems = maxItems
    def drawOne(self,index,object):
        if (index >= self._maxItems): return
        self._line.DrawLine(object.eta(),object.phi(),object.caloEta(),object.caloPhi())
    def clone(self, name):
        return copy.copy(self)

class CanvasMaker(object):
    def __init__(self,maxEvents=100,maxEta=5,vlines=[]):
        ROOT.gROOT.ProcessLine(".x %s/src/FastPUPPI/NtupleProducer/python/display/tdrstyle.cc" % os.environ['CMSSW_BASE']);
        ROOT.gStyle.SetPaperSize(25.,25*(3.14/maxEta))
        self._maxEta = maxEta
        self._canvas = ROOT.TCanvas("c1","c1",int((2*maxEta)*1.7*100), int(6.28*1.7*100));
        self._canvas.SetWindowSize(int((2*maxEta)*1.7*100), int(6.28*1.7*100))
        self._canvas.SetGridx(1)
        self._canvas.SetGridy(1)
        self._frame = ROOT.TH1F("frame","frame;#eta;#phi",1000,-maxEta,maxEta);
        self._frame.GetYaxis().SetRangeUser(-3.1416,3.1416);
        self._frame.SetStats(0);
        self._maxEvents = maxEvents
        self._currEvents = 0
        self._vlines = vlines 
        self._tvline = ROOT.TLine(0,0,1,1)
        self._tvline.SetLineWidth(3)
        self._tvline.SetLineStyle(7)
    def canvas(self): 
        return self._canvas
    def go(self,event):
        self._currEvents += 1;
        if self._currEvents > self._maxEvents: return False
        self.clear();
        return True
    def clear(self):
        self._canvas.Clear()
        self._frame.GetYaxis().SetRangeUser(-3.1416,3.1416);
        self._frame.GetXaxis().SetRangeUser(-self._maxEta,self._maxEta);
        self._frame.Draw()
        for x in self._vlines:
            self._tvline.DrawLine(x, -3.1416, x, +3.1416)
    def zoom(self,center,radius):
        xmin, xmax = max(center[0]-radius,-self._maxEta), min(center[0]+radius,self._maxEta)
        ymin, ymax = max(center[1]-radius,-3.1416), min(center[1]+radius,+3.1416)
        self._frame.Draw()
        self._frame.GetXaxis().SetRangeUser(xmin,xmax);
        self._frame.GetYaxis().SetRangeUser(ymin,ymax);
        for x in self._vlines:
            if (xmin < x and x < xmax): self._tvline.DrawLine(x, ymin, x, ymax)
    def done(self,event):
        for template in self._fileTemplate:
            self.write(event,template) 
        return True
    def write(self,event,template):
        fname = string.Formatter().vformat(template,[],{'run':event.eventAuxiliary().run(), 'lumi':event.eventAuxiliary().luminosityBlock(), 'evt':event.eventAuxiliary().event()})
        self._canvas.Print(fname)
        print fname


cMaker  = CanvasMaker(maxEta=5, vlines=[-4,-3,-1.5,1.5,3,4])
drawGenJets = [ CircleDrawer(0.4, lineColor=ROOT.kBlue+2,  lineStyle=1, lineWidth=3) ]

ptEdges = [1.0,2.5,7,9e9]
drawGenCands = [
    MultiGraphMarkerDrawerVsPt( filter = lambda d : d.charge() != 0,                          markerStyle=34, markerSizes=[0.8,1.0,1.2], ptEdges=ptEdges, markerColor=ROOT.kRed+2 , label="gen h^{#pm}"),
    MultiGraphMarkerDrawerVsPt( filter = lambda d : d.charge() == 0 and abs(d.pdgId()) == 22, markerStyle=34, markerSizes=[0.8,1.0,1.2], ptEdges=ptEdges, markerColor=ROOT.kGreen+3 , label="gen #gamma" ),
    MultiGraphMarkerDrawerVsPt( filter = lambda d : d.charge() == 0 and abs(d.pdgId()) != 22, markerStyle=34, markerSizes=[0.8,1.0,1.2], ptEdges=ptEdges, markerColor=ROOT.kBlue+2 , label="gen h^{0}" ),
]
drawTracks   = [ MultiGraphMarkerDrawerVsPt( markerStyle=20, markerSizes=[0.9,1.3,1.6], ptEdges=[2,3,10,9e9], markerColor=ROOT.kMagenta+0, label="track" ), ]
drawTracksAC = [ MultiGraphMarkerDrawerVsPt( markerStyle=20, markerSizes=[0.8,1.2,1.5], ptEdges=[2,3,10,9e9], markerColor=ROOT.kViolet+1, label="tk @ cal", getcoord = lambda o : (o.caloEta(),o.caloPhi())),
                 ExtrapToCaloDrawer(lineColor = ROOT.kViolet+1, lineWidth=1) ]
drawHcal     = [ MultiGraphMarkerDrawerVsPt( markerStyle=21, markerSizes=[0.8,1.2,1.5], ptEdges=ptEdges, markerColor=ROOT.kAzure+10, label="had cal" ), ]
drawEcal     = [ MultiGraphMarkerDrawerVsPt( markerStyle=21, markerSizes=[0.8,1.2,1.5], ptEdges=ptEdges, markerColor=ROOT.kGreen+0, label="em cal" ), ]
drawHGCalE   = [ MultiGraphMarkerDrawerVsPt( markerStyle=21, markerSizes=[0.8,1.2,1.5], ptEdges=ptEdges, markerColor=ROOT.kGreen+0 ), ]
drawHGCalH   = [ MultiGraphMarkerDrawerVsPt( markerStyle=21, markerSizes=[0.8,1.2,1.5], ptEdges=ptEdges, markerColor=ROOT.kAzure+10 ), ]
drawHGCalB   = [ MultiGraphMarkerDrawerVsPt( markerStyle=21, markerSizes=[0.8,1.2,1.5], ptEdges=ptEdges, markerColor=ROOT.kAzure+1 ), ]

drawEmCaloCands = [
    MultiGraphMarkerDrawerVsPt( markerStyle=25, markerSizes=[1.2,1.5,1.8], ptEdges=ptEdges, markerColor=ROOT.kGreen+2 ),
]
drawCaloCands = [
    MultiGraphMarkerDrawerVsPt( filter = lambda d : d.charge() == 0 and abs(d.pdgId()) == 22, markerStyle=25, markerSizes=[1.2,1.5,1.8], ptEdges=ptEdges, markerColor=ROOT.kGreen+0 ),
    MultiGraphMarkerDrawerVsPt( filter = lambda d : d.charge() == 0 and abs(d.pdgId()) != 22, markerStyle=25, markerSizes=[1.2,1.5,1.8], ptEdges=ptEdges, markerColor=ROOT.kBlue+0 ),
]
drawTkCands   = [ MultiGraphMarkerDrawerVsPt( markerStyle=25, markerSizes=[1.3,1.6,1.9], ptEdges=[2,3,10,9e9], markerColor=ROOT.kViolet+1 ), ]
drawPFCands = [
    MultiGraphMarkerDrawerVsPt( filter = lambda d : d.charge() != 0,                          markerStyle=24, markerSizes=[1.4,1.8,2.3], ptEdges=ptEdges, markerColor=ROOT.kRed+1, label="rec h^{#pm}" ),
    MultiGraphMarkerDrawerVsPt( filter = lambda d : d.charge() != 0,                          markerStyle=24, markerSizes=[1.3,1.7,2.2], ptEdges=ptEdges, markerColor=ROOT.kRed+1 ),
    MultiGraphMarkerDrawerVsPt( filter = lambda d : d.charge() == 0 and abs(d.pdgId()) == 22, markerStyle=24, markerSizes=[1.4,1.8,2.3], ptEdges=ptEdges, markerColor=ROOT.kGreen+2 , label="rec #gamma" ),
    MultiGraphMarkerDrawerVsPt( filter = lambda d : d.charge() == 0 and abs(d.pdgId()) == 22, markerStyle=24, markerSizes=[1.3,1.7,2.2], ptEdges=ptEdges, markerColor=ROOT.kGreen+2 ),
    MultiGraphMarkerDrawerVsPt( filter = lambda d : d.charge() == 0 and abs(d.pdgId()) != 22, markerStyle=24, markerSizes=[1.4,1.8,2.3], ptEdges=ptEdges, markerColor=ROOT.kBlue+1 , label="rec h^{0}" ),
    MultiGraphMarkerDrawerVsPt( filter = lambda d : d.charge() == 0 and abs(d.pdgId()) != 22, markerStyle=24, markerSizes=[1.3,1.7,2.2], ptEdges=ptEdges, markerColor=ROOT.kBlue+1 ),
]
drawPuppiCands = [
    MultiGraphMarkerDrawerVsPt( filter = lambda d : d.charge() != 0,                          markerStyle=20, markerSizes=[1.1,1.4,1.7], ptEdges=ptEdges, markerColor=ROOT.kRed-9 ),
    MultiGraphMarkerDrawerVsPt( filter = lambda d : d.charge() == 0 and abs(d.pdgId()) == 22, markerStyle=20, markerSizes=[1.1,1.4,1.7], ptEdges=ptEdges, markerColor=ROOT.kGreen+0 ),
    MultiGraphMarkerDrawerVsPt( filter = lambda d : d.charge() == 0 and abs(d.pdgId()) != 22, markerStyle=20, markerSizes=[1.1,1.4,1.7], ptEdges=ptEdges, markerColor=ROOT.kBlue+0 ),
]

drawPUCands = [
    MultiGraphMarkerDrawerVsPt( filter = lambda d : d.charge() != 0,                          markerStyle=20, markerSizes=[0.8,1.0,1.2], ptEdges=ptEdges, markerColor=ROOT.kOrange-7 ),
    MultiGraphMarkerDrawerVsPt( filter = lambda d : d.charge() == 0 and abs(d.pdgId()) == 22, markerStyle=20, markerSizes=[0.8,1.0,1.2], ptEdges=ptEdges, markerColor=ROOT.kSpring+6 ),
    MultiGraphMarkerDrawerVsPt( filter = lambda d : d.charge() == 0 and abs(d.pdgId()) != 22, markerStyle=20, markerSizes=[0.8,1.0,1.2], ptEdges=ptEdges, markerColor=ROOT.kBlue-9 ),
]

chargePdgIdPrinter = lambda p : "charge %+1d pdgId %+4d" % (p.charge(), p.pdgId())

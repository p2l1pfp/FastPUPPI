import ROOT, os
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(True)

class plotTemplate:
    def __init__(self, outdir=None, defaultExts=["png","pdf","eps"]):
        self.outdir = outdir
        if outdir:
            if not os.path.isdir(outdir):
                os.system("mkdir -p "+outdir)
            os.system("cp %s/src/FastPUPPI/NtupleProducer/python/display/index.php %s/" % (os.environ['CMSSW_BASE'], outdir));
        ROOT.gROOT.ProcessLine(".x %s/src/FastPUPPI/NtupleProducer/python/display/tdrstyle.cc" % os.environ['CMSSW_BASE']);
        ROOT.gStyle.SetErrorX(0.5)
        ROOT.gStyle.SetOptStat(0)
        ROOT.gStyle.SetOptFit(0)
        self.canvas = ROOT.TCanvas("c1","c1")
        self._defaultExts = defaultExts 
    def newCanvas(self):
        if self.canvas: del self.canvas
        self.canvas = ROOT.TCanvas("c1","c1")
        return self.canvas
    def decorations(self,energy="14 TeV", pu=200, lumi=False):
        tex = ROOT.TLatex()
        tex.SetTextSize(0.03)
        #tex.SetTextFont(42)
        #tex.DrawLatexNDC(0.17,0.95,"#bf{#scale[1.5]{CMS}} #it{Phase-2 Simulation}")
        tex.DrawLatexNDC(self.canvas.GetLeftMargin(),0.95,"#scale[1.5]{CMS} Phase-2 Simulation")
        tex.SetTextSize(0.035)
        tex.SetTextAlign(31)
        tex.DrawLatexNDC(1.0-self.canvas.GetRightMargin(),0.95,"%s%s, %s PU" % (energy, ", 3000 fb^{-1}" if lumi else "", pu))
    def addSpam(self,x1,y1,text,textSize=0.035,textAlign=21):
        tex = ROOT.TLatex()
        tex.SetTextSize(textSize)
        tex.SetTextFont(42)
        tex.DrawLatexNDC(x1,y1,text)
    def SetLogy(self, logy):
        self.canvas.SetLogy(logy)
    def Print(self,basename, exts=None):
        if exts == None: exts = self._defaultExts;
        errorLevel = ROOT.gErrorIgnoreLevel
        ROOT.gErrorIgnoreLevel = ROOT.kWarning
        if self.outdir: basename = self.outdir + "/" + basename
        for e in exts:
            if (e == "png") and "eps" in exts: continue
            self.canvas.Print(basename+"."+e)
        if "png" in exts and "eps" in exts:
            ret = os.system("LD_LIBRARY_PATH=\"/lib64:$LD_LIBRARY_PATH\" convert -density 150 -quality 100 %s.eps %s.png" % (basename,basename))
            if ret != 0 and not os.path.exists(basename+".png"):
                # fall back to plain root, even if it doesn't work great
                self.canvas.Print(basename+".png")
        ROOT.gErrorIgnoreLevel = errorLevel



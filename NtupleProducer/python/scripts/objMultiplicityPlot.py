import os
from sys import argv
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetErrorX(0.5)
ROOT.gErrorIgnoreLevel = ROOT.kWarning
from math import *
from optparse import OptionParser
from FastPUPPI.NtupleProducer.plotTemplate import plotTemplate

SAMPLE_LABEL = dict(
    TTbar_PU200 = "t#bar{t}",
    VBF_HToInvisible_PU200 = "VBF H #rightarrow inv",
    VBFHToBB_PU200 = "VBF H #rightarrow b#bar{b}",
    TTTT_PU200 = "t#bar{t}t#bar{t}",
    TTTT_PU300 = "t#bar{t}t#bar{t}",
)
SAMPLE_LABEL["SMS_T1tttt_mGluino-2100_mLSP-400_PU200"] = "SUSY tttt"
SAMPLE_LABEL["SMS_T1tttt_mGluino-2100_mLSP-400_PU300"] = "SUSY tttt"
ALL_SUBDETECTORS = ["Barrel","HF","HGCal","HGCalNoTK"]
        
parser = OptionParser("%(prog) infile [ src [ dst ] ]")
parser.add_option("--cl", type=float, dest="cl", default=0.95, help="Compute number to avoid truncations at this CL")
parser.add_option("-p", type="string", dest="particles", action="append", default=[], help="objects to count: Calo, EmCalo, Mu, TK, PF, Puppi [PF,Puppi]Charged, [PF,Puppi]Neutral, [PF,Puppi]ChargedHadron, [PF,Puppi]NeutralHadron, [PF,Puppi]Photon, [PF,Puppi]Electron, [PF,Puppi]Muon; default is all")
parser.add_option("-d", type="string", dest="detectors", action='append', default=[], help="choice of detector: %s, All (sum of everything); default is all" % (", ".join(ALL_SUBDETECTORS)))
parser.add_option("-s", dest="sample", choices=list(SAMPLE_LABEL.keys()), default='TTbar_PU200', help="choice of sample: "+", ".join(list(SAMPLE_LABEL.keys())))
parser.add_option("--CE", dest="hgcalName", default="HGCal", action="store_const", const="CE", help="Label HGCal as CE in the plots")
parser.add_option("--HGCAL", dest="hgcalName", default="HGCal", action="store_const", const="HGCAL", help="Label HGCal as HGCAL in the plots")
parser.add_option("-w", dest="what", choices=["Old","Sec","Reg"], default='Old', help="What to count: Old, Sec (input objs per sector), Ref (new objs per sector)")
parser.add_option("--nmax", dest="nmax", default=None, type=int, help="Fixed range up to nmax")
options, args = parser.parse_args()

# flatten in case commas are used
options.detectors = sum((x.split(",") for x in options.detectors), [])
options.particles = sum((x.split(",") for x in options.particles), [])

odir = args[1] 
plotter = plotTemplate(odir)
plotter.canvas.SetLeftMargin(0.20) # bigger margin 

if options.cl >= 1:
    raise RuntimeError("--cl must take an argument stricly between 0 and 1")

detectors = ALL_SUBDETECTORS+["All"]
particles = [ "Calo", "EmCalo", "Mu", "TK" ]
for Algo in "PF", "Puppi":
    particles.append(Algo)
    for Type in "Charged Neutral ChargedHadron NeutralHadron Photon Electron Muon".split():
        particles.append(Algo+Type)

tfile = ROOT.TFile.Open(args[0])
tree = tfile.Get("ntuple/tree")
if options.what == "Old":
    producer, defprefix = 'l1pfProducer', 'L1'
elif options.what in ("Sec", "Reg"):
    producer, defprefix = 'l1ctLayer1', options.what
    if options.what == "Sec":
        particles = [ p for p in particles if ("PF" not in p and "Puppi" not in p) ]


for detector in detectors:
    if options.detectors and (detector not in options.detectors): continue
    detectorLabel = detector
    detectorFull = producer + detectorLabel
    detectorLabelPlot = "%s" % detectorLabel
    detectorLabelRegionPlot = "%s region" % detectorLabel
    detectorLabelFName = detectorLabel
    if options.hgcalName != "HGCal": 
        detectorLabelPlot = detectorLabelPlot.replace("HGCalNoTK","%s (|#eta|>2.5)" % options.hgcalName)
        detectorLabelPlot = detectorLabelPlot.replace("HGCal", "%s (|#eta|#leq2.5)" % options.hgcalName)
        detectorLabelRegionPlot = detectorLabelRegionPlot.replace("HGCalNoTK region","%s region (|#eta|>2.5)" % options.hgcalName)
        detectorLabelRegionPlot = detectorLabelRegionPlot.replace("HGCal region", "%s region (|#eta|#leq2.5)" % options.hgcalName)
        detectorLabelFName = detectorLabelFName.replace("HGCal",options.hgcalName)
    for particle in particles:
        if options.particles and (particle not in options.particles): continue
        prefix = "" if (options.what == "Reg" and ("PF" in particle or "Puppi" in particle)) else defprefix;
        if options.cl > 0:
            if detectorLabel!='All':
                variable = "%smaxN%s%s"% (detectorFull,prefix,particle)
            else:
                variable = '+'.join("%stotN%s%s" % (producer + subdetectorLabel, prefix, particle) for subdetectorLabel in ALL_SUBDETECTORS)
            n = tree.Draw("min(%s,999)>>htemp(1000,-0.5,999.5)" % (variable), "mc_id == 998", "")
            if not n: continue
            h = ROOT.gROOT.FindObject("htemp")
            acc = 0
            min_obj, min_bin = -1, -1
            for b in range(0,h.GetNbinsX()+2):
                acc += h.GetBinContent(b)
                if acc > n * options.cl: 
                    min_obj = h.GetBinCenter(b)
                    print("%-20s %3.0f" % (particle, min_obj))
                    break
            ## Do some optimized setting of the x range and rebinning
            max_plot_bin = max(3,(min_bin*3/2)); acc = 0;
            for b in range(h.GetNbinsX()+1,0,-1):
                acc += h.GetBinContent(b)
                if acc >= 1e-3*n or (acc > 0 and b <= max(3,max_plot_bin)):
                    max_plot_bin = b 
                    break
            if options.nmax: max_plot_bin = options.nmax
        for x in "max","tot","vec": # IMPORTANT: max must be the first
            if detectorLabel=='All':
                if x != 'tot': continue
                variable = '+'.join("%stotN%s%s" % (producer + subdetectorLabel, prefix, particle) for subdetectorLabel in ALL_SUBDETECTORS)
            else:
                variable = "%s%sN%s%s"% (detectorFull,x,prefix,particle)
            print("plotting %s" %(variable))
            if x == "max": # we already have the histogram
                if ROOT.gROOT.FindObject("htemp2"): ROOT.gROOT.FindObject("htemp2").Delete()
                h_rebin = ROOT.TH1F("htemp2","htemp2", max_plot_bin+1, -0.5, max_plot_bin+0.5)
                # copy bins
                for b in range(1,max_plot_bin):
                    h_rebin.SetBinContent(b, h.GetBinContent(b))
                # fix overflow
                h_rebin.SetBinContent(max_plot_bin+1, n - sum(h.GetBinContent(b) for b in range(1,max_plot_bin))) 
                h = h_rebin
            else:
                if x == "vec":
                    n = tree.Draw("min(%s,%d)>>htemp(%d,-0.5,%g)" % (variable, max_plot_bin, max_plot_bin, max_plot_bin+0.5), "mc_id == 998", "")
                else:
                    n = tree.Draw(variable, "mc_id == 998", "")
                if not n: continue
                h = ROOT.gROOT.FindObject("htemp")
            if x=='tot':
                h.GetYaxis().SetTitle("Events")
                h.GetXaxis().SetTitle("Total %s in %s"%(particle, detectorLabelPlot))
            elif x=='vec':
                h.GetYaxis().SetTitle('Regions')
                h.GetXaxis().SetTitle("%s in %s"%(particle, detectorLabelRegionPlot))
            else:
                h.GetYaxis().SetTitle("Events")
                h.GetXaxis().SetTitle("Max %s per %s"%(particle, detectorLabelRegionPlot))
            h.GetYaxis().SetTitleOffset(1.5)
            h.Draw()
            plotter.SetLogy(x == 'vec')
            if ( (x=='max' and detectorLabel!='All') or (x=='tot' and detectorLabel=='All') ) and options.cl>0:
                h.SetMaximum(1.4*h.GetMaximum())
                plotter.addSpam(0.25, 0.87, SAMPLE_LABEL[options.sample], textSize=0.045)
                plotter.addSpam(0.25, 0.81, "Min. objects (%.0f%% no trunc.): %i"%(options.cl*100., min_obj), textSize=0.045)
            plotter.decorations(pu=(300 if "PU300" in options.sample else 200))
            plotter.Print(particle+"_"+detectorLabelFName+"_"+x, exts=["png","pdf","eps","root"])

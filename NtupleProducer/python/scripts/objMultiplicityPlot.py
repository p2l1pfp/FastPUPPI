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

parser = OptionParser("%(prog) infile [ src [ dst ] ]")
parser.add_option("--cl", type=float, dest="cl", default=0.95, help="Compute number to avoid truncations at this CL")
parser.add_option("-p", type="string", dest="particles", action="append", default=[], help="objects to count: Calo, EmCalo, Mu, TK, PF, Puppi [PF,Puppi]Charged, [PF,Puppi]Neutral, [PF,Puppi]ChargedHadron, [PF,Puppi]NeutralHadron, [PF,Puppi]Photon, [PF,Puppi]Electron, [PF,Puppi]Muon; default is all")
parser.add_option("-d", type="string", dest="detectors", action='append', default=[], help="choice of detector: Barrel, HGCal, HGCalNoTK, HF, All; default is all")
parser.add_option("-s", dest="sample", choices =['TTbar_PU200','VBF_HToInvisible_PU200','VBFHToBB_PU200'], default = 'TTbar_PU200', help="choice of sample: TTbar_PU200, VBF_HToInvisible_PU200, and VBFHToBB_PU200")
options, args = parser.parse_args()

odir = args[1] 
plotter = plotTemplate(odir)
plotter.canvas.SetLeftMargin(0.20) # bigger margin 

if options.cl >= 1:
    raise RuntimeError("--cl must take an argument stricly between 0 and 1")

detectors = ["Barrel","HF","HGCal","HGCalNoTK","All"]
particles = [ "Calo", "EmCalo", "Mu", "TK" ]
for Algo in "PF", "Puppi":
    particles.append(Algo)
    for Type in "Charged Neutral ChargedHadron NeutralHadron Photon Electron Muon".split():
        particles.append(Algo+Type)

tfile = ROOT.TFile.Open(args[0])
tree = tfile.Get("ntuple/tree")
for detector in detectors:
    if options.detectors and (detector not in options.detectors): continue
    detectorLabel = detector
    detectorFull = 'l1pfProducer' + detectorLabel
    for particle in particles:
        if options.particles and (particle not in options.particles): continue
        if options.cl > 0:
            if detectorLabel!='All':
                variable = "%smaxNL1%s"% (detectorFull,particle)
            else:
                variable = ''
                for subdetectorLabel  in ["Barrel","HF","HGCal","HGCalNoTK"]:
                    subdetectorFull = 'l1pfProducer' + subdetectorLabel
                    variable += "%stotNL1%s+"% (subdetectorFull,particle)
                variable= variable[:-1]
            n = tree.Draw("min(%s,999)>>htemp(1000,-0.5,999.5)" % (variable), "mc_id == 998", "")
            if not n: continue
            h = ROOT.gROOT.FindObject("htemp")
            acc = 0
            min_obj = -1
            for b in xrange(0,h.GetNbinsX()+2):
                acc += h.GetBinContent(b)
                if acc > n * options.cl: 
                    min_obj = h.GetBinCenter(b)
                    print "%-20s %3.0f" % (particle, min_obj)
                    break
        for x in "tot","max","vec":
            if detectorLabel!='All':
                variable = "%s%sNL1%s"% (detectorFull,x,particle)
            elif detectorLabel=='All' and x=='tot':
                variable = ''
                for subdetectorLabel  in ["Barrel","HF","HGCal","HGCalNoTK"]:
                    subdetectorFull = 'l1pfProducer' + subdetectorLabel
                    variable += "%stotNL1%s+"% (subdetectorFull,particle)
                variable= variable[:-1]
            print "plotting %s" %(variable)
            n = tree.Draw(variable, "mc_id == 998", "")
            if not n: continue
            h = ROOT.gROOT.FindObject("htemp")
            if x=='tot':
                h.GetYaxis().SetTitle("Events")
                if detectorLabel=='HGCal':
                    xTitle = "Total %s in CE (|#eta|#leq2.5)"%(particle)
                elif detectorLabel=='HGCalNoTK':
                    xTitle = "Total %s in CE (|#eta|>2.5)"%(particle)
                else:
                    xTitle = "Total %s in %s"%(particle, detectorLabel)
                h.GetXaxis().SetTitle(xTitle)
            elif x=='vec':
                h.GetYaxis().SetTitle('Regions')
                if detectorLabel=='HGCal':
                    xTitle = "%s in CE region (|#eta|#leq2.5)"%(particle)
                elif detectorLabel=='HGCalNoTK':
                    xTitle = "%s in CE region (|#eta|>2.5)"%(particle)
                else:
                    xTitle = "%s in %s region"%(particle, detectorLabel)
                h.GetXaxis().SetTitle(xTitle)
            else:
                h.GetYaxis().SetTitle("Events")
                if detectorLabel=='HGCal':
                    xTitle = "Max %s per CE region (|#eta|#leq2.5)"%(particle)
                elif detectorLabel=='HGCalNoTK':
                    xTitle = "Max %s per CE region (|#eta|>2.5)"%(particle)
                else:
                    xTitle = "Max %s per %s region"%(particle, detectorLabel)
                h.GetXaxis().SetTitle(xTitle)
            h.GetYaxis().SetTitleOffset(1.5)
            h.Draw()

            plotter.SetLogy(x == 'vec')
            if ( (x=='max' and detectorLabel!='All') or (x=='tot' and detectorLabel=='All') ) and options.cl>0:
                h.SetMaximum(1.4*h.GetMaximum())
                if options.sample=='TTbar_PU200':
                    spam = "t#bar{t}"
                elif options.sample=='VBF_HToInvisible_PU200':
                    spam  = "VBF H #rightarrow inv"
                else:
                    spam = "VBF H #rightarrow b#bar{b}"
                plotter.addSpam(0.25, 0.87, spam, textSize=0.045)
                plotter.addSpam(0.25, 0.81, "Min. objects (%.0f%% no trunc.): %i"%(options.cl*100., min_obj), textSize=0.045)
            plotter.decorations()
            plotter.Print(particle+"_"+detectorLabel+"_"+x, exts=["png","pdf","eps","root"])

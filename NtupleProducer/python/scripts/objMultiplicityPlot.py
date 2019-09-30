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
import sys
sys.path.insert(0,'./')
from display.examplePlot import *

parser = OptionParser("%(prog) infile [ src [ dst ] ]")
parser.add_option("--cl", type=float, dest="cl", default=0.95, help="Compute number to avoid truncations at this CL")
parser.add_option("-p", type="string", dest="particles", action="append", default=[], help="objects to count: Calo, EmCalo, Mu, TK, PF, Puppi [PF,Puppi]Charged, [PF,Puppi]Neutral, [PF,Puppi]ChargedHadron, [PF,Puppi]NeutralHadron, [PF,Puppi]Photon, [PF,Puppi]Electron, [PF,Puppi]Muon; default is all")
parser.add_option("-d", type="string", dest="detectors", action='append', default=[], help="choice of detector: Barrel, HGCal, HGCalNoTK, HF; default is all")
parser.add_option("-s", dest="sample", choices =['ttbar_200pu','vbf_200pu'], default = 'ttbar_200pu', help="choice of sample: ttbar_200pu and vbf_200pu")
options, args = parser.parse_args()

odir = args[1] 
os.system("mkdir -p "+odir)
#os.system("cp %s/src/FastPUPPI/NtupleProducer/python/display/index.php %s/" % (os.environ['CMSSW_BASE'], odir))
#ROOT.gROOT.ProcessLine(".x %s/src/FastPUPPI/NtupleProducer/python/display/tdrstyle.cc" % os.environ['CMSSW_BASE'])

if options.cl >= 1:
    raise RuntimeError("--cl must take an argument stricly between 0 and 1")

c1 = CreateCanvas("c1",Grid=False)

detectors = ["Barrel","HF","HGCal","HGCalNoTK"]
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
            n = tree.Draw("min(%smaxNL1%s,199)>>htemp(200,-0.5,199.5)" % (detectorFull,particle), "mc_id == 998", "")
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
            print "plotting %s%sNL1%s" % (detectorFull, x, particle)
            n = tree.Draw("%s%sNL1%s" % (detectorFull, x, particle), "mc_id == 998", "")
            if not n: continue
            h = ROOT.gROOT.FindObject("htemp")
            if x=='tot':
                h.GetYaxis().SetTitle("Events")
                h.GetXaxis().SetTitle("Total %s in %s"%(particle, detectorLabel))
            elif x=='vec':
                h.GetYaxis().SetTitle('Regions')
                h.GetXaxis().SetTitle("%s in %s region"%(particle, detectorLabel))
            else:
                h.GetYaxis().SetTitle("Events")
                h.GetXaxis().SetTitle("Max %s per %s region"%(particle, detectorLabel))
            h.Draw()
            if x=='vec': 
                c1.SetLogy(1)
            else:
                c1.SetLogy(0)
            if x=='max' and options.cl>0:
                h.SetMaximum(1.4*h.GetMaximum())
                leg = ROOT.TLegend(0.20,0.76,0.75,0.86)
                leg.SetBorderSize(0)
                leg.SetFillStyle(0)
                leg.SetTextSize(0.03)
                if options.sample=='ttbar_200pu':
                    leg.SetHeader("t#bar{t}, 200 PU")
                else:
                    leg.SetHeader("VBF, 200 PU")
                leg.AddEntry(h, "Min. objects (%.0f%% no trunc.): %i"%(options.cl*100., min_obj))
                leg.Draw("same")
            DrawPrelimLabel(c1)
            DrawLumiLabel(c1, Lumi="")
            out = odir+'/'+particle+"_"+detectorLabel+"_"+x
            SaveCanvas(c1,PlotName=out)
            

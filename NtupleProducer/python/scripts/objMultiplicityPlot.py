import os
from sys import argv
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(False)
ROOT.gStyle.SetErrorX(0.5)
ROOT.gErrorIgnoreLevel = ROOT.kWarning

from math import *

from optparse import OptionParser
parser = OptionParser("%(prog) infile [ src [ dst ] ]")
parser.add_option("--cl", type=float, dest="cl", default=0, help="Compute number to avoid truncations at this CL")
parser.add_option("-p", type="string", dest="particles", action="append", default=[], help="objects to count")
parser.add_option("-d", dest="detector", choices=["Barrel","HF","HGCal","HGCalNoTK"], default="Barrel", help="choice of detector: Barrel, HGCal, HGCalNoTK, HF")
options, args = parser.parse_args()
if options.cl == 0:
    odir = args[1] 
    os.system("mkdir -p "+odir)
    os.system("cp %s/src/FastPUPPI/NtupleProducer/python/display/index.php %s/" % (os.environ['CMSSW_BASE'], odir));
    ROOT.gROOT.ProcessLine(".x %s/src/FastPUPPI/NtupleProducer/python/display/tdrstyle.cc" % os.environ['CMSSW_BASE']);
elif options.cl >= 1:
    raise RuntimeError("--cl must take an argument stricly between 0 and 1")

c1 = ROOT.TCanvas("c1","c1")
particles = [ "Calo", "EmCalo", "Mu", "TK" ]
detectorLabel = options.detector
detectorFull = 'l1pfProducer' + options.detector
for Algo in "PF", "Puppi":
    particles.append(Algo)
    for Type in "Charged Neutral ChargedHadron NeutralHadron Photon Electron Muon".split():
        particles.append(Algo+Type)
tfile = ROOT.TFile.Open(args[0])
tree = tfile.Get("ntuple/tree")
ROOT.gStyle.SetOptStat("omr")
for particle in particles:
    if options.particles and (particle not in options.particles): continue
    if options.cl > 0:
        n = tree.Draw("min(%smaxNL1%s,199)>>htemp(200,-0.5,199.5)" % (detectorFull,particle), "mc_id == 998", "")
        if not n: continue
        h = ROOT.gROOT.FindObject("htemp")
        acc = 0
        for b in xrange(0,h.GetNbinsX()+2):
            acc += h.GetBinContent(b)
            if acc > n * options.cl: 
                print "%-20s %3.0f" % (particle, h.GetBinCenter(b))
                break
    else:
        for x in "tot","max":
            print "plotting %s%sNL1%s" % (detectorFull, x, particle)
            n = tree.Draw("%s%sNL1%s" % (detectorFull, x, particle), "mc_id == 998", "")
            h = ROOT.gROOT.FindObject("htemp")
            if x=='tot':
                h.GetXaxis().SetTitle("Total %s in %s"%(particle, detectorLabel))
            else:
                h.GetXaxis().SetTitle("Max %s per %s region"%(particle, detectorLabel))
            h.GetYaxis().SetTitle("Events")
            h.Draw()
            if not n: continue
            for ext in "pdf","png":
                out = odir+'/'+particle+"_"+detectorLabel+"_"+x+"."+ext
                c1.Print(out)
            tfout = ROOT.TFile.Open(odir+'/'+particle+"_"+detectorLabel+"_"+x+".root", "RECREATE");
            tfout.WriteTObject(h);
            tfout.Close()
        print "plotting %svecNL1%s" % (detectorFull, particle)
        n = tree.Draw("%svecNL1%s" % (detectorFull, particle), "mc_id == 998", "")
        h = ROOT.gROOT.FindObject("htemp")
        h.GetXaxis().SetTitle("%s in %s region"%(particle, detectorLabel))
        h.GetYaxis().SetTitle("Regions")
        h.Draw()
        c1.SetLogy()
        if not n: continue
        for ext in "pdf","png":
            out = odir+'/'+particle+"_"+detectorLabel+"_vec."+ext
            c1.Print(out)
        tfout = ROOT.TFile.Open(odir+'/'+particle+"_"+detectorLabel+"_vec.root", "RECREATE");
        tfout.WriteTObject(h);
        tfout.Close()
            

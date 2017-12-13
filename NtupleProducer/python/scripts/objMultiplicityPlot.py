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
options, args = parser.parse_args()
if options.cl == 0:
    odir = args[1] # "plots/910pre2/test"
    os.system("mkdir -p "+odir)
    os.system("cp %s/src/FastPUPPI/NtupleProducer/python/display/index.php %s/" % (os.environ['CMSSW_BASE'], odir));
    ROOT.gROOT.ProcessLine(".x %s/src/FastPUPPI/NtupleProducer/python/display/tdrstyle.cc" % os.environ['CMSSW_BASE']);
elif options.cl >= 1:
    raise RuntimeError("--cl must take an argument stricly between 0 and 1")
c1 = ROOT.TCanvas("c1","c1")
particles = [ "Calo", "EmCalo", "Mu", "TK" ]
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
        n = tree.Draw("min(maxNL1%s,199)>>htemp(200,-0.5,199.5)" % (particle), "mc_id == 998", "")
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
            print "%sNL1%s" % (x, particle)
            n = tree.Draw("%sNL1%s" % (x, particle), "mc_id == 998", "")
            if not n: continue
            out = odir+'/'+particle+"_"+x+".png"
            c1.Print(out)


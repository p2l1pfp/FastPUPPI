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
options, args = parser.parse_args()

odir = args[1] # "plots/910pre2/test"
os.system("mkdir -p "+odir)
os.system("cp %s/src/FastPUPPI/NtupleProducer/python/display/index.php %s/" % (os.environ['CMSSW_BASE'], odir));
ROOT.gROOT.ProcessLine(".x %s/src/FastPUPPI/NtupleProducer/python/display/tdrstyle.cc" % os.environ['CMSSW_BASE']);
c1 = ROOT.TCanvas("c1","c1")
particles = {
    "Calo":(10000,200),  "EmCalo":(10000,200), "Mu":(100,100), "TK":(1000,1000), "PF":(1000,1000), "PFCharged":(1000,1000), "PFNeutral":(1000,1000), "Puppi":(1000,1000), "PuppiCharged":(1000,1000), "PuppiNeutral":(1000,1000)
}.items()
tfile = ROOT.TFile.Open(args[0])
tree = tfile.Get("ntuple/tree")
ROOT.gStyle.SetOptStat("omr")
for particle, (nmax, nmaxtot) in particles:
    for x,bins in ("tot",nmaxtot),("max",nmax):
        print "%sNL1%s" % (x, particle)
        n = tree.Draw("%sNL1%s" % (x, particle), "mc_id == 998", "")
        if not n: continue
        out = odir+'/'+particle+"_"+x+".png"
        c1.Print(out)


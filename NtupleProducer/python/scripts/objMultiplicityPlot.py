import os
from sys import argv
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(True)
ROOT.gROOT.ProcessLine(".x %s/cpp/tdrstyle.cc" % os.environ['HOME']);
ROOT.gStyle.SetOptStat(False)
ROOT.gStyle.SetErrorX(0.5)
ROOT.gErrorIgnoreLevel = ROOT.kWarning

from math import *

def doRespEta(oname, tree, name, expr, cut):
    return doRespEtaMedian(oname, tree, name, expr, cut)

def quantiles(ys):
    ys.sort()
    ny = len(ys)
    median = ys[ny/2]
    u68    = ys[min(int(ceil(ny*0.84)),ny-1) ]
    l68    = ys[int(floor(ny*0.16))]
    return (median,l68,u68)

def doRespEtaMedian(oname, tree, name, expr, cut, etabins=10, etamax=5.0):
    ys = [[] for ieta in xrange(etabins)]
    npoints = tree.Draw("("+expr+")/mc_pt:abs(mc_eta)", cut, "");
    if npoints <= 0: return None
    graph = ROOT.gROOT.FindObject("Graph");
    xi, yi = graph.GetX(), graph.GetY()
    for i in xrange(graph.GetN()):
        if (xi[i] > etamax): continue
        ieta = int(floor(etabins*xi[i]/etamax))
        ys[ieta].append(yi[i])
    ret = ROOT.TGraphAsymmErrors()
    ret.SetName(name)
    for ieta in xrange(etabins):
        if len(ys[ieta]) == 0: continue
        (median,lo,hi) = quantiles(ys[ieta])
        ipoint = ret.GetN()
        ret.Set(ipoint+1)
        ret.SetPoint(ipoint, (ieta+0.5)*etamax/etabins, median)
        #ret.SetPointError(ieta, 0.5*etamax/etabins, 0.5*etamax/etabins, (median-lo),(hi-median))
        ret.SetPointError(ipoint, 0.5*etamax/etabins, 0.5*etamax/etabins, (median-lo)/sqrt(len(ys[ieta])),(hi-median)/sqrt(len(ys[ieta])))
    return ret
def doRespEtaProf(oname, tree, name, expr, cut):
    tree.Draw("("+expr+")/mc_pt:abs(mc_eta)>>"+name+"(20,0,5.0)", cut, "PROF");
    return ROOT.gROOT.FindObject(name);
def doRespPt(oname, tree, name, expr, cut, resol=False):
    if "jet" in oname: 
        ptbins = [20,25,30,35,40,45,50,55,60,70,80,90,100,120,140,160,200,250]
    else:
        ptbins = [2.5,5,7.5,10,15,20,25,30,35,40,45,50,55,60,70,80,90,100]
    ys = [[] for ipt in ptbins]
    npoints = tree.Draw("("+expr+")/mc_pt:mc_pt", cut, "");
    if npoints <= 0: return None
    graph = ROOT.gROOT.FindObject("Graph");
    xi, yi = graph.GetX(), graph.GetY()
    for i in xrange(graph.GetN()):
        for ipt,ptmax in enumerate(ptbins):
            if xi[i] < ptmax:
                ys[ipt].append(yi[i])
                break
    ret = ROOT.TGraphAsymmErrors()
    ret.SetName(name)
    for ipt,ptmax in enumerate(ptbins):
        if len(ys[ipt]) == 0: continue
        ptmin = ptbins[ipt-1] if ipt else 0
        ptc   = 0.5*(ptmin+ptmax)
        ptd   = 0.5*(ptmax-ptmin)
        (median,lo,hi) = quantiles(ys[ipt])
        ipoint = ret.GetN()
        if not resol:
            ret.Set(ipoint+1)
            ret.SetPoint(ipoint, ptc, median)
            ret.SetPointError(ipoint, ptd, ptd, (median-lo)/sqrt(len(ys[ipt])),(hi-median)/sqrt(len(ys[ipt])))
        else:
            if median <= 0.2: continue
            #print "for %s %s pt %g, median = %g high = %g lo = %g " % (oname, name, ptc, median, hi, lo)
            ret.Set(ipoint+1)
            ret.SetPoint(ipoint, ptc, (hi-lo)/2/median)
            ret.SetPointError(ipoint, ptd, ptd, 0, 0)
    if ret.GetN() <= 3: return None
    print oname, name, ret.GetN()
    if not resol:
        tf1 = ROOT.TF1(name+"_f1","1/x++1", ret.GetX()[0]-ret.GetErrorXlow(0), ret.GetX()[ret.GetN()-1]+ret.GetErrorXhigh(ret.GetN()-1) )
        ret.Fit(tf1, "WQ0C EX0 ROB=0.95")
        ret.fit = tf1
    else:
        tf1 = ROOT.TF1(name+"_f1","1/x++1", ret.GetX()[0]-ret.GetErrorXlow(0), ret.GetX()[ret.GetN()-1]+ret.GetErrorXhigh(ret.GetN()-1) )
        ret.Fit(tf1, "WQ0C EX0 ROB=0.95")
        ret.fit = tf1
    return ret


whats = [
    ('inputs',[
        ("Had",  "Hcal$", ROOT.kAzure+1,  25, 2.0),
        ("Em",   "Ecal$", ROOT.kGreen+2,  21, 1.5),
        ("Calo", "Calo$", ROOT.kViolet+2, 34, 1.5),
        ("Trk",  "TK$",   ROOT.kRed+1, 20, 1.2),
    ]),
    #('inputs-simpleCorr',[
    #    ("Had",    "Hcal$", ROOT.kAzure+1,  25, 2.0),
    #    ("Em-C",   "Ecal$/(1+0.20*(abs(mc_eta)<3)-0.15*(abs(mc_eta)<1.5))", ROOT.kGreen+2,  21, 1.5),
    #    ("Calo-C", "Hcal$ + (abs(mc_eta)<3)*Ecal$/(1+0.20*(abs(mc_eta)<3)-0.15*(abs(mc_eta)<1.5))", ROOT.kViolet+2, 34, 1.5),
    #    ("Trk",    "TK$",   ROOT.kRed+1, 20, 1.2),
    #]),
    ('l1pf',[
        ("Raw Calo",   "L1RawCalo$", ROOT.kViolet-4,  21, 1.7),
        ("Calo",       "L1Calo$",    ROOT.kViolet+2, 34, 1.5),
        ("TK",         "L1TK$",      ROOT.kRed+1, 20, 1.2),
        ("PF",         "L1PF$",      ROOT.kOrange+7, 34, 1.2),
        ("Puppi",      "L1Puppi$",   ROOT.kGray+2, 25, 1.4),
    ]),
    ('il1pf',[
        ("iCalo",       "L1ICalo$",    ROOT.kViolet+2, 34, 1.5),
        ("iTK",         "L1ITK$",      ROOT.kRed+1, 20, 1.2),
        ("iPF",         "L1IPF$",      ROOT.kOrange+7, 34, 1.2),
        ("iPuppi",      "L1IPuppi$",   ROOT.kGray+2, 25, 1.4),
    ]),
]

from optparse import OptionParser
parser = OptionParser("%(prog) infile [ src [ dst ] ]")
options, args = parser.parse_args()

odir = args[1] # "plots/910pre2/test"
os.system("mkdir -p "+odir)
os.system("cp %s/php/index.php %s/" % (os.environ['HOME'], odir));
ROOT.gROOT.ProcessLine(".x %s/cpp/tdrstyle.cc" % os.environ['HOME']);
c1 = ROOT.TCanvas("c1","c1")
particles = {
    "Calo":(10000,200), "Mu":(100,100), "TK":(1000,1000), "PF":(1000,1000), "Puppi":(1000,1000)
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


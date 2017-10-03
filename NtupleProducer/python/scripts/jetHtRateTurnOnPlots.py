
import os, sys
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(True)
ROOT.gROOT.ProcessLine(".x %s/cpp/tdrstyle.cc" % os.environ['HOME']);
ROOT.gStyle.SetOptStat(False)
ROOT.gStyle.SetErrorX(0.5)
ROOT.gErrorIgnoreLevel = ROOT.kWarning

from array import array
from bisect import bisect_left

def makeCumulativeHTEff(name, tree, expr, cut="", norm=40E3):
    tree.Draw("min(%s,1999)>>htemp(2000,0,2000)" % expr, cut);
    htemp = ROOT.gROOT.FindObject("htemp")
    tot, msum = norm/htemp.Integral(), 0
    nbins = htemp.GetNbinsX()
    for ib in xrange(0, nbins-1):
        msum += htemp.GetBinContent(nbins-ib) * tot
        htemp.SetBinContent(nbins-ib, msum)
    ret = htemp.Clone("ceff_"+name)
    ret.SetDirectory(None)
    return ret

def makeEffHist(name, tree, expr, thr, gvar, cut=""):
    tree.Draw("(%s > %s):%s>>htemp(16,0,800,2,-0.5,1.5)" % (expr,thr,gvar), cut);
    htemp = ROOT.gROOT.FindObject("htemp")
    num, den = [ ROOT.TH1F(name+"_"+x,"",16,0,800) for x in ("pass","tot") ]
    for iy,h in (2,num),(1,den):
        for ix in xrange(1,htemp.GetNbinsX()+2):
            h.SetBinContent(ix, htemp.GetBinContent(ix,iy))
    den.Add(num)
    eff = ROOT.TEfficiency(num,den)
    eff.SetStatisticOption(eff.kFCP)
    return eff

whats = [
    ('l1pf',[
        #("Raw Calo",   "RawCalo$", ROOT.kViolet-4,  21, 1.7),
        ("Calo",       "AK4CaloJets$",    ROOT.kViolet+2, 20, 1.5),
        ("TK",         "AK4TKJets$",      ROOT.kRed+0, 21, 1.5),
        ("TK #Deltaz", "AK4TKVJets$",     ROOT.kRed+1, 24, 1.5),
        ("PF",         "AK4PFJets$",      ROOT.kOrange+7, 24, 1.5),
        ("Puppi",      "AK4PuppiJets$",   ROOT.kBlue+1, 21, 1.5),
    ]),
    ('il1pf',[
        ("iCalo",       "L1Calo$",    ROOT.kViolet+2, 20, 1.5),
        ("iTK",         "L1TK$",      ROOT.kRed+1, 21, 1.5),
        ("iPF",         "L1PF$",      ROOT.kOrange+7, 24, 1.5),
        ("iPuppi",      "L1Puppi$",   ROOT.kBlue+1, 25, 1.5),
    ]),
]

from optparse import OptionParser
parser = OptionParser("%(prog) infile [ src [ dst ] ]")
parser.add_option("-w", dest="what", default=None, help="Choose set (il1pf, l1pf, ...)")
parser.add_option("-s", dest="genht",  default=300, type="float", help="Choose gen ht")
parser.add_option("-r", dest="rate",  default=20, type="float", help="Choose rate [kHz]")
parser.add_option("-p", dest="pt",  default=30, type="int", help="Choose pt cut")
parser.add_option("-e", dest="eta",  default="24", help="Choose eta")
options, args = parser.parse_args()

tfiles = [ROOT.TFile.Open(f) for f in args[:2]]
signal, background = [ f.Get("ntuple/tree") for f in tfiles ]
sels = []; fname = args[0] # "respTupleNew_D4T_NoPU.root"

odir = args[2] # "plots/910pre2/test"
os.system("mkdir -p "+odir)
os.system("cp %s/php/index.php %s/" % (os.environ['HOME'], odir));
ROOT.gROOT.ProcessLine(".x %s/cpp/tdrstyle.cc" % os.environ['HOME']);
ROOT.gStyle.SetOptStat(0)
c1 = ROOT.TCanvas("c1","c1")

for kind,things in whats:
    if options.what and (options.what != kind): continue
    plots = []
    htexpr = "E%sPt%d_ht_corr" % (options.eta, options.pt)
    for name,expr,col,msty,msiz in things:
        rexpr = expr.replace("$",htexpr)
        label = name
        if args[3] == "rate":
            plot = makeCumulativeHTEff(name, background, rexpr)
        elif args[3] == "eff":
            plot = makeEffHist(name, signal, rexpr, options.genht, "AK4GenJetsE%sPt%d_ht_raw" % (options.eta, options.pt))
        elif args[3] == "isorate":
            rateplot = makeCumulativeHTEff(name, background, rexpr)
            cut = 9999
            for ix in xrange(1,rateplot.GetNbinsX()+1):
                if rateplot.GetBinContent(ix) <= options.rate:
                    cut = rateplot.GetXaxis().GetBinLowEdge(ix)
                    break
            plot = makeEffHist(name, signal, rexpr, cut, "AK4GenJetsE%sPt%d_ht_raw" % (options.eta, options.pt))
            label = "H_{T}(%s) > %.0f" % (name,cut)
        else: raise RuntimeError
        plot.SetLineWidth(3); plot.SetLineColor(col);  plot.SetMarkerColor(col)
        plot.SetMarkerStyle(msty); plot.SetMarkerSize(msiz)
        plots.append((label,plot))
    if args[3] == "rate":
        c1.SetLogy(True)
        frame = ROOT.TH1D("",";L1 H_{T} cut (p_{T}^{corr} > %.0f, |#eta| < %.1f); Minbias rate @ PU140 [kHz]" % (options.pt, float(options.eta)/10), 100, 0, 1000)
        frame.GetYaxis().SetDecimals(True)
        frame.GetXaxis().SetNdivisions(505)
        frame.GetYaxis().SetRangeUser(0.5, 100e3)
        leg = ROOT.TLegend(0.6,0.99,0.95,0.99-0.06*len(things))
        plotname = 'ht%s-%s_eta%s_pt%d_genHT%.0f' % (args[3], kind, options.eta, options.pt, options.genht)
    elif args[3] == "eff":
        frame = ROOT.TH1D("",";Gen H_{T} (p_{T} > %.0f, |#eta| < %.1f); Eff (L1 H_{T} > %.0f)" % (options.pt, float(options.eta)/10, options.genht), 100, 0, 800)
        leg = ROOT.TLegend(0.6,0.19,0.95,0.19+0.06*len(things))
        plotname = 'ht%s-%s_eta%s_pt%d_genHT%.0f' % (args[3], kind, options.eta, options.pt, options.genht)
    elif args[3] == "isorate":
        frame = ROOT.TH1D("",";Gen H_{T} (p_{T} > %.0f, |#eta| < %.1f); Eff (L1 rate %.0f kHz)" % (options.pt, float(options.eta)/10, options.rate), 100, 0, 800)
        frame.GetYaxis().SetDecimals(True)
        leg = ROOT.TLegend(0.20,0.99,0.55,0.99-0.055*len(things))
        plotname = 'ht%s-%s_eta%s_pt%d_%.0fkHz' % (args[3], kind, options.eta, options.pt, options.rate)
    frame.Draw()
    if args[3] == "rate":
        line = ROOT.TLine()
        line.SetLineStyle(7)
        for y in 40e3, 100, 10:
            line.DrawLine(frame.GetXaxis().GetXmin(),y,frame.GetXaxis().GetXmax(),y)
    for n,p in plots: 
        p.Draw("PCX SAME" if ("TH1" not in p.ClassName()) else "C SAME")
    for n,p in plots: 
        leg.AddEntry(p, n, "LP" if args[3] != "rate" else "L")
    leg.Draw()
    
    c1.Print('%s/%s.png' % (odir, plotname))
    del frame


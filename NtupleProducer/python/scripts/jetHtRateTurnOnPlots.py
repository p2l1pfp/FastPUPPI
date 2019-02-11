
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
from math import pow

def makeCumulativeHTEff(name, tree, expr, xmax, cut="", norm=2760.0*11246/1000):
    ROOT.gROOT.cd()
    htemp = ROOT.gROOT.FindObject("htemp")
    if htemp: htemp.Delete()
    tree.Draw("min(%s,0.9998*%g)>>htemp(2000,0,%g)" % (expr,xmax,xmax), cut);
    htemp = ROOT.gROOT.FindObject("htemp")
    tot, msum = norm/htemp.Integral(), 0
    nbins = htemp.GetNbinsX()
    for ib in xrange(0, nbins):
        msum += htemp.GetBinContent(nbins-ib) * tot
        htemp.SetBinContent(nbins-ib, msum)
    ret = htemp.Clone("ceff_"+name)
    ret.SetDirectory(None)
    return ret

def makeROC(effsig,effbkg):
    graph = ROOT.TGraph(effsig.GetNbinsX())
    for i in xrange(1,graph.GetN()+1):
        graph.SetPoint(i-1, effsig.GetBinContent(i), effbkg.GetBinContent(i))
    return graph 

def makeEffHist(name, tree, expr, thr, gvar, xmax, cut="", logxbins=None):
    ROOT.gROOT.cd()
    htemp = ROOT.gROOT.FindObject("htemp")
    if htemp: htemp.Delete()
    if logxbins:
        nbins, nratio = int(logxbins[0]), float(logxbins[1])
        if nratio == 1:
            htemp = ROOT.TH2D("htemp","htemp",nbins,0,xmax,2,-0.5,1.5)
            num, den = [ ROOT.TH1F(name+"_"+x,"",nbins,0,xmax) for x in ("pass","tot") ]
        else:
            step = pow(nratio, 1.0/nbins)
            base = xmax*(step-1)/(nratio - 1)
            edges = [0]
            for i in xrange(nbins+1):
                edges.append(edges[-1] + base * pow(step, i))
            htemp = ROOT.TH2D("htemp","htemp",nbins,array('f',edges),2,array('f',[-0.5,0.5,1.5]))
            num, den = [ ROOT.TH1F(name+"_"+x,"",nbins,array('f',edges)) for x in ("pass","tot") ]
    else:
        htemp = ROOT.TH2D("htemp","htemp",20,0,xmax,2,-0.5,1.5)
        num, den = [ ROOT.TH1F(name+"_"+x,"",20,0,xmax) for x in ("pass","tot") ]
    tree.Draw("(%s > %s):%s>>htemp" % (expr,thr,gvar), cut);
    htemp = ROOT.gROOT.FindObject("htemp")
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
    ('comp',[
        #("Stage2",     "AK4Stage2Calo$", ROOT.kGreen+2, 21, 1.4),
        ("Calo",       "AK4CaloJets$",    ROOT.kViolet+1, 21, 1.3),
        #("TK",         "AK4TKJets$",      ROOT.kAzure+10, 20, 1.2),
        #("TK (tight)", "AK4TightTKJets$",  ROOT.kAzure+10, 20, 1.2),
        ("TK#Deltaz", "AK4TightTKVJets$",  ROOT.kAzure+2, 20, 1.2),
        #("PF",         "AK4PFJets$",      ROOT.kOrange+7, 20, 1.2),
        ("PF+Puppi",      "AK4PuppiJets$",   ROOT.kRed+1, 20, 1.2),
        ]),
    ('compall',[
        ("Calo",       "AK4CaloJets$",    ROOT.kViolet+1, 21, 1.3),
        ("TK",         "AK4TKJets$",      ROOT.kAzure+10, 20, 1.2),
        ("TK (tight)", "AK4TightTKJets$",  ROOT.kAzure+10, 20, 1.2),
        ("TK #Deltaz", "AK4TightTKVJets$",  ROOT.kAzure+2, 20, 1.2),
        ("PF",         "AK4PFJets$",      ROOT.kOrange+7, 20, 1.2),
        ("PF+Puppi",      "AK4PuppiJets$",   ROOT.kRed+1, 20, 1.2),
        ]),
    ('comptp',[
        #("Stage2",     "AK4Stage2Calo$", ROOT.kGray+1, 21, 1.4),
        ("TrackerVtx",     "TP_TrackerHTVtx",    ROOT.kBlue+0, 21, 1.4),
        #("TkCaloVtx",     "TP_TkCaloHTVtx",    ROOT.kGreen+2, 21, 1.4),
        ("Calo",       "AK4CaloJets$",    ROOT.kViolet+1, 21, 1.3),
        #("TK",         "AK4TKJets$",      ROOT.kAzure+10, 20, 1.2),
        #("TK", "AK4TightTKJets$",  ROOT.kAzure+10, 20, 1.2),
        ("TK#Deltaz", "AK4TightTKVJets$",  ROOT.kAzure+2, 20, 1.2),
        #("PF",         "AK4PFJets$",      ROOT.kOrange+7, 20, 1.2),
        ("PF+Puppi",      "AK4PuppiJets$",   ROOT.kRed+1, 20, 1.2),
        ]),
    ('l1met',[
        ("Calo",       "$Calo",     ROOT.kViolet+2, 20, 1.5),
        ("TK#Deltaz",  "$TightTKV", ROOT.kRed+1, 24, 1.5),
        ("PF",         "$PF",       ROOT.kOrange+7, 24, 1.5),
        ("Puppi",      "$Puppi",    ROOT.kBlue+1, 21, 1.5),
    ]),
    ('l1metc',[
        ("Calo",       "$CaloCentral",  ROOT.kViolet+2, 20, 1.5),
        ("TK#Deltaz",  "$TightTKV",     ROOT.kRed+1, 24, 1.5),
        ("PF",         "$PFCentral",    ROOT.kOrange+7, 24, 1.5),
        ("Puppi",      "$PuppiCentral", ROOT.kBlue+1, 21, 1.5),
    ]),
    ('l1metmhtc',[
        ("Calo",       "min(METCaloCentral,AK4CaloJets$)",  ROOT.kViolet+2, 20, 1.5),
        ("TK#Deltaz",  "min(METTightTKV,AK4TightTKVJets$)",     ROOT.kRed+1, 24, 1.5),
        ("PF",         "min(METPFCentral,AK4PFJets$)",    ROOT.kOrange+7, 24, 1.5),
        ("Puppi",      "min(METPuppiCentral,AK4PuppiJets$)", ROOT.kBlue+1, 21, 1.5),
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
parser.add_option("-x", dest="htvar",  default=("ht","H_{T}"), nargs=2, help="Choose variable")
parser.add_option("--xmax", dest="xmax",  default=1000., type=float, help="Choose variable")
parser.add_option("--logxbins", dest="logxbins",  default=None, nargs=2, type=float, help="--logxbins N X will make N bins, the last being a factor X larger than the first")
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
    if "METMHT" in options.htvar[0]:
        genmet = options.htvar[0].replace("MHT","").replace("MET","METGen")
        htexpr = "E%sPt%d_%s_corr" % (options.eta, options.pt, "mht")
        genexpr = "min(%s, AK4GenJetsE%sPt%d_%s_raw)" % (genmet, options.eta, options.pt, "mht")
        qualif = "p_{T}^{corr} > %.0f, |#eta| < %.1f" % (options.pt, float(options.eta)/10)
    elif "MET" in options.htvar[0]:
        htexpr ="MET" 
        genexpr = options.htvar[0].replace("MET","METGen")
        qualif = "|#eta| < 2.4" if "Central" in  options.htvar[0] else "|#eta| < 4.7" 
    else:
        pre, post = options.htvar[0], ""
        if "[" in options.htvar[0]:
            pre, post = options.htvar[0].split("[")
            post = "["+post
        htexpr = "E%sPt%d_%s_corr%s" % (options.eta, options.pt, pre, post)
        genexpr = "AK4GenJetsE%sPt%d_%s_raw%s" % (options.eta, options.pt, pre, post)
        qualif = "p_{T}^{corr} > %.0f, |#eta| < %.1f" % (options.pt, float(options.eta)/10)
    for name,expr,col,msty,msiz in things:
        if "$" not in expr and options.htvar[0] != "ht": continue
        rexpr = expr.replace("$",htexpr)
        label = name
        if args[3] == "rate":
            plot = makeCumulativeHTEff(name, background, rexpr, options.xmax)
        elif args[3] == "effc":
            cut = "%s > %g" % (genexpr, options.genht) if options.genht else ""
            plot = makeCumulativeHTEff(name, signal, rexpr, options.xmax, norm=1, cut=cut)
        #elif args[3] == "eff":
        #    plot = makeEffHist(name, signal, rexpr, options.genht, genexpr, options.xmax, logxbins=options.logxbins)
        elif args[3] == "isorate":
            rateplot = makeCumulativeHTEff(name, background, rexpr, options.xmax)
            cut = 9999
            for ix in xrange(1,rateplot.GetNbinsX()+1):
                if rateplot.GetBinContent(ix) <= options.rate:
                    cut = rateplot.GetXaxis().GetBinLowEdge(ix)
                    break
            plot = makeEffHist(name, signal, rexpr, cut, genexpr, options.xmax, logxbins=options.logxbins)
            label = "%s(%s) > %.0f" % (options.htvar[1], name,cut)
        elif args[3] == "roc":
            cut = "%s > %g" % (genexpr, options.genht) if options.genht else ""
            effsig  = makeCumulativeHTEff(name+"_s", signal, rexpr, options.xmax, norm=1, cut=cut)
            ratebkg = makeCumulativeHTEff(name+"_b", background, rexpr, options.xmax)
            plot = makeROC(effsig,ratebkg)
            msty = 0
        else: raise RuntimeError
        plot.SetLineWidth(3); plot.SetLineColor(col);  plot.SetMarkerColor(col)
        plot.SetMarkerStyle(msty); plot.SetMarkerSize(msiz)
        plots.append((label,plot))
    if args[3] == "rate":
        c1.SetLogy(True)
        frame = ROOT.TH1D("",";L1 %s cut (%s); Minbias rate @ PU200 [kHz]" % (options.htvar[1], qualif), 100, 0, options.xmax)
        frame.GetYaxis().SetDecimals(True)
        frame.GetXaxis().SetNdivisions(505)
        frame.GetYaxis().SetRangeUser(0.5, 100e3)
        leg = ROOT.TLegend(0.6,0.99,0.95,0.99-0.06*len(things))
        plotname = '%s%s-%s_eta%s_pt%d' % (options.htvar[0], args[3], kind, options.eta, options.pt)
    elif args[3] == "effc":
        gentext, genpost = "iciency", "" 
        if options.genht > 0:
            gentext = " (Gen %s > %s)" % (options.htvar[1], options.genht)
            gentpost = "_gen%.0f" % (options.genht)
        frame = ROOT.TH1D("",";L1 %s thresh (%s); Eff%s" % (options.htvar[1], qualif, gentext), 100, 0, options.xmax)
        leg = ROOT.TLegend(0.6,0.19,0.95,0.19+0.06*len(things))
        plotname = '%s%s-%s_eta%s_pt%d%s' % (options.htvar[0], args[3], kind, options.eta, options.pt, genpost)
    #elif args[3] == "eff":
    #   gentext, genpost = "", "" 
    #   if options.genht > 0:
    #       gentext = "(Gen %s > %s)" % (options.htvar[1], options.genht)
    #       gentpost = "_gen%.0f" % (options.genht)
    #   frame = ROOT.TH1D("",";Gen %s (%s); Eff %s" % (options.htvar[1], qualif, gentext), 100, 0, options.xmax)
    #   leg = ROOT.TLegend(0.6,0.19,0.95,0.19+0.06*len(things))
    #   plotname = '%s%s-%s_eta%s_pt%d%s' % (options.htvar[0], args[3], kind, options.eta, options.pt, genpost)
    elif args[3] == "isorate":
        frame = ROOT.TH1D("",";Gen %s (%s); Eff (L1 rate %.0f kHz)" % (options.htvar[1], qualif, options.rate), 100, 0, options.xmax)
        frame.GetYaxis().SetDecimals(True)
        #leg = ROOT.TLegend(0.20,0.99,0.55,0.99-0.055*len(things))
        leg = ROOT.TLegend(0.65,0.19,0.99,0.19+0.065*len(things))
        plotname = '%s%s-%s_eta%s_pt%d_%.0fkHz' % (options.htvar[0], args[3], kind, options.eta, options.pt, options.rate)
    elif args[3] == "roc":
        c1.SetLogy(True)
        gentext, genpost = "iciency", "" 
        if options.genht > 0:
            gentext = " (Gen %s > %s)" % (options.htvar[1], options.genht)
            gentpost = "_gen%.0f" % (options.genht)
        frame = ROOT.TH1D("",";Eff%s; Minbias rate @ PU200 [kHz]" % (gentext), 100, 0, 1)
        frame.GetYaxis().SetDecimals(True)
        frame.GetXaxis().SetNdivisions(505)
        frame.GetYaxis().SetRangeUser(0.5, 100e3)
        leg = ROOT.TLegend(0.2,0.98,0.55,0.98-0.06*len(things))
        plotname = '%s%s-%s_eta%s_pt%d%s' % (options.htvar[0], args[3], kind, options.eta, options.pt, genpost)
    frame.Draw()
    if args[3] in ("rate","roc"):
        line = ROOT.TLine()
        line.SetLineStyle(7)
        for y in 40e3, 100, 10:
            line.DrawLine(frame.GetXaxis().GetXmin(),y,frame.GetXaxis().GetXmax(),y)
    for n,p in plots: 
        p.Draw("PCX SAME" if ("TH1" not in p.ClassName()) else "C SAME")
    for n,p in plots: 
        leg.AddEntry(p, n, "L" if args[3] in ("rate","roc") else "LP")
    leg.Draw()
    
    c1.Print('%s/%s.png' % (odir, plotname))
    del frame


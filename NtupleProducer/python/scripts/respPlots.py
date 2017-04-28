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
parser.add_option("-w", dest="what",     default=None, help="Choose set (inputs, l1pf, ...)")
parser.add_option("-p", dest="particle", default=None, help="Choose particle (electron, ...)")
options, args = parser.parse_args()

sels = []; fname = args[0] # "respTupleNew_D4T_NoPU.root"
for (particle, pdgIdCut, minPt, maxEta) in [ 
        ("pion", "abs(mc_id) == 211", 10, 5),
        ("photon", "abs(mc_id) == 22", 10, 5),
        ("electron", "abs(mc_id) == 11", 10, 5),
        ("muon", "abs(mc_id) == 13", 10, 5),
        ("tau", "(abs(mc_id) == 15 || abs(mc_id) == 211)", 20, 5),
        ("jet", "abs(mc_id) == 0", 30, 5),
        ("null", "abs(mc_id) == 999", 0, 5)
    ]:
    if options.particle and (options.particle not in particle): 
        continue
    sels.append(("%s_pt_%2d_inf" % (particle, minPt), fname, "mc_pt > %g && %s" % (minPt, pdgIdCut)))
    if "null" in particle: continue; # not point in profiling random cones vs pt
    sels.append(("%s_eta_00_15"     % (particle),        fname, "abs(mc_eta) < 1.5                      && %s" % (pdgIdCut)))
    sels.append(("%s_eta_00_15_res" % (particle),        fname, "abs(mc_eta) < 1.5                      && %s" % (pdgIdCut)))
    sels.append(("%s_eta_15_25"  % (particle),        fname, "abs(mc_eta) > 1.5 && abs(mc_eta) < 2.5 && %s" % (pdgIdCut)))
    sels.append(("%s_eta_15_25_res"  % (particle),        fname, "abs(mc_eta) > 1.5 && abs(mc_eta) < 2.5 && %s" % (pdgIdCut)))
    sels.append(("%s_eta_25_30"  % (particle),        fname, "abs(mc_eta) > 2.5 && abs(mc_eta) < 3.0 && %s" % (pdgIdCut)))
    sels.append(("%s_eta_25_30_res"  % (particle),        fname, "abs(mc_eta) > 2.5 && abs(mc_eta) < 3.0 && %s" % (pdgIdCut)))
    sels.append(("%s_eta_30_50"  % (particle),        fname, "abs(mc_eta) > 3.0 && abs(mc_eta) < 5.0 && %s" % (pdgIdCut)))
    sels.append(("%s_eta_30_50_res"  % (particle),        fname, "abs(mc_eta) > 3.0 && abs(mc_eta) < 5.0 && %s" % (pdgIdCut)))

odir = args[1] # "plots/910pre2/test"
os.system("mkdir -p "+odir)
os.system("cp %s/php/index.php %s/" % (os.environ['HOME'], odir));
ROOT.gROOT.ProcessLine(".x %s/cpp/tdrstyle.cc" % os.environ['HOME']);
c1 = ROOT.TCanvas("c1","c1")
for oname,fname,cut in sels:
    tfile = ROOT.TFile.Open(fname)
    tree = tfile.Get("ntuple/tree")
    if "electron" in oname or "muon" in oname or "pion" in oname:
        cut += " && abs(mc_iso04) < 0.05" # isolated
    for kind,things in whats:
        if options.what and (options.what != kind): 
            continue
        if "tau" in oname:
            ptdefs = [ "pt", "pt02" ]
        elif "jet" in oname:
            ptdefs = [ "pt" ]
        elif "muon" in oname:
            ptdefs = [ "pthighest" ]
        else:
            ptdefs = [ "pt02", "pthighest" ]   
        for ptdef in ptdefs:
            resps = [] 
            for name,expr,col,msty,msiz in things:
                exprptdef = expr.replace("$","_"+ptdef)
                if "eta_25" in oname or "eta_30" in oname:
                    if "TK" in expr: continue
                if "Puppi" in name and "PU0" in odir: continue
                if "pt" in oname:
                    prof = doRespEta(oname,tree,name,exprptdef,cut)
                    if not prof: continue
                else:
                    prof = doRespPt(oname,tree,name,exprptdef,cut,resol=("res" in oname))
                    if not prof: continue
                    if getattr(prof,'fit',None):
                        prof.fit.SetLineWidth(2); prof.fit.SetLineColor(col)
                prof.SetLineWidth(3); prof.SetLineColor(col);  prof.SetMarkerColor(col)
                prof.SetMarkerStyle(msty); prof.SetMarkerSize(msiz)
                resps.append((name,prof))
            if not resps: 
                print oname, ptdef 
                continue
            c1.SetLogy(False)
            if "pt" in oname: 
                frame = ROOT.TH1F("stk","stk",100,0.0,5.0)
                frame.GetYaxis().SetRangeUser(0,2.2)
                frame.GetXaxis().SetTitle("|#eta|")
                leg = ROOT.TLegend(0.6,0.99,0.95,0.99-0.05*len(things))
                frame.GetYaxis().SetTitle("median p_{T}^{rec}/p_{T}^{gen}")
            else:
                frame = ROOT.TH1F("stk","stk",100,0.0,250.0 if "jet" in oname else 100.0)
                if "res" in oname:
                    frame.GetYaxis().SetTitle("#sigma_{eff}(p_{T}^{rec})/p_{T}^{rec}")
                    frame.GetYaxis().SetRangeUser(0.008,10.0)
                    c1.SetLogy(True)
                else:
                    frame.GetYaxis().SetTitle("median p_{T}^{rec}/p_{T}^{gen}")
                    frame.GetYaxis().SetRangeUser(0,2.2)
                frame.GetXaxis().SetTitle("p_{T} (GeV)")
                leg = ROOT.TLegend(0.2,0.99,0.95,0.99-0.05*len(things))
                leg.SetNColumns(2)
            frame.GetYaxis().SetDecimals(True)
            frame.Draw()
            line = ROOT.TLine()
            line.SetLineStyle(7)
            if "res" not in oname:
                line.DrawLine(0.0,1,frame.GetXaxis().GetXmax(),1)
            if "pt" in oname: 
                line.DrawLine(1.5,0,1.5,2.2)
                line.DrawLine(2.5,0,2.5,2.2)
                line.DrawLine(3.0,0,3.0,2.2)
            for n,p in resps: 
                p.Draw("P SAME" if "TGraph" in p.ClassName() else "SAME")
                if hasattr(p,'fit'): p.fit.Draw("SAME")
            for n,p in resps: 
                leg.AddEntry(p, n, "LP")
                if hasattr(p,'fit'): 
                    if "res" not in oname:
                        eq = ("%.2f #timesp_{T} %+ 6.1f" % (p.fit.GetParameter(1), p.fit.GetParameter(0))).replace("-","#minus")
                    else:
                        eq = ("%.2f #timesp_{T} %+ 6.1f" % (p.fit.GetParameter(1), p.fit.GetParameter(0))).replace("-","#minus")
                    leg.AddEntry(p.fit, eq, "L")
            leg.Draw()
            out = odir+'/'+oname+"-"+kind+"_"+ptdef+".png"
            c1.Print(out)
            del frame


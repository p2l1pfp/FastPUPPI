import os
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kWarning

from math import ceil, floor, sqrt

def doRespEta(oname, tree, name, expr, cut):
    return doRespEtaMedian(oname, tree, name, expr, cut)

def quantiles(ys):
    ys.sort()
    ny = len(ys)
    median = ys[ny/2]
   #if ny > 400e9:
   #    u95 = ys[min(int(ceil(ny*0.975)),ny-1) ]
   #    l95 = ys[int(floor(ny*0.025))]
   #    u68 = 0.5*(median+u95)
   #    l68 = 0.5*(median+l95)
    if ny > 20:
        u68 = ys[min(int(ceil(ny*0.84)),ny-1) ]
        l68 = ys[int(floor(ny*0.16))]
    else:
        rms = sqrt(sum((y-median)**2 for y in ys)/ny)
        u68 = median + rms
        l68 = median - rms
    return (median,l68,u68)

def doRespEtaMedian(oname, tree, name, expr, cut, etabins=10, etamax=5.0):
    ys = [[] for ieta in xrange(etabins)]
    npoints = tree.Draw("("+expr+")/mc_pt:abs(mc_eta)", cut, "");
    if npoints <= 0: return None
    graph = ROOT.gROOT.FindObject("Graph");
    xi, yi = graph.GetX(), graph.GetY()
    for i in xrange(graph.GetN()):
        if yi[i] == 0: continue
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
def doRespPt(oname, tree, name, expr, cut, mcpt="mc_pt", fitopt="WQ0C EX0 ROB=0.95"):
    if "jet" in oname: 
        ptbins = [20,25,30,35,40,45,50,55,60,70,80,90,100,120,140,160,200,250]
    else:
        ptbins = [2.5,5,7.5,10,15,20,25,30,35,40,45,50,55,60,70,80,90,100]
    ys = [[] for ipt in ptbins]
    npoints = tree.Draw("("+expr+")/("+mcpt+"):"+mcpt, cut, "");
    if npoints <= 0: return (None,None)
    graph = ROOT.gROOT.FindObject("Graph");
    xi, yi = graph.GetX(), graph.GetY()
    for i in xrange(graph.GetN()):
        if yi[i] == 0: continue
        for ipt,ptmax in enumerate(ptbins):
            if xi[i] < ptmax:
                ys[ipt].append(yi[i])
                break
    ret = ROOT.TGraphAsymmErrors()
    ret.SetName(name)
    resols = []
    for ipt,ptmax in enumerate(ptbins):
        if len(ys[ipt]) == 0: continue
        ptmin = ptbins[ipt-1] if ipt else 0
        ptc   = 0.5*(ptmin+ptmax)
        ptd   = 0.5*(ptmax-ptmin)
        (median,lo,hi) = quantiles(ys[ipt])
        ipoint = ret.GetN()
        ret.Set(ipoint+1)
        ret.SetPoint(ipoint, ptc, median)
        ret.SetPointError(ipoint, ptd, ptd, (median-lo)/sqrt(len(ys[ipt])),(hi-median)/sqrt(len(ys[ipt])))
        if median > 0.2:
            resols.append((ptc,ptd,(hi-lo)/2,0))
            #print "for %s %s at pt ~ %6.1f  reco pt ~ %6.1f  [ %6.1f , %6.1f ]" % (oname,name, ptc, median, lo, hi)
    if ret.GetN() <= 3: return (None,None)
    tf1 = ROOT.TF1(name+"_f1","1/x++1", ret.GetX()[0]-ret.GetErrorXlow(0), ret.GetX()[ret.GetN()-1]+ret.GetErrorXhigh(ret.GetN()-1) )
    ret.Fit(tf1, fitopt)
    ret.fit = tf1
    ## Now we do the approximate response corrected resolution
    scale = ret.fit.GetParameter(1)
    if scale < 0.1 or len(resols) < 3: return (ret,None)
    resol = ROOT.TGraphAsymmErrors(len(resols))
    resol.SetName(name+"_res")
    for ipoint,(ptc,ptd,res,err) in enumerate(resols):
        resol.SetPoint(ipoint, ptc, res/scale)
        resol.SetPointError(ipoint, ptd, ptd, 0, 0)
    if ("Trk" in name or "TK" in name) and "jet" not in oname and "ele" not in oname:
        rtf1 = ROOT.TF1(name+"_rf1","sqrt([0]*[0]+[1]*[1]*(x/1000)*(x/1000))", ret.GetX()[0]-ret.GetErrorXlow(0), ret.GetX()[ret.GetN()-1]+ret.GetErrorXhigh(ret.GetN()-1) )
        rtf1.SetParameters(0.01, 1)
    else:
        rtf1 = ROOT.TF1(name+"_rf1","1/x++1", ret.GetX()[0]-ret.GetErrorXlow(0), ret.GetX()[ret.GetN()-1]+ret.GetErrorXhigh(ret.GetN()-1) )
    resol.Fit(rtf1, fitopt+"W")
    resol.fit = rtf1
    return (ret,resol)

whats = [
    ('inputs',[ # On RelVals
        ("Had",  "Hcal$", ROOT.kAzure+1,  25, 2.0),
        ("Em",   "Ecal$", ROOT.kGreen+2,  21, 1.5),
        ("Calo", "Calo$", ROOT.kViolet+2, 34, 1.5),
        ("Trk",  "TK$",   ROOT.kRed+1, 20, 1.2),
    ]),
    ('debug-1',[
        ("EcalT",  "EcalT$", ROOT.kAzure+1,  25, 2.0),
        ("EcalC",  "EcalC$", ROOT.kAzure+2,  21, 1.6),
        ("HGC/E",  "HGCalE$", ROOT.kAzure+3,  21, 1.5),
        ("Hcal",  "Hcal$", ROOT.kGreen+1,  25, 2.0),
        ("HGC/H",  "HGCalH$", ROOT.kGreen+2,  25, 2.0),
        ("HF", "HF$", ROOT.kViolet+2, 34, 1.5),
    ]),
    ('debug-2',[
        ("EcalT",  "TPEcalT$", ROOT.kAzure+1,  25, 2.0),
        ("EcalC",  "TPEcalC$", ROOT.kAzure+2,  21, 1.6),
        ("Hcal",  "TPHcal$", ROOT.kGreen+3,  25, 2.0),
        ("HGC/TC",  "TPHGCalTC$", ROOT.kGreen+2,  25, 2.0),
        ("HGC/3D",  "TPHGCal3D$", ROOT.kGreen+1,  21, 1.5),
        ("TK", "TPL1Tk$", ROOT.kRed+2, 34, 1.5),
    ]),
    ('OfflineInputs',[
        ("Had",  "Hcal$+HGCalH$+HF$", ROOT.kAzure+1,  25, 2.0),
        ("Em",   "Em$", ROOT.kGreen+2,  21, 1.5),
        ("Calo", "Calo$", ROOT.kViolet+2, 34, 1.5),
    ]),
    ('TPs',[
        ("Had",  "TPHcal$", ROOT.kAzure+1,  25, 2.0),
        ("Em",   "TPEcal$", ROOT.kGreen+2, 21, 1.5),
        ("Calo", "TPCalo$", ROOT.kViolet+2, 34, 1.5),
        ("Trk",  "TPTK$",   ROOT.kRed+1, 20, 1.2),
    ]),
    ('l1pf',[
        ("Gen #times Acc",        "GenAcc$",    ROOT.kAzure+1,  20, 1.2),
        ("Raw Calo",   "L1RawCalo$", ROOT.kViolet-4,  21, 1.7),
        ("Calo",       "L1Calo$",    ROOT.kViolet+2, 34, 1.5),
        ("TK",         "L1TK$",      ROOT.kRed+1, 20, 1.2),
        ("PF",         "L1PF$",      ROOT.kOrange+7, 34, 1.2),
        ("Puppi",      "L1Puppi$",   ROOT.kGray+2, 25, 1.4),
    ]),
    ('il1pf',[
        ("Gen #times Acc",         "GenAcc$",     ROOT.kAzure+1,  20, 1.2),
        ("iCalo",       "L1ICalo$",    ROOT.kViolet+2, 34, 1.5),
        ("iTK",         "L1ITK$",      ROOT.kRed+1, 20, 1.2),
        ("iPF",         "L1IPF$",      ROOT.kOrange+7, 34, 1.2),
        ("iPuppi",      "L1IPuppi$",   ROOT.kGray+2, 25, 1.4),
    ]),
]

if __name__ == "__main__":
    from optparse import OptionParser
    parser = OptionParser("%(prog) infile [ src [ dst ] ]")
    parser.add_option("-w", dest="what",     default=None, help="Choose set (inputs, l1pf, ...)")
    parser.add_option("-p", dest="particle", default=None, help="Choose particle (electron, ...)")
    options, args = parser.parse_args()
    selparticles = options.particle.split(",") if options.particle else []

    sels = []; 
    for (particle, pdgIdCut, minPt, maxEta) in [ 
            ("pion", "abs(mc_id) == 211", 10, 5),
            ("pizero", "abs(mc_id) == 111", 10, 5),
            ("pimix", "(abs(mc_id) == 211 || (abs(mc_id) == 111 && (event % 2) == 1))", 10, 5),
            ("photon", "abs(mc_id) == 22", 10, 5),
            ("electron", "abs(mc_id) == 11", 10, 5),
            ("muon", "abs(mc_id) == 13", 10, 5),
            ("tau", "(abs(mc_id) == 15 || abs(mc_id) == 211)", 20, 5),
            ("jet", "abs(mc_id) == 0", 30, 5),
            ("null", "abs(mc_id) == 999", 0, 5)
        ]:
        if options.particle and (particle not in selparticles): 
            continue
        sels.append(("%s_pt_%2d_inf" % (particle, minPt), "mc_pt > %g && %s" % (minPt, pdgIdCut)))
        if "null" in particle: continue; # not point in profiling random cones vs pt
        sels.append(("%s_eta_00_15"     % (particle),     "abs(mc_eta) < 1.5                      && %s" % (pdgIdCut)))
        sels.append(("%s_eta_15_25"  % (particle),        "abs(mc_eta) > 1.5 && abs(mc_eta) < 2.5 && %s" % (pdgIdCut)))
        sels.append(("%s_eta_25_30"  % (particle),        "abs(mc_eta) > 2.5 && abs(mc_eta) < 3.0 && %s" % (pdgIdCut)))
        sels.append(("%s_eta_30_50"  % (particle),        "abs(mc_eta) > 3.0 && abs(mc_eta) < 5.0 && %s" % (pdgIdCut)))
        sels.append(("%s_eta_30_45"  % (particle),        "abs(mc_eta) > 3.0 && abs(mc_eta) < 4.5 && %s" % (pdgIdCut)))

    tree = ROOT.TChain("ntuple/tree")
    print args
    for a in args[:]:
        if a.endswith(".root"): 
            tree.Add(a)
            args.remove(a)
    print args
    odir = args[0] # "plots/910pre2/test"
    os.system("mkdir -p "+odir)
    os.system("cp %s/src/FastPUPPI/NtupleProducer/python/display/index.php %s/" % (os.environ['CMSSW_BASE'], odir));
    ROOT.gROOT.ProcessLine(".x %s/src/FastPUPPI/NtupleProducer/python/display/tdrstyle.cc" % os.environ['CMSSW_BASE']);
    ROOT.gStyle.SetOptStat(False)
    ROOT.gStyle.SetOptFit(False)
    ROOT.gStyle.SetErrorX(0.5)
    c1 = ROOT.TCanvas("c1","c1")
    for oname,cut in sels:
        print "Plotting ",oname
        if "electron" in oname or "muon" in oname or "pi" in oname:
            cut += " && abs(mc_iso04) < 0.05" # isolated
        for kind,things in whats:
            if options.what and (kind not in options.what.split(",")): 
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
                resps = []; resols = [] 
                for name,expr,col,msty,msiz in things:
                    exprptdef = expr.replace("$","_"+ptdef)
                    if "eta_25" in oname or "eta_30" in oname:
                        if "TK" in expr: continue
                    if "Puppi" in name and "PU0" in odir: continue
                    if name == "Gen" and  "jet" not in oname: continue
                    if "pt" in oname:
                        prof, pres = doRespEta(oname,tree,name,exprptdef,cut), None
                    else:
                        prof, pres = doRespPt(oname,tree,name,exprptdef,cut)
                        #print oname, kind, ptdef, name, prof, pres
                    for (p,ps) in (prof,resps), (pres,resols):
                        if not p: continue
                        if getattr(p,'fit',None):
                            p.fit.SetLineWidth(2); p.fit.SetLineColor(col)
                        p.SetLineWidth(3); p.SetLineColor(col);  p.SetMarkerColor(col)
                        p.SetMarkerStyle(msty); p.SetMarkerSize(msiz)
                        ps.append((name,p))
                for plots,ptype,pfix in (resps,"response",""),(resols,"resolution","_res"):
                    if not plots: 
                        if "_pt" in oname and ptype == "resolution": continue # not implemented
                        print "No ",ptype," plot for ", oname, ptdef 
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
                        if "resolution" in ptype:
                            frame.GetYaxis().SetTitle("#sigma_{eff}(p_{T}^{corr})/p_{T}^{corr}")
                            frame.GetYaxis().SetRangeUser(0.008,10.0)
                            c1.SetLogy(True)
                        else:
                            frame.GetYaxis().SetTitle("median p_{T}^{rec}/p_{T}^{gen}")
                            frame.GetYaxis().SetRangeUser(0,2.2)
                        frame.GetXaxis().SetTitle("p_{T} (GeV)")
                        leg = ROOT.TLegend(0.2,0.99,0.95,0.99-0.05*(len(plots)+0.5))
                        leg.SetTextSize(0.04);
                        leg.SetNColumns(2)
                    frame.GetYaxis().SetDecimals(True)
                    frame.Draw()
                    line = ROOT.TLine()
                    line.SetLineStyle(7)
                    if "resolution" not in ptype:
                        line.DrawLine(0.0,1,frame.GetXaxis().GetXmax(),1)
                    if "pt" in oname: 
                        line.DrawLine(1.5,0,1.5,2.2)
                        line.DrawLine(2.5,0,2.5,2.2)
                        line.DrawLine(3.0,0,3.0,2.2)
                    for n,p in plots: 
                        p.Draw("P SAME" if "TGraph" in p.ClassName() else "SAME")
                        if hasattr(p,'fit'): p.fit.Draw("SAME")
                    for n,p in plots: 
                        leg.AddEntry(p, n, "LP")
                        if hasattr(p,'fit'): 
                            if "resolution" in ptype:
                                if ("Trk" in n or "TK" in n) and ("jet" not in oname) and ("ele" not in oname):
                                    eq = ("%.1f #times p_{T}^{2} #oplus %.3f #times p_{T} [TeV]" % (p.fit.GetParameter(1), p.fit.GetParameter(0)))
                                else:
                                    eq = ("%.2f #timesp_{T} %+.1f" % (p.fit.GetParameter(1), p.fit.GetParameter(0))).replace("+","+ ").replace(" -"," #minus ")
                            else:
                                eq = ("%.2f #timesp_{T} %+.1f" % (max(0,p.fit.GetParameter(1)), p.fit.GetParameter(0))).replace("+","+ ").replace(" -"," #minus ")
                            leg.AddEntry(p.fit, eq, "L")
                    leg.Draw()
                    out = odir+'/'+oname+pfix+"-"+kind+"_"+ptdef+".png"
                    c1.Print(out)
                    del frame


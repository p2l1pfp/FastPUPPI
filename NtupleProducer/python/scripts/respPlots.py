import os
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kWarning

from math import ceil, floor, sqrt
from array import array

def doRespEta(oname, tree, name, expr, cut):
    return doRespEtaMedian(oname, tree, name, expr, cut)

def quantiles(yswz):
    ys = [y for y in yswz if y > 0]
    if len(ys) < 3: return (0,0,0)
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
        #if yi[i] == 0: continue
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
def doEtaProf(oname, tree, name, expr, cut):
    tree.Draw(expr+":abs(mc_eta)>>"+name+"(20,0,5.0)", cut, "PROF");
    return ROOT.gROOT.FindObject(name);
def doRespEtaProf(oname, tree, name, expr, cut):
    return doEtaProf("("+expr+")/mc_pt");
def ptBins(oname):
    if "jet" in oname: 
        return [20,25,30,35,40,45,50,55,60,70,80,90,100,120,140,160,200,250]
    else:
        return [2.5,5,7.5,10,15,20,25,30,35,40,45,50,55,60,70,80,90,100]
def doPtProf(oname, tree, name, expr, cut):
    ptbins = ptBins(oname)
    ROOT.gROOT.cd()
    prof = ROOT.gROOT.FindObject(name)
    if prof: del prof
    prof = ROOT.TProfile(name,name,len(ptbins),array('f',[0]+ptbins))
    tree.Draw(expr+":mc_pt>>"+name, cut, "PROF");
    return prof;
def doRespPt(oname, tree, name, expr, cut, mcpt="mc_pt", fitopt="WQ0C EX0 ROB=0.95"):
    ptbins = ptBins(oname)
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
        if not median: continue
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
    ('debug-ecal',[
        ("EcalT",  "EcalT$+HGCalTE$", ROOT.kAzure+1,  25, 2.0),
        ("EcalC",  "EcalC$+HGCal3DE$", ROOT.kAzure+2,  21, 1.6),
        ("EcalC2",  "EcalC$+HGCal3DE2$", ROOT.kBlue+1,  21, 1.6),
        ("L1T", "L1RawEcal$", ROOT.kRed+1, 34, 1.5),
        ("L1C", "L1RawEcalC$", ROOT.kViolet+1, 34, 1.5),
    ]),
    ('debug-hcal',[
        ("HcalT",  "Hcal$+HGCalTH$", ROOT.kGreen+3,  25, 2.0),
        ("HcalC",  "Hcal$+HGCal3D$", ROOT.kGreen+1,  21, 1.5),
        ("HcalC2", "Hcal$+HGCal3DH2$", ROOT.kBlue+1,  21, 1.5),
        ("BH",   "HGCalTBH$", ROOT.kGreen+3,  21, 1.5),
        ("L1T", "L1RawCalo$", ROOT.kRed+1, 34, 1.5),
        ("L1C", "L1RawCaloC$", ROOT.kViolet+1, 34, 1.5),
    ]),
    ('debug-calo',[
        ("CalT",  "EcalT$+Hcal$+HGCalTE$+HGCalTH$", ROOT.kGreen+3,  25, 2.0),
        ("CalC",  "EcalC$+Hcal$+HGCal3D$", ROOT.kGreen+1,  21, 1.5),
        ("L1T", "L1RawCalo$", ROOT.kRed+1, 34, 1.5),
        ("L1C", "L1RawCaloC$", ROOT.kViolet+1, 34, 1.5),
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
        ("TK",         "L1TK$",      ROOT.kRed+0, 20, 1.2),
        ("TK #Deltaz", "L1TKV$",     ROOT.kRed+2, 20, 1.2),
        ("PF",         "L1PF$",      ROOT.kOrange+7, 34, 1.2),
        ("Puppi",      "L1Puppi$",   ROOT.kGray+2, 25, 1.4),
    ]),
    ('pfdebug',[
        ("Gen #times Acc", "GenAcc$",           ROOT.kAzure+1,  20, 1.2),
        ("CH #times Acc",  "ChGenAcc$",         ROOT.kAzure+2,  20, 1.2),
        ("NE #times Acc",  "GenAcc$-ChGenAcc$", ROOT.kAzure+3,  20, 1.2),
        ("Calo",       "L1Calo$",               ROOT.kViolet+1, 34, 1.5),
        ("TK",         "L1TK$",                 ROOT.kRed+1,    20, 1.2),
        ("PF",         "L1PF$",                 ROOT.kOrange+7, 34, 1.2),
        ("PFCh",       "L1PFCharged$",          ROOT.kGreen+1,  34, 1.2),
        ("PFNh",       "L1PF$-L1PFCharged$",    ROOT.kGreen+3,  34, 1.2),
    ]),
    ('pfdebug2',[
        ("Gen #times Acc", "GenAcc$",           ROOT.kAzure+1,  20, 1.2),
        ("PH #times Acc",  "PhGenAcc$",         ROOT.kAzure+2,  20, 1.2),
        ("NE #times Acc",  "GenAcc$-ChGenAcc$-PhGenAcc$", ROOT.kAzure+3,  20, 1.2),
        ("Calo",       "L1Calo$",               ROOT.kViolet+1, 34, 1.5),
        ("PF",         "L1PF$",                 ROOT.kOrange+7, 34, 1.2),
        ("PFPh",       "L1PFPhoton$",          ROOT.kGreen+1,  34, 1.2),
        ("PFNh",       "L1PF$-L1PFCharged$-L1PFPhoton$",    ROOT.kGreen+3,  34, 1.2),
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
    parser.add_option("-m", dest="more", default=False, action="store_true", help="make more plots (multiplicity, distance)")
    parser.add_option("-E", "--etaMax", dest="etaMax",  default=5.0, type=float)
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
        etas = [ (0.0,1.3), (1.3,1.7), (1.7,2.5), (2.5,3.0), (3.0,5.0) ]
        for etamin, etamax in etas:
            if etamax > options.etaMax: break
            sels.append(("%s_eta_%02.0f_%02.0f" % (particle,10*etamin,10*etamax), "abs(mc_eta) > %s && abs(mc_eta) < %s && %s" % (etamin,etamax,pdgIdCut)))

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
            cut = cut + " && abs(mc_iso04) < 0.05" # isolated
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
            if options.more:
                ptdefs += [ "ptbest", "mindr025", "n025", "n010" ]
            for ptdef in ptdefs:
                resps = []; resols = [] 
                for name,expr,col,msty,msiz in things:
                    exprptdef = expr.replace("$","_"+ptdef)
                    cutptdef = cut
                    if "mindr025" in ptdef: 
                        if "+" in expr:
                            exprptdef = "min(%s)" % (",".join( exprptdef.split("+") ))
                        cutptdef = cut +" && %s > 0" % expr.replace("$","_n025")
                    if "eta_25" in oname or "eta_30" in oname:
                        if "TK" in expr: continue
                    if ("Puppi" in name or "TKV" in expr) and "PU0" in odir: continue
                    if "Gen" in name and  "jet" not in oname: continue
                    if ptdef.startswith("pt"):
                        if "pt" in oname:
                            prof, pres = doRespEta(oname,tree,name,exprptdef,cutptdef), None
                        else:
                            prof, pres = doRespPt(oname,tree,name,exprptdef,cutptdef)
                    else:
                        if "pt" in oname:
                            prof, pres = doEtaProf(oname,tree,name,exprptdef,cutptdef), None
                        else:
                            prof, pres = doPtProf(oname,tree,name,exprptdef,cutptdef), None
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
                        if "_pt" in oname and ptype == "resolution" and ptdef.startswith("pt"): continue # not implemented
                        print "No ",ptype," plot for ", oname, ptdef 
                        continue
                    if "pt" in oname: 
                        frame = ROOT.TH1F("stk","stk",100,0.0,5.0)
                        frame.GetXaxis().SetTitle("|#eta|")
                        leg = ROOT.TLegend(0.6,0.99,0.95,0.99-0.05*len(things))
                        if ptdef.startswith("pt"):
                            frame.GetYaxis().SetRangeUser(0,2.2)
                            frame.GetYaxis().SetTitle("median p_{T}^{rec}/p_{T}^{gen}")
                        else:
                            frame.GetYaxis().SetTitle("< "+ptdef+" >")
                            frame.GetYaxis().SetRangeUser(0,2*max(h.GetMaximum() for (k,h) in plots))
                    else:
                        frame = ROOT.TH1F("stk","stk",100,0.0,ptBins(oname)[-1])
                        if "resolution" in ptype:
                            frame.GetYaxis().SetTitle("#sigma_{eff}(p_{T}^{corr})/p_{T}^{corr}")
                            frame.GetYaxis().SetRangeUser(0.0,0.8)
                        else:
                            if ptdef.startswith("pt"):
                                frame.GetYaxis().SetTitle("median p_{T}^{rec}/p_{T}^{gen}")
                                frame.GetYaxis().SetRangeUser(0,2.2)
                            else:
                                frame.GetYaxis().SetTitle("< "+ptdef+" >")
                                frame.GetYaxis().SetRangeUser(0,2*max(h.GetMaximum() for (k,h) in plots))
                        frame.GetXaxis().SetTitle("p_{T} (GeV)")
                        leg = ROOT.TLegend(0.2,0.99,0.95,0.99-0.05*(len(plots)+0.5))
                        leg.SetTextSize(0.04);
                        leg.SetNColumns(2)
                    frame.GetYaxis().SetDecimals(True)
                    frame.Draw()
                    line = ROOT.TLine()
                    line.SetLineStyle(7)
                    if "resolution" not in ptype and ptdef.startswith("pt"):
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


import os, sys
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kWarning

from math import ceil, floor, sqrt
from array import array


def quantiles(yswz, zeroSuppress=True):
    ys = [y for y in yswz if y > 0] if zeroSuppress else yswz[:]
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

def doRespEta(oname, tree, name, expr, cut, mcpt="mc_pt", maxEntries=999999999):
    if oname.startswith("null"):
        return doRespEtaProf(oname, tree, name, expr, cut, mcpt=mcpt, maxEntries=maxEntries)
    return doRespEtaMedian(oname, tree, name, expr, cut, mcpt=mcpt, maxEntries=maxEntries)
def doRespEtaMedian(oname, tree, name, expr, cut, mcpt="mc_pt", etabins=25, etamax=5.0, maxEntries=999999999):
    ys = [[] for ieta in xrange(etabins)]
    npoints = tree.Draw("("+expr+")/"+mcpt+":abs(mc_eta)", cut, "", maxEntries);
    #print "Draw(("+expr+")/"+mcpt+":abs(mc_eta), "+cut+") --> ",npoints
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
        (median,lo,hi) = quantiles(ys[ieta], zeroSuppress=not(oname.startswith("null")))
        if not median: continue
        ipoint = ret.GetN()
        ret.Set(ipoint+1)
        ret.SetPoint(ipoint, (ieta+0.5)*etamax/etabins, median)
        #ret.SetPointError(ieta, 0.5*etamax/etabins, 0.5*etamax/etabins, (median-lo),(hi-median))
        ret.SetPointError(ipoint, 0.5*etamax/etabins, 0.5*etamax/etabins, (median-lo)/sqrt(len(ys[ieta])),(hi-median)/sqrt(len(ys[ieta])))
    return ret
def doEtaProf(oname, tree, name, expr, cut, maxEntries=999999999):
    npoints = tree.Draw(expr+":abs(mc_eta)>>"+name+"(50,0,5.0)", cut, "PROF", maxEntries);
    if npoints == 0: return None
    #print "Draw("+expr+":abs(mc_eta), "+cut+") --> ",npoints
    return ROOT.gROOT.FindObject(name);
def doRespEtaProf(oname, tree, name, expr, cut, mcpt="mc_pt", maxEntries=999999999):
    return doEtaProf(oname, tree, name, "("+expr+")/"+mcpt, cut, maxEntries=maxEntries);
def ptBins(oname):
    if "jet" in oname: 
        return [20,25,30,35,40,45,50,55,60,70,80,90,100,120,140,160,200,250,300]
    else:
        return [1,2,3,4,5,7.5,10,12.5,15,17.5,20,25,30,35,40,45,50,55,60,70,80,90,100,120,150,175,200,250]
def doPtProf(oname, tree, name, expr, cut, maxEntries=999999999):
    ptbins = ptBins(oname)
    ROOT.gROOT.cd()
    prof = ROOT.gROOT.FindObject(name)
    if prof: del prof
    prof = ROOT.TProfile(name,name,len(ptbins),array('f',[0]+ptbins))
    tree.Draw(expr+":mc_pt>>"+name, cut, "PROF", maxEntries);
    return prof;
def doRespPt(oname, tree, name, expr, cut, mcpt="mc_pt", xpt="mc_pt", fitopt="WQ0C EX0 ROB=0.95", fitrange=(0,999), maxEntries=999999999, respCorr="simple", fitfunc="lin",ptbins=None,keepGraph=False,zeroIsOk=False,fromGraph=None):
    if not ptbins: ptbins = ptBins(oname)
    ys = [[] for ipt in ptbins]
    if fromGraph == None:
        npoints = tree.Draw("("+expr+")/("+mcpt+"):"+xpt, cut, "", maxEntries);
        if npoints <= 0: return (None,None)
        graph = ROOT.gROOT.FindObject("Graph");
    else:
        graph = fromGraph
        npoints = graph.GetN()
    xi, yi = graph.GetX(), graph.GetY()
    for i in xrange(graph.GetN()):
        if yi[i] == 0 and not zeroIsOk: continue
        for ipt,ptmax in enumerate(ptbins):
            if xi[i] < ptmax:
                ys[ipt].append(yi[i])
                break
    ret = ROOT.TGraphAsymmErrors()
    ret.SetName(name)
    if keepGraph: ret.graph = graph
    retg = ROOT.TGraphAsymmErrors()
    if keepGraph: retg.graph = graph
    retg.SetName(name+"_gauss")
    resols, resolsg = [], []
    for ipt,ptmax in enumerate(ptbins):
        if len(ys[ipt]) == 0: continue
        ptmin = ptbins[ipt-1] if ipt else 0
        ptc   = 0.5*(ptmin+ptmax)
        ptd   = 0.5*(ptmax-ptmin)
        (median,lo,hi) = quantiles(ys[ipt], zeroSuppress=not(zeroIsOk))
        #print "for %s %s at pt ~ %6.1f  reco pt ~ %6.1f  [ %6.1f , %6.1f ]" % (oname,name, ptc, median, lo, hi)
        if not median and not zeroIsOk: continue
        ipoint = ret.GetN()
        ret.Set(ipoint+1)
        ret.SetPoint(ipoint, ptc, median)
        ret.SetPointError(ipoint, ptd, ptd, (median-lo)/sqrt(len(ys[ipt])),(hi-median)/sqrt(len(ys[ipt])))
        ## Now we also try an approximate Gaussian fit
        avg = median; rms2 = (hi - lo);
        for niter in xrange(3):
            truncated = [y for y in ys[ipt] if abs(y-avg) < rms2]
            if len(truncated) <= 2: break
            avg = sum(truncated)/len(truncated)
            rms2 = 2*sqrt(sum((t-avg)**2 for t in truncated)/(len(truncated)-1))
        retg.Set(ipoint+1)
        retg.SetPoint(ipoint, ptc, avg)
        retg.SetPointError(ipoint, ptd, ptd, 0.5*rms2/sqrt(len(ys[ipt])), 0.5*rms2/sqrt(len(ys[ipt])))
        if median > 0.2 or zeroIsOk:
            resols.append((ptc,ptd,(hi-lo)/2,0))
        if avg >  0.2 or zeroIsOk:
            resolsg.append((ptc,ptd,0.5*rms2,0))
    #print "Got ",ret.GetN()
    if ret.GetN() <= 3: return (None,None)
    respformula = "1/x++1"
    if fitfunc == "linsq": respformula = "1/x++1/sqrt(x)++1"
    tf1 = ROOT.TF1(name+"_f1", respformula, 0, ret.GetX()[ret.GetN()-1]+ret.GetErrorXhigh(ret.GetN()-1) )
    #tf1 = ROOT.TF1(name+"_f1","1/x++1++1/(x*x)", ret.GetX()[0]-ret.GetErrorXlow(0), ret.GetX()[ret.GetN()-1]+ret.GetErrorXhigh(ret.GetN()-1) )
    ret.Fit(tf1, fitopt, "", fitrange[0], fitrange[1])
    ret.fit = tf1
    tf1g = ROOT.TF1(name+"_f1g",respformula, 0, ret.GetX()[ret.GetN()-1]+ret.GetErrorXhigh(ret.GetN()-1) )
    retg.Fit(tf1g, fitopt, "", fitrange[0], fitrange[1])
    retg.fit = tf1
    ret.gaus = retg
    retc = ROOT.TGraphAsymmErrors()
    retc.SetName(name+"_corr")
    retcg = ROOT.TGraphAsymmErrors()
    retcg.SetName(name+"_corr_gauss")
    #print oname, expr
    if respCorr in ("gen","exact","fit"):
        if respCorr == "exact":
            invret = ROOT.TGraph(ret.GetN())
            for i in xrange(ret.GetN()):
                ptgen = ret.GetX()[i]
                ptrec = ret.GetY()[i]*ptgen
                invret.SetPoint(i, ptrec, ptgen)
                #if i % 2 == 1: print "JEC: %5.1f -> %5.1f" % (ptrec, ptgen)
            invret.Sort()
            ret.inv = invret
        elif respCorr == "fit":
            if respformula == "1/x++1": 
                scale, offs = tf1.GetParameter(1), tf1.GetParameter(0)
            else:
                raise RuntimeError("Only linear fits supported") 
        ycs = [[] for ipt in ptbins]
        #iprint = 3
        for i in xrange(graph.GetN()):
            if yi[i] == 0: continue
            if respCorr == "gen":
                yic = yi[i] / ret.Eval(xi[i]);
            elif respCorr == "exact":
                yic = invret.Eval(yi[i]*xi[i]) / xi[i];
                #if i == iprint:
                #    print "%8d   xi %5.1f  yi %7.3f   ptrec %5.1f   ptcorr %5.1f  yic %7.3f" % (i, xi[i], yi[i], yi[i]*xi[i], invret.Eval(yi[i]*xi[i]), yic)
                #    iprint *= 3
            elif respCorr == "fit":
                if respformula == "1/x++1": 
                    yic = (yi[i]*xi[i] - offs) / scale / xi[i];
            for ipt,ptmax in enumerate(ptbins):
                if xi[i] < ptmax:
                    ycs[ipt].append(yic)
                    break
        resols, resolsg = [], []
        for ipt,ptmax in enumerate(ptbins):
            if len(ycs[ipt]) == 0: continue
            ptmin = ptbins[ipt-1] if ipt else 0
            ptc   = 0.5*(ptmin+ptmax)
            ptd   = 0.5*(ptmax-ptmin)
            (median,lo,hi) = quantiles(ycs[ipt])
            if not median: continue
            ipoint = retc.GetN()
            retc.Set(ipoint+1)
            retc.SetPoint(ipoint, ptc, median)
            retc.SetPointError(ipoint, ptd, ptd, (median-lo)/sqrt(len(ycs[ipt])),(hi-median)/sqrt(len(ycs[ipt])))
            ## Now we also try an approximate Gaussian fit
            avg = median; rms2 = (hi - lo);
            for niter in xrange(3):
                truncated = [y for y in ycs[ipt] if abs(y-avg) < rms2]
                if len(truncated) <= 2: break
                avg = sum(truncated)/len(truncated)
                rms2 = 2*sqrt(sum((t-avg)**2 for t in truncated)/(len(truncated)-1))
            retcg.Set(ipoint+1)
            retcg.SetPoint(ipoint, ptc, avg)
            retcg.SetPointError(ipoint, ptd, ptd, 0.5*rms2/sqrt(len(ycs[ipt])), 0.5*rms2/sqrt(len(ycs[ipt])))
            if median > 0.2:
                resols.append((ptc,ptd,(hi-lo)/2,0))
                #print "for %s %s at pt ~ %6.1f  reco pt ~ %6.1f  [ %6.1f , %6.1f ]" % (oname,name, ptc, median, lo, hi)
            if avg >  0.2:
                resolsg.append((ptc,ptd,0.5*rms2,0))
        ret.corr = retc
        retg.corr = retcg
        scale = 1
    elif respCorr == "none":
        scale = 1
    elif respCorr == "simple":
        ## Now we do the approximate response corrected resolution
        scale = ret.fit.GetParameter(1)
    else:
        raise RuntimeError("Unsuppored respCorr = "+respCorr)
    if scale < 0.1 or len(resols) < 3: return (ret,None)
    resol = ROOT.TGraphAsymmErrors(len(resols))
    resol.SetName(name+"_res")
    for ipoint,(ptc,ptd,res,err) in enumerate(resols):
        if respCorr == "divide": scale = ret.Eval(ptc)
        resol.SetPoint(ipoint, ptc, res/scale)
        resol.SetPointError(ipoint, ptd, ptd, 0, 0)
    if ("Trk" in name or "TK" in name) and "jet" not in oname and "ele" not in oname:
        rtf1 = ROOT.TF1(name+"_rf1","sqrt([0]*[0]+[1]*[1]*(x/1000)*(x/1000))", 0, ret.GetX()[ret.GetN()-1]+ret.GetErrorXhigh(ret.GetN()-1) )
        rtf1.SetParameters(0.01, 1)
    else:
        rtf1 = ROOT.TF1(name+"_rf1","1/x++1", 0, ret.GetX()[ret.GetN()-1]+ret.GetErrorXhigh(ret.GetN()-1) )
    resol.Fit(rtf1, fitopt+"W", "", fitrange[0], fitrange[1])
    resol.fit = rtf1
    if len(resolsg) > 3:
        resolg = ROOT.TGraphAsymmErrors(len(resolsg))
        resolg.SetName(name+"_gauss_res")
        for ipoint,(ptc,ptd,res,err) in enumerate(resolsg):
            if respCorr == "divide": scale = ret.Eval(ptc)
            resolg.SetPoint(ipoint, ptc, res/scale)
            resolg.SetPointError(ipoint, ptd, ptd, 0, 0)
        resol.gaus = resolg
        rtf1g = rtf1.Clone(name+"_rf1g")
        resolg.Fit(rtf1g, fitopt+"W", "", fitrange[0], fitrange[1])
        resolg.fit = rtf1g
    return (ret,resol)

whats = [
    ('debug-ecal',[
        ("RawOld",  "L1RawBarrelEcalOldL1EG$", ROOT.kGreen+2,  25, 2.0),
        ("Raw",  "L1RawBarrelEcal$",           ROOT.kGreen+2,  21, 1.6),
        ("CorrOld",  "L1BarrelEcalOldL1EG$",   ROOT.kAzure+2,  24, 1.3),
        ("Corr",  "L1BarrelEcal$",             ROOT.kAzure+2,  20, 1.0),
    ]),
    ('debug-hf',[
        ("Cells",  "L1RawHFCells$",  ROOT.kGreen+2,  20, 1.1),
        ("Raw",  "L1RawHFCalo$",     ROOT.kAzure+1,  25, 2.0),
        ("Corr",  "L1HFCalo$",       ROOT.kAzure+2,  21, 1.6),
    ]),
    ('debug-hcal',[
        ("Raw_OldTow",   "L1RawBarrelCaloOldTowers$",        ROOT.kGreen+1,  25, 2.0),
        ("Raw_NewTow",   "L1RawBarrelCalo$",   ROOT.kAzure+2,  21, 1.6),
        ("Raw_UncOnly",   "L1RawBarrelCaloUnclust$",   ROOT.kRed-7,    21, 1.4),
        ("Corr_Old",  "L1BarrelCaloOldTowers$",           ROOT.kGreen+1,  24, 1.3),
        ("Corr_New",  "L1BarrelCalo$",      ROOT.kAzure+2,  20, 1.0),
        ("RawEM_Old", "L1RawBarrelCaloEMOldTowers$",      ROOT.kGreen+3,  20, 0.9),
        ("RawEM_New", "L1RawBarrelCaloEM$", ROOT.kViolet+1, 20, 0.9),
        ("RawEM_Unc", "L1RawBarrelCaloEMUnclust$", ROOT.kRed+2,    20, 0.9),
    ]),
    ('debug-hgc',[
        ("Raw",     "L1RawHGCal$",   ROOT.kGreen+2,   25, 2.0),
        ("RawEM",   "L1RawHGCalEM$", ROOT.kRed+1,     21, 1.5),
        ("Corr",    "L1HGCal$",      ROOT.kAzure+1,   20, 1.3),
        ("CorrEM",  "L1HGCalEM$",    ROOT.kViolet+1,  20, 0.9),
    ]),
    ('l1pfold',[
        ("Gen #times Acc",        "GenAcc$",    ROOT.kAzure+1,  20, 1.2),
        ("Raw Ecal",   "L1OldRawEcal$", ROOT.kGreen+3,  21, 1.7),
        ("Raw Calo",   "L1OldRawCalo$", ROOT.kViolet-4,  21, 1.7),
        ("Ecal",       "L1OldEcal$",    ROOT.kGreen+1,  21, 1.7),
        ("Calo",       "L1OldCalo$",    ROOT.kViolet+2, 34, 1.5),
        ("TK",         "L1TK$",      ROOT.kRed+0, 20, 1.2),
        ("PF",         "L1OldPF$",      ROOT.kOrange+7, 20, 1.2),
    ]),
    ('l1pfoldpu',[
        ("Gen #times Acc",        "GenAcc$",    ROOT.kAzure+1,  20, 1.2),
        ("Calo",       "L1OldCalo$",    ROOT.kViolet+1, 21, 1.5),
        ("TK #Deltaz", "L1TKV5$",       ROOT.kRed+1, 34, 1.2),
        ("PF",         "L1OldPF$",      ROOT.kOrange+7, 21, 1.4),
        ("Puppi",      "L1OldPuppi$",   ROOT.kGreen+1, 20, 1.5),
        ("Puppi4MET",  "L1OldPuppiForMET$",   ROOT.kGreen+3, 20, 1.1),
    ]),
    ('l1pfnew',[
        ("Gen #times Acc",        "GenAcc$",    ROOT.kAzure+1,  20, 1.2),
        ("Raw Calo",   "L1RawBarrelCalo$+L1RawHGCal$", ROOT.kViolet-4,  21, 1.7),
        ("Ecal",       "L1BarrelEcal$+L1RawHGCalEM$",    ROOT.kGreen+1,  21, 1.7),
        ("Calo",       "L1Calo$",    ROOT.kViolet+2, 34, 1.5),
        ("TK",         "L1TK$",      ROOT.kRed+0, 20, 1.2),
        ("PF",         "L1PF$",      ROOT.kOrange+7, 20, 1.2),
        ("Puppi",      "L1Puppi$",   ROOT.kGray+2, 20, 1.2),
    ]),
    ('l1pfnewpu',[
        ("Gen #times Acc",        "GenAcc$",    ROOT.kAzure+1,  20, 1.2),
        ("Calo",       "L1Calo$",    ROOT.kViolet+1, 34, 1.5),
        ("TK #Deltaz", "L1TKV5$",    ROOT.kRed+1, 20, 1.2),
        ("PF",         "L1PF$",      ROOT.kOrange+7, 20, 1.2),
        ("Puppi",      "L1Puppi$",   ROOT.kGreen+1, 20, 1.5),
        ("Puppi4MET",  "L1PuppiForMET$",   ROOT.kGreen+3, 20, 1.1),
    ]),
    ('l1pfcomp',[
        ("Ecal",      "L1BarrelEcal$+L1HGCalEM$", ROOT.kGreen+1,  21, 1.7),
        ("OldCalo",   "L1OldCalo$", ROOT.kViolet+1, 25, 1.6),
        ("NewCalo",   "L1Calo$",    ROOT.kViolet+2, 21, 1.1),
        ("OldPF",     "L1OldPF$",   ROOT.kOrange-3, 24, 1.5),
        ("NewPF",     "L1PF$",      ROOT.kRed+1,    20, 1.0),
    ]),
    ('l1pfcomppu',[
        ("OldPuppi",       "L1OldPuppi$",       ROOT.kGreen+0,  24, 1.4),
        ("NewPuppi",       "L1Puppi$",          ROOT.kGreen+2,  20, 1.1),
        ("OldPuppiForMET", "L1OldPuppiForMET$", ROOT.kAzure+10, 24, 1.4),
        ("NewPuppiForMET", "L1PuppiForMET$",    ROOT.kAzure+2,  20, 1.1),
    ]),
    ('pfdebug',[
        ("Gen #times Acc", "GenAcc$",           ROOT.kGray+1,  21, 1.2),
        ("PF",         "L1PF$",                 ROOT.kBlack,    20, 1.0),
        ("Gen_Charged",  "ChGenAcc$",         ROOT.kAzure+10,  21, 1.5),
        ("PF_Charged",       "L1PFCharged$",          ROOT.kBlue+2,   20, 1.3),
        ("Gen_Neutral",  "GenAcc$-ChGenAcc$", ROOT.kGreen+0,  21, 1.2),
        ("PF_Neutral",       "L1PFNeutral$",    ROOT.kGreen+2,  20, 1.0),
        ("RawTK",      "RawTK$",                ROOT.kRed+2,    34, 1.0),
        ("TK",         "L1TK$",                 ROOT.kRed+0,    34, 1.0),
        ("Calo",       "L1Calo$",               ROOT.kViolet+1, 34, 1.0),
    ]),
    ('pfdebug1',[
        ("Gen #times Acc", "GenAcc$",           ROOT.kGray+1,  21, 1.2),
        ("PF",         "L1PF$",                 ROOT.kBlack,    20, 1.0),
        ("Gen_Charged",  "ChGenAcc$",         ROOT.kAzure+10,  21, 1.5),
        ("PF_Charged",   "L1PFCharged$",          ROOT.kBlue+2,   20, 1.3),
        ("RawTK",      "RawTK$",                ROOT.kRed+2,    34, 1.0),
        ("TK",         "L1TK$",                 ROOT.kRed+0,    34, 1.0),
    ]),
    ('pfdebug2',[
        ("Gen #times Acc", "GenAcc$",           ROOT.kGray+1,  21, 1.2),
        ("PF",         "L1PF$",                 ROOT.kBlack,    20, 1.0),
        ("Gen_Photon",  "PhGenAcc$",         ROOT.kGreen+0,  21, 1.2),
        ("PF_Photon",       "L1PFPhoton$",           ROOT.kGreen+2,  20, 1.0),
        ("Gen_NeutralHad",  "GenAcc$-ChGenAcc$-PhGenAcc$", ROOT.kRed-4,  21, 1.2),
        ("PF_NeutralHad",   "L1PFNeutralHad$",  ROOT.kRed+2,  20, 1.0),
    ]),
    ('puppidebug',[
        ("Gen #times Acc", "GenAcc$",         ROOT.kGray+1,  21, 1.2),
        ("PF",             "L1PF$",           ROOT.kGray+2,   24, 1.0),
        ("Puppi",          "L1Puppi$",        ROOT.kBlack,    20, 1.0),
        ("Gen_Charged",  "ChGenAcc$",         ROOT.kAzure+10,  21, 1.5),
        ("PF_Charged",       "L1PFCharged$",          ROOT.kBlue+2,   24, 1.3),
        ("Puppi_Charged",       "L1PuppiCharged$",    ROOT.kBlue+2,   20, 1.3),
        ("Gen_Neutral",    "GenAcc$-ChGenAcc$", ROOT.kGreen+0,  21, 1.2),
        ("PF_Neutral",     "L1PFNeutral$",        ROOT.kGreen+2,  24, 1.0),
        ("Puppi_Neutral",  "L1PuppiNeutral$",  ROOT.kGreen+2,  20, 1.0),
        ("TK",         "L1TK$",                 ROOT.kRed+0,    34, 1.0),
        ("TKV",        "L1TKV$",                ROOT.kRed+2,    34, 1.0),
    ]),
     ('puppidebug2',[
        ("Gen #times Acc", "GenAcc$",           ROOT.kGray+1,  21, 1.2),
        ("PF",         "L1PF$",                 ROOT.kBlack,    20, 1.0),
        ("Gen_Photon",  "PhGenAcc$",         ROOT.kGreen+0,  21, 1.2),
        ("PF_Photon",       "L1PFPhoton$",           ROOT.kGreen+2,  20, 1.0),
        ("Gen_NeutralHad",  "GenAcc$-ChGenAcc$-PhGenAcc$", ROOT.kRed-4,  21, 1.2),
        ("PF_NeutralHad",   "L1PFNeutralHad$",  ROOT.kRed+2,  20, 1.0),
    ]),

    ('pfcomp4',[
        ("Global",    "L1PFGlobal$",   ROOT.kGreen+1,  20, 1.2),
        ("Regional",  "L1PFRegional$", ROOT.kAzure+1,  24, 1.2),
        ("Simpler",   "L1PFSimpler$",  ROOT.kViolet+1, 21, 1.1),
        ("Firmware",  "L1PF$",         ROOT.kOrange+7, 34, 1.2),
    ]),



]

if __name__ == "__main__":
    from optparse import OptionParser
    parser = OptionParser("%(prog) infile [ src [ dst ] ]")
    parser.add_option("-w", dest="what",     default=None, help="Choose set (inputs, l1pf, ...)")
    parser.add_option("-p", dest="particle", default=None, help="Choose particle (electron, ...)")
    parser.add_option("--ptdef", dest="ptdef", default=None, help="Pt definition")
    parser.add_option("--ptmax", dest="ptmax", nargs=1, default=9999, type="float")
    parser.add_option("--ptminFit", dest="ptminFit", default=0, type="float", help="Pt definition")
    parser.add_option("--cut", dest="extracut", default=None, help="Extra cut", nargs=2)
    parser.add_option("--label", dest="extralabel", default=None, help="Extra label")
    parser.add_option("-m","--more", dest="more", default=False, action="store_true", help="make more plots (multiplicity, distance)")
    parser.add_option("-g","--gauss", dest="gauss", default=False, action="store_true", help="make also gaussian estimates")
    parser.add_option("-G","--gauss-only", dest="gaussOnly", default=False, action="store_true", help="make only gaussian estimates")
    parser.add_option("--fit", dest="fit", default="lin", type="string", help="fit: lin, linsq, none")
    parser.add_option("--no-fit", dest="fit", action="store_const", const="none", help="skip fits")
    parser.add_option("--no-resol", dest="noResol", default=False, action="store_true", help="skip resolution plots")
    parser.add_option("--fitrange", dest="fitrange", nargs=2, default=(0,999), type="float")
    parser.add_option("--no-eta", dest="noEtaPlot", default=False, action="store_true", help="skip resolution plots")
    parser.add_option("--corr-resol", dest="respCorr", default="exact", help="response correction for resolution: exact (compute & apply JECs as function of pt reco, then get resolution), gen (JECs vs pt gen), divide (divide by response), simple (divide by response at infinity)")
    parser.add_option("-E", "--etaMax", dest="etaMax",  default=5.0, type=float)
    parser.add_option("--ymax", dest="yMax",  default=1.6, type=float)
    parser.add_option("--yrange", "--yr", dest="yRange",  default=None, type=float, nargs=2, help="Y axis range (for profile plots)")
    parser.add_option("--ymaxRes", dest="yMaxRes",  default=0.0, type=float)
    parser.add_option("--eta", dest="eta", nargs=2, default=None, type="float")
    parser.add_option("--eff-threshold", dest="efft", nargs=2, default=(0.2, 0.5), type="float", help="Minimum reconstructed pT for efficiency plots. Takes two argument A, B to get reco_pt > A * mc_pt + B")
    parser.add_option("-M", "--max-entries", dest="maxEntries",  default=999999999, type=int)
    parser.add_option("--mc", "--mc", dest="mcpt",  default="mc_pt")
    parser.add_option("--xpt", "--xpt", dest="xpt",  default=("mc_pt","p_{T}^{gen}"), nargs=2)
    options, args = parser.parse_args()
    selparticles = options.particle.split(",") if options.particle else []

    sels = []; 
    for (particle, pdgIdCut, minPt, maxEta) in [ 
            ("pion", "abs(mc_id) == 211", 10, 5),
            ("pizero", "abs(mc_id) == 111", 10, 5),
            ("klong", "abs(mc_id) == 130", 10, 5),
            ("pimix", "(abs(mc_id) == 211 || abs(mc_id) == 111)", 10, 5),
            ("mixmix", "(abs(mc_id) == 211 || abs(mc_id) == 22)", 10, 5),
            ("mix", "(abs(mc_id) == 211 || abs(mc_id) == 111 || abs(mc_id) == 22 || abs(mc_id) == 130)", 2, 5),
            ("photon", "abs(mc_id) == 22", 10, 5),
            ("electron", "abs(mc_id) == 11", 10, 5),
            ("muon", "abs(mc_id) == 13", 10, 5),
            ("tau", "(abs(mc_id) == 15 || abs(mc_id) == 211)", 20, 5),
            ("jet", "abs(mc_id) == 0", 30, 5),
            ("null", "abs(mc_id) == 999", 0, 5)
        ]:
        if options.particle and (particle not in selparticles): 
            continue
        if options.extracut:
            particle += "_"+options.extracut[0]
            pdgIdCut += " && ("+options.extracut[1]+")"
        if options.extralabel:
            particle += "_"+options.extralabel
        if not options.noEtaPlot:
            sels.append(("%s_pt_%02d_inf" % (particle, minPt), "mc_pt > %g && %s" % (minPt, pdgIdCut)))
        if "null" in particle: continue; # not point in profiling random cones vs pt
        etas = [ (0.0,1.3), (1.3,1.7), (1.7,2.5), (2.5,3.0), (3.0,5.0) ]
        if options.eta: etas = [ options.eta ]
        for etamin, etamax in etas:
            if etamax > options.etaMax: break
            sels.append(("%s_eta_%02.0f_%02.0f" % (particle,10*etamin,10*etamax), "abs(mc_eta) > %s && abs(mc_eta) < %s && %s" % (etamin,etamax,pdgIdCut)))

    tree = ROOT.TChain("ntuple/tree")
    for a in args[:]:
        if a.endswith(".root"): 
            tree.Add(a)
            args.remove(a)
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
        isNeutrino = oname.startswith("null")
        if isNeutrino and options.mcpt == "mc_pt": options.mcpt = "1"
        if "jet" not in oname:
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
            elif isNeutrino:
                ptdefs = [ "pthighest", "pt", "ptbest", "pt02" ]
            else:
                ptdefs = [ "pt02", "ptbest" ]   
            if options.more:
                ptdefs += [ "ptbest", "mindr025", "n025", "n010" ]
            if options.ptdef: ptdefs = options.ptdef.split(",")
            for ptdef in ptdefs:
                resps = []; resols = []; gresps = []; gresols = []; cresps = []; cgresps = []
                for name,expr,col,msty,msiz in things:
                    exprptdef = expr.replace("$","_"+ptdef)
                    if ptdef.startswith("eff"):
                        mypt = expr.replace("$","_"+ptdef.replace("eff","pt"))
                        exprptdef = "%s > %g*%s + %g" % (mypt, options.efft[0], options.mcpt, options.efft[1])
                    elif ptdef == "psplit":
                        mypt = expr.replace("$","_pt")
                        exprptdef = "(%s02 - %sbest) > %g*%s + %g" % (mypt, mypt, options.efft[0], options.mcpt, options.efft[1])
                    cutptdef = cut
                    if "mindr025" in ptdef: 
                        if "+" in expr:
                            exprptdef = "min(%s)" % (",".join( exprptdef.split("+") ))
                        cutptdef = cut +" && %s > 0" % expr.replace("$","_n025")
                    if "eta_25" in oname or "eta_30" in oname:
                        if "TK" in expr: continue
                    if "photon" in oname and "TK" in expr: continue
                    if ("Puppi" in name) and "PU0" in odir: continue
                    if ("TKV" in expr)   and "PU0" in odir: continue
                    if "Gen" in name and  "jet" not in oname: continue
                    if "#times Acc" in name and  "jet" not in oname: continue
                    if ptdef.startswith("pt"):
                        if "pt" in oname:
                            prof, pres = doRespEta(oname,tree,name,exprptdef,cutptdef,mcpt=options.mcpt,maxEntries=options.maxEntries), None
                        else:
                            prof, pres = doRespPt(oname,tree,name,exprptdef,cutptdef,mcpt=options.mcpt,xpt=options.xpt[0],maxEntries=options.maxEntries,respCorr=options.respCorr,fitfunc=options.fit,fitrange=options.fitrange)
                    else:
                        if "pt" in oname:
                            prof, pres = doEtaProf(oname,tree,name,exprptdef,cutptdef,maxEntries=options.maxEntries), None
                        else:
                            prof, pres = doPtProf(oname,tree,name,exprptdef,cutptdef,maxEntries=options.maxEntries), None
                    #print oname, kind, ptdef, name, prof, pres
                    allplots = [(prof,resps), (pres,resols)]
                    if options.gauss or options.gaussOnly:
                        if options.gaussOnly and prof and getattr(prof,'gaus',None) != None:
                            allplots = []
                        if prof: allplots.append( (getattr(prof,'gaus',None),gresps) )
                        if pres: allplots.append( (getattr(pres,'gaus',None),gresols) )
                    #if options.respCorrExact:
                    if prof: allplots.append( (getattr(prof,'corr',None), cresps) )
                    for (p,ps) in allplots:
                        if not p: continue
                        if p.InheritsFrom("TGraph") and p.GetN() == 0: 
                            print "no points in non-null plot "+p.GetName()
                            continue
                        if getattr(p,'fit',None):
                            p.fit.SetLineWidth(2); p.fit.SetLineColor(col)
                            p.fit.SetRange(options.ptminFit, p.fit.GetXmax())
                        p.SetLineWidth(3); p.SetLineColor(col);  p.SetMarkerColor(col)
                        p.SetMarkerStyle(msty); p.SetMarkerSize(msiz)
                        ps.append((name,p))
                allplots = [ (resps,"response",""),(resols,"resolution","_res") ]
                if options.gauss or options.gaussOnly:
                    allplots += [ (gresps, "response_gauss", "_gauss") ]
                    allplots += [ (gresols, "resolution_gauss", "_res_gauss") ]
                #allplots += [ (cresps, "response_corr", "_corr") ]
                for plots,ptype,pfix in allplots:
                    if "resolution" in ptype  and options.noResol: continue
                    if not plots: 
                        if "_pt" in oname and ("resolution" in ptype) and ptdef.startswith("pt"): continue # not implemented
                        print "No ",ptype," plot for ", oname, ptdef 
                        continue
                    if "pt" in oname:
                        if plots[0][1].InheritsFrom("TGraph") and not isNeutrino:
                            minEta = min(min(p.GetX()[i] - p.GetErrorXlow(i)  for i in xrange(p.GetN())) for (n,p) in plots)
                            maxEta = max(max(p.GetX()[i] + p.GetErrorXhigh(i) for i in xrange(p.GetN())) for (n,p) in plots)
                            if abs(minEta-0) < 0.05: minEta = 0
                            if abs(maxEta-5) < 0.05: maxEta = 5
                        else:
                            (minEta,maxEta) = plots[0][1].GetXaxis().GetXmin(), plots[0][1].GetXaxis().GetXmax()
                            if options.eta: minEta, maxEta = options.eta
                        frame = ROOT.TH1F("stk","stk",100,minEta,maxEta)
                        frame.GetXaxis().SetTitle("|#eta|")
                        ymax = options.yMax
                        leg = ROOT.TLegend(0.2,0.99,0.95,0.99-0.025*len(things))
                        leg.SetTextSize(0.04);
                        leg.SetNColumns(2)
                        if ptdef.startswith("pt"):
                            frame.GetYaxis().SetRangeUser(0,options.yMax)
                            frame.GetYaxis().SetTitle("median p_{T}^{rec}/p_{T}^{gen}")
                            if isNeutrino:
                                if "--ymax" in sys.argv:
                                    ymax = options.yMax
                                elif "TGraph" in plots[0][1].ClassName():
                                    ymax = 1.3*max(max(p.GetY(i) for i in xrange(p.GetN())) for (n,p) in plots)
                                else:
                                    ymax = 1.3*max(p.GetMaximum() for (n,p) in plots)
                                frame.GetYaxis().SetRangeUser(0,ymax)
                                frame.GetYaxis().SetTitle("mean p_{T}^{rec} (GeV)") 
                        else:
                            frame.GetYaxis().SetTitle("< "+ptdef+" >")
                            if options.yRange:
                                frame.GetYaxis().SetRangeUser(options.yRange[0], options.yRange[1])
                            else:
                                frame.GetYaxis().SetRangeUser(0,2*max(h.GetMaximum() for (k,h) in plots))
                        hasfitlegend = False
                    else:
                        frame = ROOT.TH1F("stk","stk",100,0.0,min(ptBins(oname)[-1],options.ptmax))
                        if "resolution" in ptype:
                            hasfitlegend = (options.fit != "none")
                            frame.GetYaxis().SetTitle("#sigma(p_{T}^{corr})/p_{T}^{corr}")
                            ymax = options.yMaxRes if options.yMaxRes else  (1.2 if "PU200" in odir else 0.8) 
                            frame.GetYaxis().SetRangeUser(0.0, ymax)
                        else:
                            if ptdef.startswith("pt"):
                                hasfitlegend = (options.fit == "lin")
                                frame.GetYaxis().SetTitle("median p_{T}^{rec}/p_{T}^{gen}")
                                frame.GetYaxis().SetRangeUser(0, options.yMax)
                            else:
                                hasfitlegend = False
                                frame.GetYaxis().SetTitle("< "+ptdef+" >")
                                if options.yRange:
                                    frame.GetYaxis().SetRangeUser(options.yRange[0], options.yRange[1])
                                else:
                                    frame.GetYaxis().SetRangeUser(0,2*max(h.GetMaximum() for (k,h) in plots))
                        frame.GetXaxis().SetTitle(options.xpt[1]+" (GeV)")
                        legsize = 0.042*(len(plots)+0.5)
                        if not hasfitlegend: legsize /= 2;
                        leg = ROOT.TLegend(0.2,0.99,0.95,0.99-legsize)
                        leg.SetTextSize(0.035);
                        leg.SetNColumns(2)
                    frame.GetYaxis().SetDecimals(True)
                    frame.Draw()
                    line = ROOT.TLine()
                    line.SetLineStyle(7)
                    if "resolution" not in ptype and ptdef.startswith("pt") and options.mcpt != "1":
                        if "pt" in oname: 
                            line.DrawLine(minEta,1,maxEta,1)
                        else:
                            line.DrawLine(0.0,1,frame.GetXaxis().GetXmax(),1)
                    if "pt" in oname: 
                        for etaLine in 1.5, 2.5, 3.0:
                            if minEta < etaLine and etaLine < maxEta: 
                                line.DrawLine(etaLine,0,etaLine,ymax)
                    for n,p in plots: 
                        p.Draw("P SAME" if "TGraph" in p.ClassName() else "SAME")
                        if hasattr(p,'fit') and (options.fit != "none"): p.fit.Draw("SAME")
                    for n,p in plots: 
                        leg.AddEntry(p, n.replace("_"," "), "LP")
                        if hasfitlegend: 
                            if "resolution" in ptype:
                                if ("Trk" in n or "TK" in n) and ("jet" not in oname) and ("ele" not in oname):
                                    eq = ("%.1f #times p_{T}^{2} #oplus %.3f #times p_{T} [TeV]" % (p.fit.GetParameter(1), p.fit.GetParameter(0)))
                                else:
                                    eq = ("%.3f #timesp_{T} %+.1f" % (p.fit.GetParameter(1), p.fit.GetParameter(0))).replace("+","+ ").replace(" -"," #minus ")
                            else:
                                eq = ("%.2f #timesp_{T} %+.1f" % (max(0,p.fit.GetParameter(1)), p.fit.GetParameter(0))).replace("+","+ ").replace(" -"," #minus ")
                            leg.AddEntry(p.fit, eq, "L")
                    leg.Draw()
                    out = odir+'/'+oname+pfix+"-"+kind+"_"+ptdef+".png"
                    c1.Print(out)
                    del frame


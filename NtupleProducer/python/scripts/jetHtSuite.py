import re
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kWarning

from array import array
from math import pow, sin, cos, hypot, sqrt
import itertools
from FastPUPPI.NtupleProducer.scripts.makeJecs import _progress
from FastPUPPI.NtupleProducer.plotTemplate import plotTemplate


def calcHT(jets):
    return sum(j[0] for j in jets) if jets else 0
def calcMHT(jets):
    return hypot(sum(j[0]*sin(j[2]) for j in jets),sum(j[0]*cos(j[2]) for j in jets))  if jets else 0
def calcMJJ(jets):
    if len(jets) <= 1: return 0
    lvclass = ROOT.Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<double>")
    jp4s = [ lvclass(j[0],j[1],j[2],0) for j in jets ]
    return max((j1+j2).M() for (j1,j2) in itertools.combinations(jp4s,2))

class CalcJ:
    def __init__(self,index):
        self._index = index
    def __call__(self,jets):
        if len(jets) <= self._index: return 0
        jsort = sorted(jets, key = lambda j : -j[0])
        return jsort[self._index][0]
class CalcJ2_MJJcut:
    def __init__(self,mjj):
        self._mjj = mjj
    def __call__(self,jets):
        if len(jets) <= 1: return 0
        lvclass = ROOT.Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<double>")
        jp4s = [ lvclass(j[0],j[1],j[2],0) for j in jets ]
        jjs = [ jj for jj in itertools.combinations(jp4s,2) if (jj[0]+jj[1]).M() > self._mjj ]
        if len(jjs) <= 0: return 0
        return max(min(jj[0].Pt(), jj[1].Pt()) for jj in jjs)

def makeCalc(what):
    if what == "ht": 
        return calcHT
    if what == "mht": 
        return calcMHT
    if what == "mjj": 
        return calcMJJ
    if re.match(r"jet\d+$", options.var): 
        return CalcJ(int(options.var.replace("jet",""))-1)
    if what.startswith("ptj-mjj"): 
        return CalcJ2_MJJcut(float(what.replace("ptj-mjj","")))

def makeGenArray(tree, what, ptCut, etaCut, _cache={}):
    _key = (id(tree),what,int(ptCut*100),int(etaCut*1000))
    if _key in _cache: return _cache[_key]
    if what == "metmht":
        met = makeGenArray(tree, "met",     0,    5.0, _cache=_cache)
        mht = makeGenArray(tree, "mht", ptCut, etaCut, _cache=_cache)
        ret = map(min, zip(met,mht))
        _cache[_key] = ret
        return ret
    if "met" in what:
        ret = makeGenMETArray(tree, what, etaCut)
        _cache[_key] = ret
        return ret
    progress = _progress("Reading GenJets ...")
    calc = makeCalc(what)
    ret = []
    tree.SetBranchStatus("*",0);
    tree.SetBranchStatus("nGenJets",1);
    tree.SetBranchStatus("GenJets_*",1);
    for i in xrange(tree.GetEntries()):
        tree.GetEntry(i)
        pt,eta,phi = tree.GenJets_pt, tree.GenJets_eta, tree.GenJets_phi
        jets = [ (pt[j],eta[j],phi[j]) for j in xrange(tree.nGenJets) if pt[j] > ptCut and abs(eta[j]) < etaCut ]
        ret.append(calc(jets))
    _cache[_key] = ret
    progress.done("done, %d entries" % len(ret))
    return ret
def makeGenMETArray(tree, what, etaCut):
    if   etaCut <= 1.5: post = "MetBarrel_pt" 
    elif etaCut <= 2.4: post = "MetCentral_pt"
    else:               post = "Met_pt"
    progress = _progress("Reading gen"+post+" ...")
    tree.SetBranchStatus("*",0);
    tree.SetBranchStatus("gen"+post,1);
    ret = []
    for i in xrange(tree.GetEntries()):
        tree.GetEntry(i)
        ret.append(getattr(tree,"gen"+post))
    progress.done("done, %d entries" % len(ret))
    return ret
def makeCorrArray(tree, what, obj, ptCorrCut, etaCut, corr, _cache={}):
    _key = (id(tree),what,obj,int(ptCorrCut*100),int(etaCut*1000))
    if _key in _cache: return _cache[_key]
    if what == "metmht":
        met = makeCorrArray(tree, "met", obj, ptCorrCut,    5.0, corr, _cache=_cache)
        mht = makeCorrArray(tree, "mht", obj, ptCorrCut, etaCut, corr, _cache=_cache)
        ret = map(min, zip(met,mht))
        _cache[_key] = ret
        return ret
    if "met" in what:
        ret = makeRecoMETArray(tree, what, obj, etaCut)
        _cache[_key] = ret
        return ret
    calc = makeCalc(what)
    ret = []
    if not tree.GetBranch("n"+obj+"Jets"): 
        return None
    progress = _progress("Reading "+obj+"Jets ...")
    tree.SetBranchStatus("*",0);
    tree.SetBranchStatus("n"+obj+"Jets",1);
    tree.SetBranchStatus(obj+"Jets_pt",1);
    tree.SetBranchStatus(obj+"Jets_eta",1);
    tree.SetBranchStatus(obj+"Jets_phi",1);
    for i in xrange(tree.GetEntries()):
        tree.GetEntry(i)
        number = getattr(tree, "n"+obj+"Jets")
        rawpt,eta,phi = getattr(tree, obj+"Jets_pt"), getattr(tree, obj+"Jets_eta"), getattr(tree, obj+"Jets_phi")
        jets = [ ]
        for j in xrange(number):
            if abs(eta[j]) > etaCut: continue
            if corr:
                pt = corr.correctedPt(rawpt[j], eta[j])
            else:
                pt = rawpt[j]
            if pt > ptCorrCut: 
                jets.append( (pt,eta[j],phi[j]) ) 
        ret.append(calc(jets))
    _cache[_key] = ret
    progress.done("done, %d entries" % len(ret))
    return ret
def makeRecoMETArray(tree, what, obj, etaCut):
    if   etaCut <= 1.5: post = "MetBarrel_pt" 
    elif etaCut <= 2.4: post = "MetCentral_pt"
    else:               post = "Met_pt"
    if not tree.GetBranch(obj+post): 
        return None
    progress = _progress("Reading "+obj+post+" ...")
    tree.SetBranchStatus("*",0);
    tree.SetBranchStatus(obj+post,1);
    ret = []
    for i in xrange(tree.GetEntries()):
        tree.GetEntry(i)
        ret.append(getattr(tree,obj+post))
    progress.done("done, %d entries" % len(ret))
    return ret


def genCutCorrArray(corrArray, genArray, genThr):
    if len(genArray) != len(corrArray): raise RuntimeError("Mismatch")
    return [ r for (g,r) in zip(genArray,corrArray) if g > genThr ]

def makeCumulativeHTEff(name, corrArray, xmax, norm=2760.0*11246/1000):
    if len(corrArray) == 0: return None
    nbins = 2000
    ret = ROOT.TH1F("ceff_"+name.replace(" ","_"), "", nbins, 0., xmax);
    for val in corrArray:
        ret.Fill(min(0.9998*xmax, val))
    tot, msum = float(norm)/len(corrArray), 0
    for ib in xrange(0, nbins):
        msum += ret.GetBinContent(nbins-ib) * tot
        ret.SetBinContent(nbins-ib, msum)
    ret.SetDirectory(None)
    return ret

def makeCumulativeHTEffGenCut(name, corrArray, genArray, genThr, xmax, norm):
    return makeCumulativeHTEff(name, genCutCorrArray(corrArray,genArray,genThr), xmax, norm=norm)

def makeROC(effsig,effbkg):
    graph = ROOT.TGraph(effsig.GetNbinsX())
    for i in xrange(1,graph.GetN()+1):
        graph.SetPoint(i-1, effsig.GetBinContent(i), effbkg.GetBinContent(i))
    return graph 

def makeEffHist(name, refArr, corrArr, corrThr, xmax, logxbins=None):
    if logxbins:
        nbins, nratio = int(logxbins[0]), float(logxbins[1])
        if nratio == 1:
            ret = ROOT.TEfficiency(name.replace(" ","_")+"_eff","",nbins,0,xmax)
        else:
            step = pow(nratio, 1.0/nbins)
            base = xmax*(step-1)/(nratio - 1)
            edges = [0]
            for i in xrange(nbins+1):
                edges.append(edges[-1] + base * pow(step, i))
            ret = ROOT.TEfficiency(name.replace(" ","_")+"_eff","",nbins,array('d',edges))
    else:
        ret = ROOT.TEfficiency(name+"_eff","",20,0,xmax)
    for ref,corr in zip(refArr,corrArr):
        ret.Fill( corr >= corrThr , ref )
    ret.SetStatisticOption(ret.kFCP)
    return ret

from FastPUPPI.NtupleProducer.scripts.respPlots import whats as WHATS
whats = WHATS + [
    ('oldcomp',[
        ("Calo",      "L1OldCalo",        ROOT.kViolet+2, 20, 1.5),
        ("TK 5s",     "L1TKV5",           ROOT.kRed+1, 24, 1.5),
        ("PF",        "L1OldPF",          ROOT.kOrange+7, 24, 1.5),
        ("Puppi",     "L1OldPuppi",       ROOT.kBlue+1, 21, 1.5),
        ("Puppi4MET", "L1OldPuppiForMET", ROOT.kAzure+10, 21, 1.5),
    ]),
    ('newcomp',[
        ("Calo",      "L1Calo",        ROOT.kViolet+2, 20, 1.5),
        ("TK 5s",     "L1TKV5",        ROOT.kRed+1, 24, 1.5),
        ("PF",        "L1PF",          ROOT.kOrange+7, 24, 1.5),
        ("Puppi",     "L1Puppi",       ROOT.kBlue+1, 21, 1.5),
        ("Puppi4MET", "L1PuppiForMET", ROOT.kAzure+10, 21, 1.5),
    ]),
    ('l1pfpu_metref',[
        ("Calo",       "L1Calo$",     ROOT.kViolet+1, 21, 1.5),
        ("TK #Deltaz", "L1TKV5$",     ROOT.kOrange-7, 34, 1.2),
        ("RefTK",      "RefL1TrackerEtMiss$",     ROOT.kGreen+2, 34, 1.2),
        ("PF",         "L1PF$",       ROOT.kOrange+7, 21, 1.5),
        ("Puppi",      "L1Puppi$",    ROOT.kRed+1, 20, 1.1),
    ]),
    ('l1pfpu_jetref',[
        ("Calo",       "L1Calo$",     ROOT.kViolet+1, 21, 1.5),
        ("TK #Deltaz", "L1TKV5$",     ROOT.kOrange-7, 34, 1.2),
        ("RefTK",      "RefTwoLayerJets$",     ROOT.kGreen+2, 34, 1.2),
        ("PF",         "L1PF$",       ROOT.kOrange+7, 21, 1.5),
        ("PFCHS",      "L1CHS$",      ROOT.kMagenta-6, 21, 1.2),
        ("Puppi",      "L1Puppi$",    ROOT.kRed+1, 20, 1.1),
    ]),

]

from optparse import OptionParser
parser = OptionParser("%(prog) infile [ src [ dst ] ]")
parser.add_option("-w", dest="what", default=None, help="Choose set (il1pf, l1pf, ...)")
parser.add_option("-W", dest="what_reg",     default=None, help="Choose set (inputs, l1pf, ...)")
parser.add_option("-P","--plots", dest="plots", default="rate,isorate,roc,effc,plateff,platroc", help="Choose plot or comma-separated list of plots")
parser.add_option("-j","--jecs", dest="jecs", default="jecs.root", help="Choose JEC file")
parser.add_option("--jm","--jec-method", dest="jecMethod", default="", help="Choose JEC method")
parser.add_option("-R","--raw", dest="rawJets", default=False, action="store_true", help="Don't appy JECs")
parser.add_option("-s", dest="genht",  default=None, type="float", help="Choose gen ht")
parser.add_option("-E", dest="eff",  default=None, type="float", help="Choose plateau efficiency")
parser.add_option("-r", dest="rate",  default="10,20,50", type="string", help="Choose rate [kHz] (for isorate plots, can specify more than one)")
parser.add_option("-l","--label", dest="label",  default=None, type="string", help="Extra label for plots")
parser.add_option("-p", "--pt", dest="pt",  default=30, type="float", help="Choose pt cut")
parser.add_option("-e", "--eta", dest="eta",  default=2.4, type="float", help="Choose eta")
parser.add_option("-v", dest="var",  default="ht", help="Choose variable (ht, met, metCentral, mht, jet<N>, mjj, ptj-mjj<M>)")
parser.add_option("--xlabel","--varlabel", dest="varlabel", default=None, help="X axis label for the variable")
parser.add_option("--xmax", dest="xmax",  default=None, type=float, help="Choose variable")
parser.add_option("--logxbins", dest="logxbins",  default=None, nargs=2, type=float, help="--logxbins N X will make N bins, the last being a factor X larger than the first")
options, args = parser.parse_args()

tfiles = [ROOT.TFile.Open(f) for f in args[:2]]

odir = args[2] 
plotter = plotTemplate(odir)

ROOT.gSystem.Load("libL1TriggerPhase2L1ParticleFlow")
ROOT.gInterpreter.ProcessLine('#include "L1Trigger/Phase2L1ParticleFlow/src/corrector.h"')

if options.var == "ht":
    if options.varlabel is None: options.varlabel = "H_{T}"
    if options.genht    is None: options.genht    = 300
    if options.xmax     is None: options.xmax     = 1000
    if options.eff      is None: options.eff      = 0.95
    qualif = "p_{T}^{corr} > %.0f, |#eta| < %.1f" % (options.pt, options.eta)
    what = options.var
elif options.var == "mht":
    if options.varlabel is None: options.varlabel = "H_{T}^{miss}"
    if options.genht    is None: options.genht    = 150
    if options.xmax     is None: options.xmax     = 500
    if options.eff      is None: options.eff      = 0.95
    qualif = "p_{T}^{corr} > %.0f, |#eta| < %.1f" % (options.pt, options.eta)
    what = options.var
elif options.var == "metmht":
    if options.varlabel is None: options.varlabel = "E_{T}^{miss} - H_{T}^{miss}"
    if options.genht    is None: options.genht    = 150
    if options.xmax     is None: options.xmax     = 500
    if options.eff      is None: options.eff      = 0.95
    qualif = "p_{T}^{corr} > %.0f, |#eta| < %.1f" % (options.pt, options.eta)
    what = options.var
elif options.var == "mjj":
    if options.varlabel is None: options.varlabel = "m(jj)"
    if options.genht    is None: options.genht    = 750
    if options.xmax     is None: options.xmax     = 2000
    if options.eff      is None: options.eff      = 0.95
    qualif = "p_{T}^{corr} > %.0f, |#eta| < %.1f" % (options.pt, options.eta)
    what = options.var
elif re.match(r"jet\d+$", options.var):
    ijet = int(options.var.replace("jet",""))
    if options.varlabel is None: options.varlabel = "jet %d p_{T}" % ijet
    if options.genht    is None: options.genht    = 200 if ijet <= 2 else ( 45 if ijet <= 4 else  30)
    if options.xmax     is None: options.xmax     = 400 if ijet <= 2 else (150 if ijet <= 4 else 100)
    if options.eff      is None: options.eff      = 0.95 if ijet <= 2 else 0.90
    options.pt = 10
    qualif = "|#eta| < %.1f" % (options.eta)
    what = options.var
elif options.var.startswith("ptj-mjj"):
    if options.varlabel is None: options.varlabel = "jet p_{T}"
    if options.genht    is None: options.genht    = 100
    if options.xmax     is None: options.xmax     = 300
    if options.eff      is None: options.eff      = 0.95
    options.pt = 10
    qualif = "|#eta| < %.1f, m(jj) > %s" % (options.eta, options.var.replace("ptj-mjj",""))
    what = options.var
elif options.var.startswith("met"):
    if options.varlabel is None: options.varlabel = "E_{T}^{miss}"
    if options.genht    is None: options.genht    = 150
    if options.xmax     is None: options.xmax     = 500
    if options.eff      is None: options.eff      = 0.95
    what = "met"
    options.eta = 5.0 
    if "Central" in options.var:  options.eta = 2.4 
    elif "Barrel" in options.var: options.eta = 1.5
    qualif = "|#eta| < %.1f" % options.eta
else:
    raise RuntimeError("Unknown variable "+options.var)

signal, background = [ f.Get("Events") for f in tfiles ]
jecfile = ROOT.TFile.Open(options.jecs)

for plotkind in options.plots.split(","):
  print "Make plot "+plotkind
  if plotkind != "rate":
    genArray = makeGenArray(signal, what, options.pt, options.eta)
  rates = map(float, options.rate.split(","))
  if plotkind != "isorate": rates = rates[:1]
  for targetrate in rates:
      for objset,things in whats:
          if options.what and (objset not in options.what.split(",")): continue
          if options.what_reg:
              if not any(re.match(p+"$",objset) for p in options.what_reg.split(",")): 
                  continue
          plots = []
          for name,obj,col,msty,msiz in things:
              if "GenAcc$" in obj: continue
              obj = obj.replace("$","")
              if options.var.startswith("met") or "RefTwoLayerJets" in obj:
                  jecs = None
              else:
                  jecdirname = obj+"Jets"+( "_"+options.jecMethod if options.jecMethod else "")
                  jecdir = jecfile.GetDirectory(jecdirname)
                  if not jecdir: 
                      print "Missing JECs "+jecdirname+" in "+options.jecs
                      continue
                  jecs = ROOT.l1tpf.corrector(jecdir)
              label = name
              ptcut = options.pt
              if "RefTwoLayerJets" in obj: ptcut = 5
              if plotkind == "rate":
                  recoArrayB = makeCorrArray(background, what, obj, ptcut, options.eta, jecs)
                  if not recoArrayB: continue
                  plot = makeCumulativeHTEff(name, recoArrayB, options.xmax)
              elif plotkind == "effc":
                  recoArrayS = makeCorrArray(signal, what, obj, ptcut, options.eta, jecs)
                  if not recoArrayS: continue
                  plot = makeCumulativeHTEffGenCut(name, recoArrayS, genArray, options.genht, options.xmax, norm=1)
              elif plotkind == "isorate":
                  recoArrayB = makeCorrArray(background, what, obj, ptcut, options.eta, jecs)
                  if not recoArrayB: continue
                  rateplot = makeCumulativeHTEff(name, recoArrayB, options.xmax)
                  cut = 9999
                  for ix in xrange(1,rateplot.GetNbinsX()+1):
                      if rateplot.GetBinContent(ix) <= targetrate:
                          cut = rateplot.GetXaxis().GetBinLowEdge(ix)
                          break
                  recoArrayS = makeCorrArray(signal, what, obj, ptcut, options.eta, jecs)
                  if not recoArrayS: continue
                  plot = makeEffHist(name, genArray, recoArrayS, cut, options.xmax, logxbins=options.logxbins)
                  label = "%s(%s) > %.0f" % (options.varlabel, name,cut)
              elif plotkind == "roc":
                  recoArrayB = makeCorrArray(background, what, obj, ptcut, options.eta, jecs)
                  if not recoArrayB: continue
                  recoArrayS = makeCorrArray(signal,     what, obj, ptcut, options.eta, jecs)
                  if not recoArrayS: continue
                  effsig  = makeCumulativeHTEffGenCut(name+"_s", recoArrayS, genArray, options.genht, options.xmax, norm=1)
                  ratebkg = makeCumulativeHTEff(name+"_b", recoArrayB, options.xmax)
                  if not effsig or not ratebkg: continue
                  plot = makeROC(effsig,ratebkg)
                  msty = 0
              elif plotkind == "plateff":
                  recoArrayB = makeCorrArray(background, what, obj, ptcut, options.eta, jecs)
                  if not recoArrayB: continue
                  rateplot = makeCumulativeHTEff(name, recoArrayB, options.xmax)
                  def effForRate(rate):
                    cut = 9999
                    for ix in xrange(1,rateplot.GetNbinsX()+1):
                        if rateplot.GetBinContent(ix) <= rate:
                            cut = rateplot.GetXaxis().GetBinLowEdge(ix)
                            break
                    #print "Cut for %s @ rate %g: %g" % (what+obj, rate, cut)
                    recoArrayS = makeCorrArray(signal, what, obj, ptcut, options.eta, jecs)
                    if not recoArrayS: return None
                    plot = makeEffHist(name, genArray, recoArrayS, cut, options.xmax, logxbins=options.logxbins)
                    eff = plot.GetEfficiency(plot.FindFixBin(options.genht))
                    #print "Eff at %g: %g" % (options.genht, eff)
                    return eff, plot, cut
                  rate = 50
                  while True:
                    eff, plot, cut = effForRate(rate)
                    if not eff: break
                    if eff < options.eff:
                        if rate >= 30e3: 
                            eff = None; break
                        rate = min(rate*5, 30e3)
                    else:
                        break
                  if not eff: continue
                  #print "Upper bound: eff %g at rate %g" % (eff, rate) 
                  maxrate = rate; rate = rate / 5;
                  while True:
                    eff, plot, cut = effForRate(rate)
                    if not eff: break
                    if eff > options.eff:
                        rate /= 2; 
                    else:
                        break
                  if not eff: continue                 
                  #print "Lower bound: eff %g at rate %g" % (eff, rate) 
                  minrate = rate
                  while True:
                    rate = sqrt(maxrate * minrate)
                    eff, plot, cut = effForRate(rate)
                    if not eff: 
                        #print "Bisection failed?"
                        break
                    if eff < options.eff:
                        #print "New lower bound: eff %g at rate %g" % (eff, rate)
                        minrate = rate;
                    else:
                        #print "New upper bound: eff %g at rate %g" % (eff, rate)
                        maxrate = rate
                    if maxrate/minrate < 1.1: 
                        break
                  if not eff: continue
                  if rate < 10:    ratestr = "%.1fkHz" % rate
                  elif rate < 500: ratestr = "%.0fkHz" % rate
                  else:            ratestr = "%.1fMHz" % (rate/1000)
                  label = "%s @ %s" % (name, ratestr)
              elif plotkind == "platroc":
                  recoArrayB = makeCorrArray(background, what, obj, ptcut, options.eta, jecs)
                  if not recoArrayB: continue
                  rateplot = makeCumulativeHTEff(name, recoArrayB, options.xmax)
                  def platForRate(rate):
                    cut = None
                    for ix in xrange(1,rateplot.GetNbinsX()+1):
                        if rateplot.GetBinContent(ix) <= rate:
                            cut = rateplot.GetXaxis().GetBinLowEdge(ix)
                            break
                    if cut is None:
                        return None
                    #print "Cut for %s @ rate %g: %g" % (what+obj, rate, cut)
                    recoArrayS = makeCorrArray(signal, what, obj, ptcut, options.eta, jecs)
                    if not recoArrayS: return None
                    plot = makeEffHist(name, genArray, recoArrayS, cut, options.xmax, logxbins=options.logxbins)
                    hist = plot.GetTotalHistogram(); xaxis = hist.GetXaxis()
                    for i in xrange(2,hist.GetNbinsX()+1):
                        if plot.GetEfficiency(i) > options.eff:
                            graph = plot.CreateGraph()
                            xmin, xmax = xaxis.GetBinLowEdge(i-1), xaxis.GetBinCenter(i)
                            if graph.Eval(xmin,None,"S") > options.eff: 
                                #print "ERROR xmin for %s @ %g (%g): i = %d" % (what+obj, rate, cut, i)
                                pass
                            elif graph.Eval(xmax,None,"S") < options.eff:
                                #print "ERROR xmax for %s @ %g (%g): i = %d" % (what+obj, rate, cut, i)
                                pass
                            else:
                                xmid = 0.5*(xmax+xmin)
                                while abs(xmax-xmin) > 0.05*(abs(xmax)+abs(xmin)):
                                    if graph.Eval(xmid,None,"S") < options.eff:
                                        xmin = xmid
                                    else:
                                        xmax = xmid
                                    xmid = 0.5*(xmax+xmin)
                                #print "Plateau at %g eff for %s @ %g: %g" % (options.eff, what+obj, rate, xmid)
                                return xmid
                            break
                    return None
                  points = []
                  rate = 1
                  while rate <= 30e3:
                    plat = platForRate(rate)
                    if plat: 
                        if len(points) == 0 or points[-1][0] != plat:
                            points.append((plat,rate))
                    rate *= 1.2
                  if not points: continue
                  plot = ROOT.TGraph(len(points))
                  for i,(x,y) in enumerate(points):
                    plot.SetPoint(i,x,y)
                  label = name
              else: raise RuntimeError
              if not plot: continue
              plot.SetLineWidth(3); plot.SetLineColor(col);  plot.SetMarkerColor(col)
              plot.SetMarkerStyle(msty); plot.SetMarkerSize(msiz)
              plots.append((label,plot))
          if not plots: 
              print "   nothing to plot!"
              continue
          plotter.SetLogy(False)
          if plotkind == "rate":
              plotter.SetLogy(True)
              frame = ROOT.TH1D("",";L1 %s cut (%s); Minbias rate @ PU200 [kHz]" % (options.varlabel, qualif), 100, 0, options.xmax)
              frame.GetYaxis().SetDecimals(True)
              frame.GetXaxis().SetNdivisions(505)
              frame.GetYaxis().SetRangeUser(0.5, 100e3)
              leg = ROOT.TLegend(0.56,0.93,0.93,0.93-0.055*len(things))
              plotname = '%s%s-%s_eta%s_pt%d' % (options.var, plotkind, objset, options.eta, options.pt)
          elif plotkind == "effc":
              gentext, genpost = "iciency", "" 
              if options.genht > 0:
                  gentext = " (Gen %s > %s)" % (options.varlabel, options.genht)
                  genpost = "_gen%.0f" % (options.genht)
              frame = ROOT.TH1D("",";L1 %s thresh (%s); Eff%s" % (options.varlabel, qualif, gentext), 100, 0, options.xmax)
              leg = ROOT.TLegend(0.6,0.93,0.93,0.93-0.055*len(things))
              plotname = '%s%s-%s_eta%s_pt%d%s' % (options.var, plotkind, objset, options.eta, options.pt, genpost)
          elif plotkind == "isorate":
              frame = ROOT.TH1D("",";Gen %s (%s); Eff (L1 rate %.0f kHz)" % (options.varlabel, qualif, targetrate), 100, 0, options.xmax)
              frame.GetYaxis().SetDecimals(True)
              leg = ROOT.TLegend(0.50,0.18,0.93,0.18+0.055*len(things))
              plotname = '%s%s-%s_eta%s_pt%d_%.0fkHz' % (options.var, plotkind, objset, options.eta, options.pt, targetrate)
          elif plotkind == "plateff":
              frame = ROOT.TH1D("",";Gen %s (%s); Eff" % (options.varlabel, qualif), 100, 0, options.xmax)
              frame.GetYaxis().SetDecimals(True)
              leg = ROOT.TLegend(0.30,0.18,0.93,0.18+0.055*len(things))
              plotname = '%s%s-%s_eta%s_pt%d_eff%.2fat%g' % (options.var, plotkind, objset, options.eta, options.pt, options.eff, options.genht)
          elif plotkind == "roc":
              plotter.SetLogy(True)
              gentext, genpost = "iciency", "" 
              if options.genht > 0:
                  gentext = " (Gen %s > %s)" % (options.varlabel, options.genht)
                  genpost = "_gen%.0f" % (options.genht)
              frame = ROOT.TH1D("",";Eff%s; Minbias rate @ PU200 [kHz]" % (gentext), 100, 0, 1)
              frame.GetYaxis().SetDecimals(True)
              frame.GetXaxis().SetNdivisions(505)
              frame.GetYaxis().SetRangeUser(0.5, 100e3)
              leg = ROOT.TLegend(0.2,0.93,0.55,0.93-0.055*len(things))
              plotname = '%s%s-%s_eta%s_pt%d%s' % (options.var, plotkind, objset, options.eta, options.pt, genpost)
          elif plotkind == "platroc":
              plotter.SetLogy(True)
              xtit = "Gen threshold at %g%% eff" % (options.eff*100)
              frame = ROOT.TH1D("",";%s; Minbias rate @ PU200 [kHz]" % (xtit), 100, 0, options.xmax)
              frame.GetYaxis().SetDecimals(True)
              #frame.GetXaxis().SetNdivisions(505)
              frame.GetYaxis().SetRangeUser(0.5, 100e3)
              leg = ROOT.TLegend(0.6,0.93,0.93,0.93-0.055*len(things))
              plotname = '%s%s-%s_eta%s_pt%d_eff%.3f' % (options.var, plotkind, objset, options.eta, options.pt, options.eff)
          frame.Draw()
          if plotkind in ("rate","roc","platroc"):
              line = ROOT.TLine()
              line.SetLineStyle(7)
              for y in 40e3, 100, 10:
                  line.DrawLine(frame.GetXaxis().GetXmin(),y,frame.GetXaxis().GetXmax(),y)
          for n,p in plots: 
              if plotkind == "platroc":
                  p.Draw("PLX SAME")
              else:
                  p.Draw("PCX SAME" if ("TH1" not in p.ClassName()) else "C SAME")
          for n,p in plots: 
              leg.AddEntry(p, n, "L" if plotkind in ("rate","roc") else "LP")
          leg.Draw()
          plotter.decorations()
          plotter.Print('%s%s' % (plotname, ("_"+options.label) if options.label else ""))
          fout = ROOT.TFile.Open('%s/%s%s.root' % (odir, plotname, ("_"+options.label) if options.label else ""), "RECREATE")
          fout.WriteTObject(frame,"frame")
          for n,p in plots: 
              p.SetTitle(n)
              fout.WriteTObject(p)
          fout.Close()
          del frame



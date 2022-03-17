import re, os
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kWarning

from array import array
from math import pow, sqrt
from FastPUPPI.NtupleProducer.scripts.makeJecs import _progress
from FastPUPPI.NtupleProducer.plotTemplate import plotTemplate

ROOT.gROOT.ProcessLine('#include "%s/src/FastPUPPI/NtupleProducer/python/scripts/jetHtSuite.h"' % os.environ['CMSSW_BASE']);

def makeCalcCpp(what):
    if what == "ht": 
        return ROOT.CalcHT()
    if what == "mht": 
        return ROOT.CalcMHT()
    if what == "mjj": 
        return ROOT.CalcMJJ()
    if re.match(r"jet\d+$", what): 
        return ROOT.CalcJ(int(what.replace("jet","")))
    if what.startswith("ptj-mjj"): 
        return ROOT.CalcJ2_MJJcut(float(what.replace("ptj-mjj","")))
    raise RuntimeError("CalcCpp: can't parse "+what)

def makeLepCalcCpp(what):
    if re.match(r"lep\d*_(pt|iso|id|eta|abseta)", what): 
        (a,b) = what.split("_")
        if a == "lep": a = "lep0"
        return ROOT.CalcL(int(a[3:]), {"pt":0,"iso":1,"id":2,"eta":3,"abseta":4}[b])
    raise RuntimeError("LepCalcCpp: can't parse "+what)

def makeGenArray(tree, what, ptCut, etaCut):
    return makeCorrArray(tree, what, "Gen", ptCut, etaCut, ROOT.nullptr)

def makeCorrArray(tree, what, obj, ptCorrCut, etaCut, corr, _cache={}):
    _key = (id(tree),what,obj,int(ptCorrCut*100),int(etaCut*1000))
    if _key in _cache: return _cache[_key]
    if what == "metmht":
        met = makeCorrArray(tree, "met", obj,     0,      5.0,   corr, _cache=_cache)
        mht = makeCorrArray(tree, "mht", obj, ptCorrCut, etaCut, corr, _cache=_cache)
        ret = ROOT.makeMinimum(met,mht)
        _cache[_key] = ret
        return ret
    if "met" in what:
        ret = makeMETArray(tree, what, obj, etaCut)
        _cache[_key] = ret
        return ret
    if not tree.GetBranch("n"+obj+"Jets"):
        if obj == "Gen": raise RuntimeError("Missing GenJets!");
        return None
    cppcalc = makeCalcCpp(what)
    progress = _progress("  Reading "+obj+"Jets in C++...")
    ret = ROOT.makeJetArray(tree, obj, ptCorrCut, etaCut, cppcalc, corr);
    _cache[_key] = ret
    progress.done("done, %d entries" % len(ret))
    return ret

def makeMETArray(tree, what, obj, etaCut):
    if obj == "Gen": obj = "gen" # fix issue with naming convention
    if   etaCut <= 1.5: post = "MetBarrel" 
    elif etaCut <= 2.4: post = "MetCentral"
    else:               post = "Met"
    if not tree.GetBranch(obj+post+"_pt"): 
        if obj == "gen": raise RuntimeError("Missing gen"+post);
        return None
    progress = _progress("  Reading "+obj+post+" in C++...")
    ret = ROOT.makeMetArray(tree, obj+post);
    progress.done("done, %d entries" % len(ret))
    return ret

def makeGenLepArray(tree, what, genobj, ptCut, etaCut, prompt=True):
    return makeRecoLepArray(tree, what, genobj, ptCut, etaCut, lepid=("prompt", 2 if prompt else 0, False))

def parseLepRequirement(obj,ptCut,lepid,lepiso):
    if ":" not in obj: return (obj,ptCut,lepid,lepiso)
    obj, reqs = obj.split(":")
    for req in reqs.split(";"):
        m = re.match(r"id=(\w+),([0-9.+-]+)([if]?)", req)
        if m:
            lepid = (m.group(1), float(m.group(2)), m.group(3) == "f")
            continue
        m = re.match(r"pt>=([0-9.]+)", req)
        if m:
            ptCut = float(m.group(1))
            continue
        raise RuntimeError("Can't parse requirement %r in obj %r" % (obj,req))
    return (obj,ptCut,lepid,lepiso)

def makeRecoLepArray(tree, what, obj, ptCut, etaCut, lepid=("",-1.,True), lepiso=("", -1., False), _cache={}):
    (obj,ptCut,lepid,lepiso) = parseLepRequirement(obj,ptCut,lepid,lepiso)
    _key = (id(tree),what,obj,int(ptCut*100),int(etaCut[0]*1000),int(etaCut[1]*1000),lepid[0],int(lepid[1]*1000),lepiso[0],int(lepiso[1]*1000))
    if _key in _cache: return _cache[_key]
    if not tree.GetBranch("n"+obj):
        if "Gen" in obj: raise RuntimeError("Missing "+obj+"!");
        return None
    cppcalc = makeLepCalcCpp(what)
    progress = _progress("  Reading "+obj+" in C++...")
    ret = ROOT.makeLepArray(tree, obj, ptCut, etaCut[0], etaCut[1], cppcalc, lepid[0], lepid[1], lepid[2], lepiso[0], lepiso[1], lepiso[2]);
    _cache[_key] = ret
    progress.done("done, %d entries" % len(ret))
    return ret

def makeMatchedRecoLepArray(tree, what, 
                            obj, ptCut, etaCut, lepid, lepiso,
                            genObj, genPtCut, genEtaCut, genPrompt, genDr, _cache={}):
    (obj,ptCut,lepid,lepiso) = parseLepRequirement(obj,ptCut,lepid,lepiso)
    _key = (id(tree),what,obj,int(ptCut*100),int(etaCut[0]*1000),int(etaCut[1]*1000),lepid[0],int(lepid[1]*1000),lepiso[0],int(lepiso[1]*1000))
    if _key in _cache: return _cache[_key]
    if not tree.GetBranch("n"+obj):
        if "Gen" in obj: raise RuntimeError("Missing "+obj+"!");
        return None
    cppcalc = makeLepCalcCpp(what)
    progress = _progress("  Reading "+obj+" in C++...")
    ret = ROOT.makeMatchedLepArray(tree, cppcalc, 
                            obj, ptCut, etaCut[0], etaCut[1], lepid[0], lepid[1], lepid[2], lepiso[0], lepiso[1], lepiso[2],
                            genObj, genPtCut, genEtaCut[0], genEtaCut[1], 2 if genPrompt else 0, genDr);
    _cache[_key] = ret
    progress.done("done, %d entries" % len(ret))
    return ret

def makeCumulativeHTEff(name, corrArray, xmax, norm=2760.0*11246/1000):
    return makeCumulativeHTEffGenCut(name, corrArray, None, None, xmax, norm)

def makeCumulativeHTEffGenCut(name, corrArray, genArray, genThr, xmax, norm):
    if len(corrArray) == 0: return None
    nbins = 2000
    ret = ROOT.TH1F("ceff_"+name.replace(" ","_"), "", nbins, 0., xmax);
    if genArray:
        nsel = ROOT.fillTH1FastGenCut(ret, corrArray, genArray, genThr)
    else:
        nsel = ROOT.fillTH1Fast(ret, corrArray)
    if nsel == 0: return None
    tot, msum = float(norm)/nsel, 0
    for ib in range(0, nbins):
        msum += ret.GetBinContent(nbins-ib) * tot
        ret.SetBinContent(nbins-ib, msum)
    ret.SetDirectory(ROOT.nullptr)
    return ret

def makeInclusiveEffGenHTCut(name, corrArray, corrThr, genArray, genThr):
    if len(corrArray) == 0: return None
    pair = ROOT.inclusiveEffFastGenCut(genArray, genThr, corrArray, corrThr)
    num, den = pair.first, pair.second
    eff = num/float(den); err = sqrt(eff*(1-eff)/den)
    print("%s, L1 cut > %g (gen cut > %g) --> %d/%d = %.4f +/- %.4f\n" % (name, corrThr, genThr, num, den, eff, err))

def makeInclusiveEffRate(name, corrArray, corrThr, norm=2760.0*11246/1000):
    if len(corrArray) == 0: return None
    pair = ROOT.inclusiveEffFast(corrArray, corrThr)
    num, den = pair.first, pair.second
    eff = num/float(den); err = sqrt(eff*(1-eff)/den)
    print("%s, L1 cut > %g --> %d/%d = %.4f +/- %.4f, rate %.2f +/- %.2f kHz\n" % (name, corrThr, num, den, eff, err, norm*eff, norm*err))

def makeEffHist(name, refArr, corrArr, corrThr, xmax, logxbins=None):
    if logxbins:
        nbins, nratio = int(logxbins[0]), float(logxbins[1])
        if nratio == 1:
            ret = ROOT.TEfficiency(name.replace(" ","_")+"_eff","",nbins,0,xmax)
        else:
            step = pow(nratio, 1.0/nbins)
            base = xmax*(step-1)/(nratio - 1)
            edges = [0]
            for i in range(nbins+1):
                edges.append(edges[-1] + base * pow(step, i))
            ret = ROOT.TEfficiency(name.replace(" ","_")+"_eff","",nbins,array('d',edges))
    else:
        ret = ROOT.TEfficiency(name+"_eff","",20,0,xmax)
    ROOT.fillTEffFast(ret, refArr, corrArr, corrThr)
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
        ("TK #Deltaz", "L1TKV5$",     ROOT.kGreen+1, 34, 1.2),
        ("RefTK",      "RefL1TrackerEtMiss$",     ROOT.kGreen+3, 34, 1.2),
        ("Puppi",      "L1Puppi$",    ROOT.kRed+1, 20, 1.1),
    ]),
    ('l1pfpu_metnoref',[
        ("Calo",       "L1Calo$",     ROOT.kViolet+1, 21, 1.5),
        ("TK #Deltaz", "L1TKV5$",     ROOT.kGreen+2, 34, 1.2),
        ("Puppi",      "L1Puppi$",    ROOT.kRed+1, 20, 1.1),
    ]),
    ('l1pfpu_metrefonly',[
        ("TK",      "RefL1TrackerEtMiss$",  ROOT.kGreen+2, 34, 1.2),
        ("Puppi",   "L1Puppi$",             ROOT.kRed+1, 20, 1.1),
    ]),
    ('l1pfpu_jetnoref',[
        ("Calo",       "L1Calo$",    ROOT.kViolet+1, 21, 1.5),
        ("TK #Deltaz", "L1TKV5$",    ROOT.kGreen+2, 34, 1.2),
        ("ak4Puppi",   "L1Puppi$",   ROOT.kRed+1, 20, 1.1),
    ]),
    ('l1pfpu_jetref',[
        ("Calo",       "L1Calo$",              ROOT.kViolet+1, 21, 1.5),
        ("RefCalo",    "RefCaloJets$",         ROOT.kViolet+2, 21, 1.5),
        ("TK #Deltaz", "L1TKV5$",              ROOT.kGreen+1, 34, 1.2),
        ("RefTK",      "RefTwoLayerJets$",     ROOT.kGreen+3, 34, 1.2),
        ("ak4Puppi",   "L1Puppi$",             ROOT.kRed+0, 20, 1.1),
        ("RefPuppi",   "RefPhase1PuppiJets$",  ROOT.kRed+2, 20, 0.9), 
    ]),
    ('l1pfpu_jetrefonly',[
        ("Calo",    "RefCaloJets$",         ROOT.kViolet+1, 21, 1.5),
        ("TK",      "RefTwoLayerJets$",     ROOT.kGreen+2, 34, 1.2),
        ("Puppi",   "RefPhase1PuppiJets$",  ROOT.kRed+1, 20, 0.9), 
    ]),
    ('l1pfpu_bitwise',[
        ("CMSSW",       "L1CMSSWPuppi$",      ROOT.kGray+2, 20, 2.0),
        ("Tuned",       "L1Puppi$",           ROOT.kRed+1, 20, 1.7),
        ("Regional",    "L1PuppiRegional$",   ROOT.kViolet+1, 20, 1.4),
        ("BitwisePF",   "L1PuppiBitwisePF$",  ROOT.kAzure+1, 20, 1.1),
        ("BitwisePuppi","L1PuppiBitwise$",    ROOT.kGreen+2, 20, 0.8),
    ]),
    ('l1pfpu_tkemu',[
        ("Default",    "L1Puppi$",       ROOT.kGray+2, 20, 2.0),
        ("Emu Flt",    "L1PuppiTkFlt$",  ROOT.kRed+1, 20, 1.7),
        ("Emu Int",    "L1PuppiTkEmu$",  ROOT.kAzure+1, 20, 1.4),
    ]),
    ('l1pfmu_test',[
        ("Gen",  "GenMu",   ROOT.kGray+2, 20, 2.0),
        ("Sta",  "StaMu",   ROOT.kRed+1,  20, 1.8),
        ("PF",   "PFMu",    ROOT.kGreen+1, 20, 1.5),
        #("PF q12", "PFMu:id=quality,12",    ROOT.kBlue+1, 20, 1.2), # example of quality requirement
        #("PFpt30", "PFMu:pt>=30",    ROOT.kRed+1, 20, 1.0),
    ]),
    ('l1tkeg_test_both',[
        ("Gen",      "GenEl",   ROOT.kGray+2, 20, 2.0),
        ("TkEM EB",  "TkEmEB",   ROOT.kRed+0,  20, 1.7),
        ("TkEM EE",  "TkEmEE",   ROOT.kRed+2,  20, 1.7),
        ("TkEle EB",  "TkEleEB",   ROOT.kGreen+0,  20, 1.3),
        ("TkEle EE",  "TkEleEE",   ROOT.kGreen+2,  20, 1.3),
    ]),
    ('l1tkeg_test_EB',[
        ("Gen",      "GenEl",   ROOT.kGray+2, 20, 2.0),
        ("TkEM EB",  "TkEmEB",   ROOT.kRed+0,  20, 1.7),
        ("TkEle EB",  "TkEleEB",   ROOT.kGreen+0,  20, 1.3),
    ]),
    ('l1tkeg_test_EE',[
        ("Gen",      "GenEl",   ROOT.kGray+2, 20, 2.0),
        ("TkEM EE",  "TkEmEE",   ROOT.kRed+2,  20, 1.7),
        ("TkEle EE",  "TkEleEE",   ROOT.kGreen+2,  20, 1.3),
    ]),
    ('l1tkeg_test_both_tkemu',[
        ("Gen",      "GenEl",   ROOT.kGreen+0, 20, 2.0),
        ("EleEB Sim",  "TkEleEB",        ROOT.kGray+1,    20, 1.8),
        ("EleEB Flt",  "TkEleEBTkFlt",   ROOT.kRed+0,     20, 1.4),
        ("EleEB Emu",  "TkEleEBTkEmu",   ROOT.kAzure+10,  20, 1.0),
        ("EleEE Sim",  "TkEleEE",        ROOT.kGray+3,    20, 1.8),
        ("EleEE Flt",  "TkEleEETkFlt",   ROOT.kRed+2,     20, 1.4),
        ("EleEE Emu",  "TkEleEETkEmu",   ROOT.kAzure+2,   20, 1.0),
    ]),
    ('l1tkeg_test_EB_tkemu',[
        ("Gen",      "GenEl",   ROOT.kGreen+0, 20, 2.0),
        ("EleEB Sim",  "TkEleEB",        ROOT.kGray+1,    20, 1.8),
        ("EleEB Flt",  "TkEleEBTkFlt",   ROOT.kRed+0,     20, 1.4),
        ("EleEB Emu",  "TkEleEBTkEmu",   ROOT.kAzure+10,  20, 1.0),
    ]),
    ('l1tkeg_test_EE_tkemu',[
        ("Gen",      "GenEl",   ROOT.kGreen+0, 20, 2.0),
        ("EleEE Sim",  "TkEleEE",        ROOT.kGray+3,    20, 1.8),
        ("EleEE Flt",  "TkEleEETkFlt",   ROOT.kRed+2,     20, 1.4),
        ("EleEE Emu",  "TkEleEETkEmu",   ROOT.kAzure+2,   20, 1.0),
    ]),
    ('l1pfpu_scjetemu',[
        ("Ak4",       "L1Puppi$",      ROOT.kGray+2, 20, 2.0),
        ("scPuppi",   "scPuppi$",      ROOT.kRed+1, 20, 1.7),
        ("scEmu",     "scEmuPuppi$",   ROOT.kViolet+1, 20, 1.4),
        ("scDereg",   "scDeregPuppi$", ROOT.kAzure+1, 20, 1.1),
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
parser.add_option("-E", dest="eff",  default=None, type="string", help="Choose plateau efficiency")
parser.add_option("-r", dest="rate",  default="10,20,50", type="string", help="Choose rate [kHz] (for isorate plots, can specify more than one)")
parser.add_option("-l","--label", dest="label",  default=None, type="string", help="Extra label for plots")
parser.add_option("-p", "--pt", dest="pt",  default=30, type="float", help="Choose pt cut")
parser.add_option("-e", "--eta", "--etamax", dest="eta",  default=2.4, type="float", help="Choose max |eta|")
parser.add_option("--etamin", dest="etamin",  default=0, type="float", help="Choose minimum |eta|")
parser.add_option("--genprompt", dest="genprompt",  default=2, type="int", help="Gen prompt: 2 = direct, 1 = also tau, 0 = any")
parser.add_option("--gendr", dest="gendr",  default=0.1, type="float", help="Delta(R) match between gen and reco")
parser.add_option("--lepid", dest="lepid",  default=None, help="Choose lep id")
parser.add_option("--lepiso", dest="lepiso",  default=None, help="Choose lep iso")
parser.add_option("-v", dest="var",  default="ht", help="Choose variable (ht, met, metCentral, mht, jet<N>, mjj, ptj-mjj<M>)")
parser.add_option("--xv", dest="xvar",  default=None, help="Choose a gen variable for the x axis in lepton efficiency plots)")
parser.add_option("--xlabel","--varlabel", dest="varlabel", default=None, help="X axis label for the variable")
parser.add_option("--xmax", dest="xmax",  default=None, type=float, help="X axis maximum")
parser.add_option("--ymin", dest="ymin",  default=None, type=float, help="Y axis minimum vale")
parser.add_option("--logxbins", dest="logxbins",  default=None, nargs=2, type=float, help="--logxbins N X will make N bins, the last being a factor X larger than the first")
parser.add_option("--print", dest="printQualityPlots", action="store_true", default=False, help="Make print-quality plots (incl. pdf & eps)")
options, args = parser.parse_args()

tfiles = [ROOT.TFile.Open(f) for f in args[:2]]

odir = args[2] 
plotter = plotTemplate(odir, defaultExts = (["png","eps","pdf"] if options.printQualityPlots else ["png"]))

ROOT.gSystem.Load("libL1TriggerPhase2L1ParticleFlow")
ROOT.gInterpreter.ProcessLine('#include "L1Trigger/Phase2L1ParticleFlow/src/corrector.h"')

isJetMet = True
if options.var == "ht":
    if options.varlabel is None: options.varlabel = "H_{T}"
    if options.genht    is None: options.genht    = 300
    if options.xmax     is None: options.xmax     = 1000
    if options.eff      is None: options.eff      = "0.5,0.9,0.95"
    qualif = "p_{T}^{corr} > %.0f, |#eta| < %.1f" % (options.pt, options.eta)
    what = options.var
elif options.var == "mht":
    if options.varlabel is None: options.varlabel = "H_{T}^{miss}"
    if options.genht    is None: options.genht    = 150
    if options.xmax     is None: options.xmax     = 500
    if options.eff      is None: options.eff      = "0.5,0.9,0.95"
    qualif = "p_{T}^{corr} > %.0f, |#eta| < %.1f" % (options.pt, options.eta)
    what = options.var
elif options.var == "metmht":
    if options.varlabel is None: options.varlabel = "E_{T}^{miss} - H_{T}^{miss}"
    if options.genht    is None: options.genht    = 150
    if options.xmax     is None: options.xmax     = 500
    if options.eff      is None: options.eff      = "0.5,0.9,0.95"
    qualif = "p_{T}^{corr} > %.0f, |#eta| < %.1f" % (options.pt, options.eta)
    what = options.var
elif options.var == "mjj":
    if options.varlabel is None: options.varlabel = "m(jj)"
    if options.genht    is None: options.genht    = 750
    if options.xmax     is None: options.xmax     = 2000
    if options.eff      is None: options.eff      = "0.9,0.95"
    qualif = "p_{T}^{corr} > %.0f, |#eta| < %.1f" % (options.pt, options.eta)
    what = options.var
elif re.match(r"jet\d+$", options.var):
    ijet = int(options.var.replace("jet",""))
    if options.varlabel is None: options.varlabel = "jet %d p_{T}" % ijet
    if options.genht    is None: options.genht    = 200 if ijet <= 2 else ( 45 if ijet <= 4 else  30)
    if options.xmax     is None: options.xmax     = 400 if ijet <= 2 else (150 if ijet <= 4 else 100)
    if options.eff      is None: options.eff      = "0.95" if ijet <= 2 else "0.9,0.95"
    options.pt = 10
    qualif = "|#eta| < %.1f" % (options.eta)
    what = options.var
elif options.var.startswith("ptj-mjj"):
    if options.varlabel is None: options.varlabel = "jet p_{T}"
    if options.genht    is None: options.genht    = 100
    if options.xmax     is None: options.xmax     = 300
    if options.eff      is None: options.eff      = "0.9,0.95"
    options.pt = 20 # slightly increase the cut, to avoid too large combinatoric.
    qualif = "|#eta| < %.1f, m(jj) > %s" % (options.eta, options.var.replace("ptj-mjj",""))
    what = options.var
elif options.var.startswith("met"):
    if options.varlabel is None: options.varlabel = "E_{T}^{miss}"
    if options.genht    is None: options.genht    = 150
    if options.xmax     is None: options.xmax     = 500
    if options.eff      is None: options.eff      = "0.5,0.9,0.95"
    what = "met"
    options.eta = 5.0 
    if "Central" in options.var:  options.eta = 2.4 
    elif "Barrel" in options.var: options.eta = 1.5
    qualif = "|#eta| < %.1f" % options.eta
elif options.var.startswith("lep"):
    isJetMet = False
    what = options.var
    if options.xvar is None: 
        options.xvar = options.var
    if options.xvar.endswith("pt"):
        if options.varlabel is None: options.varlabel = "p_{T}"
        if options.genht    is None: options.genht    = 18
        if options.xmax     is None: options.xmax     = 80
        if options.eff      is None: options.eff      = "0.95"
        options.genleppt = 0
        qualif = "|#eta| < %.1f" % options.eta
    elif options.xvar.endswith("abseta"):
        options.eta = 3.2
        if options.varlabel is None: options.varlabel = "|#eta|"
        if options.genht    is None: options.genht    = 25
        if options.xmax     is None: options.xmax     = options.eta
        if options.eff      is None: options.eff      = "0.95"
        options.genleppt = options.genht
        qualif = "p_{T}^{L1} > %.1f, p_{T}^{gen} > %.1f" % (options.pt, options.genht)
    else:
        raise RuntimeError("Unknown variable "+options.var)
    if options.var.endswith("pt"):
        options.pt = 0
    if options.lepid is None:
        options.lepid_ = ("", -1., True)
    if options.lepiso is None:
        options.lepiso_ = ("", -1., False)
    # finally, convert this into a tuple for convenience
    options.etas = (options.etamin, options.eta)

else:
    raise RuntimeError("Unknown variable "+options.var)

signal, background = [ f.Get("Events") for f in tfiles ]
jecfile = ROOT.TFile.Open(options.jecs) if isJetMet else None


def makePlatEffPlot(signal, background, what, obj, ptcut, jecs, plotparam, _cache={}):
    _key = (id(signal),id(background),what,obj,str(ptcut),str(plotparam))
    if _key in _cache: 
        #print "  retrieved plateff for %s, %s, %s from cache" % (what, obj, plotparam)
        return _cache[_key]
    # ok there we go
    recoArrayB = makeCorrArray(background, what, obj, ptcut, options.eta, jecs)
    if not recoArrayB: return (None, None)
    rateplot = makeCumulativeHTEff(name, recoArrayB, options.xmax)
    recoArrayS = makeCorrArray(signal, what, obj, ptcut, options.eta, jecs)
    if not recoArrayS: return (None, None)
    #platprogress = _progress("  making plateff for %s, %s, %s" % (what, obj, plotparam))
    def effForRate(rate):
      cut = 9999
      for ix in range(1,rateplot.GetNbinsX()+1):
          if rateplot.GetBinContent(ix) <= rate:
              cut = rateplot.GetXaxis().GetBinLowEdge(ix)
              break
      #print "Cut for %s @ rate %g: %g" % (what+obj, rate, cut)
      plot = makeEffHist(name, genArray, recoArrayS, cut, options.xmax, logxbins=options.logxbins)
      eff = plot.GetEfficiency(plot.FindFixBin(options.genht))
      #print "Eff at %g: %g" % (options.genht, eff)
      return eff, plot, cut
    rate = 50
    while True:
      eff, plot, cut = effForRate(rate)
      if not eff: break
      if eff < plotparam:
          if rate >= 30e3: 
              eff = None; break
          rate = min(rate*5, 30e3)
      else:
          break
    if not eff: return (None, None)
    #print "Upper bound: eff %g at rate %g" % (eff, rate) 
    maxrate = rate; rate = rate / 5;
    while True:
      eff, plot, cut = effForRate(rate)
      if not eff: break
      if eff > plotparam:
          rate /= 2; 
      else:
          break
    if not eff: return (None, None)                 
    #print "Lower bound: eff %g at rate %g" % (eff, rate) 
    minrate = rate
    while True:
      rate = sqrt(maxrate * minrate)
      eff, plot, cut = effForRate(rate)
      if not eff: 
          #print "Bisection failed?"
          break
      if eff < plotparam:
          #print "New lower bound: eff %g at rate %g" % (eff, rate)
          minrate = rate;
      else:
          #print "New upper bound: eff %g at rate %g" % (eff, rate)
          maxrate = rate
      if maxrate/minrate < 1.1: 
          break
    #platprogress.done()
    if not eff: return (None, None)
    if rate < 10:    ratestr = "%.1fkHz" % rate
    elif rate < 500: ratestr = "%.0fkHz" % rate
    else:            ratestr = "%.1fMHz" % (rate/1000)
    label = "%s @ %s" % (name, ratestr)
    _cache[_key] = (plot,label)
    return (plot,label)

def makePlatRocPlot(signal, background, what, obj, ptcut, jecs, plotparam, _cache={}):
    _key = (id(signal),id(background),what,obj,str(ptcut),str(plotparam))
    if _key in _cache: 
        #print "  retrieved platroc for %s, %s, %s from cache" % (what, obj, plotparam)
        return _cache[_key]
    # ok there we go
    #platprogress = _progress("  making platroc for %s, %s, %s" % (what, obj, plotparam))
    if isJetMet:
        recoArrayB = makeCorrArray(background, what, obj, ptcut, options.eta, jecs)
        if not recoArrayB: return (None,None)
        recoArrayS = makeCorrArray(signal, what, obj, ptcut, options.eta, jecs)
        if not recoArrayS: return (None,None)
    else:
        recoArrayB = makeRecoLepArray(background, what, obj, ptcut, options.etas, lepid=options.lepid_, lepiso=options.lepiso_) 
        if not recoArrayB: return (None,None)
        #print "Lepeff %s %s > %s ptCut %s gen %s ptCut %s" % (obj, what, plotparam, ptcut, genObjName, options.genleppt)
        recoArrayS = makeMatchedRecoLepArray(signal, what, 
                      obj, ptcut, options.etas, options.lepid_, options.lepiso_,
                      genObjName, options.genleppt, options.etas, options.genprompt, options.gendr)
        if not recoArrayS: return (None,None)
    rateplot = makeCumulativeHTEff(name, recoArrayB, options.xmax)
    def platForRate(rate):
      cut = None
      for ix in range(1,rateplot.GetNbinsX()+1):
          if rateplot.GetBinContent(ix) <= rate:
              cut = rateplot.GetXaxis().GetBinLowEdge(ix)
              break
      if cut is None:
          return None
      #print "Cut for %s @ rate %g: %g" % (what+obj, rate, cut)
      plot = makeEffHist(name, genArray, recoArrayS, cut, options.xmax, logxbins=options.logxbins)
      hist = plot.GetTotalHistogram(); xaxis = hist.GetXaxis()
      for i in range(2,hist.GetNbinsX()+1):
          if plot.GetEfficiency(i) > plotparam:
              graph = plot.CreateGraph()
              xmin, xmax = xaxis.GetBinLowEdge(i-1), xaxis.GetBinCenter(i)
              if graph.Eval(xmin,ROOT.nullptr,"S") > plotparam: 
                  #print "ERROR xmin for %s @ %g (%g): i = %d" % (what+obj, rate, cut, i)
                  pass
              elif graph.Eval(xmax,ROOT.nullptr,"S") < plotparam:
                  #print "ERROR xmax for %s @ %g (%g): i = %d" % (what+obj, rate, cut, i)
                  pass
              else:
                  xmid = 0.5*(xmax+xmin)
                  while abs(xmax-xmin) > 0.05*(abs(xmax)+abs(xmin)):
                      if graph.Eval(xmid,ROOT.nullptr,"S") < plotparam:
                          xmin = xmid
                      else:
                          xmax = xmid
                      xmid = 0.5*(xmax+xmin)
                  #print "Plateau at %g eff for %s @ %g: %g" % (plotparam, what+obj, rate, xmid)
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
    #platprogress.done(" done, with %d points" % len(points))
    if not points: return (None,None)
    plot = ROOT.TGraph(len(points))
    for i,(x,y) in enumerate(points):
      plot.SetPoint(i,x,y)
    label = name
    _cache[_key] = (plot,label)
    return (plot,label)

print("Plotting for %s (%s)" % (options.var, options.varlabel))
for plotkind in options.plots.split(","):
  progress = _progress("Make plot %s: \n" % plotkind)
  if plotkind not in ("rate","inclrate"):
    if isJetMet:
        genArray = makeGenArray(signal, what, options.pt, options.eta)
    else:
        genArray = False; genObjName = None
        for objset,things in whats:
          if options.what and (objset not in options.what.split(",")): continue
          if options.what_reg:
              if not any(re.match(p+"$",objset) for p in options.what_reg.split(",")): 
                  continue
          for name,obj,col,msty,msiz in things:
            obj = obj.rstrip("$")
            if name == "Gen":
                if obj == genObjName: continue
                if genObjName != None: raise RuntimeError("Duplicate gen definition: %s, %s" % (genObjName, obj))
                genObjName = obj
                genArray = makeGenLepArray(signal, options.xvar, obj, options.genleppt, options.etas)
  if plotkind in ("isorate", "incleff","inclrate","lepeff") :
      plotparams = list(map(float, options.rate.split(",")))
  elif plotkind in ("platroc", "plateff"):
      plotparams = list(map(float, options.eff.split(",")))
  else:
      plotparams = [None]
  for plotparam in plotparams:
      for objset,things in whats:
          if options.what and (objset not in options.what.split(",")): continue
          if options.what_reg:
              if not any(re.match(p+"$",objset) for p in options.what_reg.split(",")): 
                  continue
          plots = []
          for name,obj,col,msty,msiz in things:
              obj = obj.rstrip("$")
              if (name == "Gen") or ("GenAcc" in obj): continue
              if options.var.startswith("met") or obj.startswith("Ref") or not isJetMet:
                  jecs = None
              else:
                  jecdirname = obj+"Jets"+( "_"+options.jecMethod if options.jecMethod else "")
                  jecdir = jecfile.GetDirectory(jecdirname)
                  if not jecdir: 
                      print("Missing JECs "+jecdirname+" in "+options.jecs)
                      continue
                  jecs = ROOT.l1tpf.corrector(jecdir)
              label = name
              ptcut = options.pt
              if "RefTwoLayerJets" in obj: ptcut = 5
              if plotkind == "rate":
                  if isJetMet:
                      recoArrayB = makeCorrArray(background, what, obj, ptcut, options.eta, jecs)
                  else:
                      recoArrayB = makeRecoLepArray(background, what, obj, ptcut, options.etas, lepid=options.lepid_, lepiso=options.lepiso_) 
                  if not recoArrayB: continue
                  plot = makeCumulativeHTEff(name, recoArrayB, options.xmax)
              elif plotkind == "effc":
                  if isJetMet:
                      recoArrayS = makeCorrArray(signal, what, obj, ptcut, options.eta, jecs)
                  else:
                      recoArrayS = makeRecoLepArray(signal, what, obj, ptcut, options.etas, lepid=options.lepid_, lepiso=options.lepiso_)
                  if not recoArrayS: continue
                  plot = makeCumulativeHTEffGenCut(name, recoArrayS, genArray, options.genht, options.xmax, norm=1)
              elif plotkind == "lepeff":
                  if isJetMet: raise RuntimeError()
                  #print "Lepeff vs %s for %s %s > %s ptCut %s gen %s ptCut %s eta %s" % (options.xvar, obj, what, plotparam, ptcut, genObjName, options.genleppt, options.eta)
                  recoArrayS = makeMatchedRecoLepArray(signal, what, 
                                obj, ptcut, options.etas, options.lepid_, options.lepiso_,
                                genObjName, options.genleppt, options.etas, options.genprompt, options.gendr)
                  plot = makeEffHist(name, genArray, recoArrayS, plotparam, options.xmax, logxbins=options.logxbins)
              elif plotkind == "isorate":
                  targetrate = plotparam
                  if isJetMet:
                      recoArrayB = makeCorrArray(background, what, obj, ptcut, options.eta, jecs)
                      if not recoArrayB: continue
                      recoArrayS = makeCorrArray(signal, what, obj, ptcut, options.eta, jecs)
                      if not recoArrayS: continue
                  else:
                      recoArrayB = makeRecoLepArray(background, what, obj, ptcut, options.etas, lepid=options.lepid_, lepiso=options.lepiso_) 
                      if not recoArrayB: continue
                      #print "Lepeff %s %s > %s ptCut %s gen %s ptCut %s" % (obj, what, plotparam, ptcut, genObjName, options.genleppt)
                      recoArrayS = makeMatchedRecoLepArray(signal, what, 
                                    obj, ptcut, options.etas, options.lepid_, options.lepiso_,
                                    genObjName, options.genleppt, options.etas, options.genprompt, options.gendr)
                      if not recoArrayS: continue
                  rateplot = makeCumulativeHTEff(name, recoArrayB, options.xmax)
                  if not rateplot: continue
                  cut = 9999
                  for ix in range(1,rateplot.GetNbinsX()+1):
                      if rateplot.GetBinContent(ix) <= targetrate:
                          cut = rateplot.GetXaxis().GetBinLowEdge(ix)
                          break
                  plot = makeEffHist(name, genArray, recoArrayS, cut, options.xmax, logxbins=options.logxbins)
                  label = "%s(%s) > %.0f" % (options.varlabel, name,cut)
              elif plotkind == "roc":
                  if isJetMet:
                      recoArrayB = makeCorrArray(background, what, obj, ptcut, options.eta, jecs)
                      if not recoArrayB: continue
                      recoArrayS = makeCorrArray(signal,     what, obj, ptcut, options.eta, jecs)
                      if not recoArrayS: continue
                  else:
                      recoArrayB = makeRecoLepArray(background, what, obj, ptcut, options.etas, lepid=options.lepid_, lepiso=options.lepiso_) 
                      if not recoArrayB: continue
                      #print "Lepeff %s %s > %s ptCut %s gen %s ptCut %s" % (obj, what, plotparam, ptcut, genObjName, options.genleppt)
                      recoArrayS = makeMatchedRecoLepArray(signal, what, 
                                    obj, ptcut, options.etas, options.lepid_, options.lepiso_,
                                    genObjName, options.genleppt, options.etas, options.genprompt, options.gendr)
                      if not recoArrayS: continue
                  effsig  = makeCumulativeHTEffGenCut(name+"_s", recoArrayS, genArray, options.genht, options.xmax, norm=1)
                  ratebkg = makeCumulativeHTEff(name+"_b", recoArrayB, options.xmax)
                  if not effsig or not ratebkg: continue
                  plot = ROOT.makeROCFast(effsig,ratebkg)
                  msty = 0
              elif plotkind == "plateff":
                  (plot, label) = makePlatEffPlot(signal, background, what, obj, ptcut, jecs, plotparam)
              elif plotkind == "platroc":
                  (plot, label) = makePlatRocPlot(signal, background, what, obj, ptcut, jecs, plotparam)
              elif plotkind == "incleff":
                  if isJetMet:
                      recoArrayS = makeCorrArray(signal, what, obj, ptcut, options.eta, jecs)
                  else:
                      recoArrayS = makeMatchedRecoLepArray(signal, what, 
                                obj, ptcut, options.etas, options.lepid_, options.lepiso_,
                                genObjName, options.genleppt, options.etas, options.genprompt, options.gendr)

                  if not recoArrayS: continue
                  makeInclusiveEffGenHTCut(name, recoArrayS, plotparam, genArray, options.genht)
                  continue
              elif plotkind == "inclrate":
                  if isJetMet:
                    recoArrayB = makeCorrArray(background, what, obj, ptcut, options.eta, jecs)
                  else:
                    recoArrayB = makeRecoLepArray(background, what, obj, ptcut, options.etas, lepid=options.lepid_, lepiso=options.lepiso_) 
                  if not recoArrayB: continue
                  makeInclusiveEffRate(name, recoArrayB, plotparam)
                  continue
              else: raise RuntimeError
              if not plot: continue
              plot.SetLineWidth(3); plot.SetLineColor(col);  plot.SetMarkerColor(col)
              plot.SetMarkerStyle(msty); plot.SetMarkerSize(msiz)
              plots.append((label,plot))
          if not plots: 
              print("   nothing to plot!")
              continue
          plotter.SetLogy(False)
          plotstyle = "PCX"
          ROOT.gStyle.SetErrorX(0)
          if plotkind == "rate":
              plotter.SetLogy(True)
              frame = ROOT.TH1D("",";L1 %s cut (%s); Minbias rate @ PU200 [kHz]" % (options.varlabel, qualif), 100, 0, options.xmax)
              frame.GetYaxis().SetDecimals(True)
              frame.GetXaxis().SetNdivisions(505)
              frame.GetYaxis().SetRangeUser(options.ymin if options.ymin else 0.5, 100e3)
              leg = ROOT.TLegend(0.56,0.93,0.93,0.93-0.055*len(things))
              plotname = '%s%s-%s_eta%s_pt%d' % (options.var, plotkind, objset, options.eta, options.pt)
          elif plotkind == "effc":
              gentext, genpost = "iciency", "" 
              if options.genht > 0:
                  gentext = " (Gen %s > %s)" % (options.varlabel, options.genht)
                  genpost = "_gen%.0f" % (options.genht)
              frame = ROOT.TH1D("",";L1 %s thresh (%s); Eff%s" % (options.varlabel, qualif, gentext), 100, 0, options.xmax)
              frame.GetYaxis().SetRangeUser(options.ymin if options.ymin else 0.0, 1.05)
              leg = ROOT.TLegend(0.6,0.93,0.93,0.93-0.055*len(things))
              plotname = '%s%s-%s_eta%s_pt%d%s' % (options.var, plotkind, objset, options.eta, options.pt, genpost)
          elif plotkind == "isorate":
              frame = ROOT.TH1D("",";Gen %s (%s); Eff (L1 rate %.0f kHz)" % (options.varlabel, qualif, plotparam), 100, 0, options.xmax)
              frame.GetYaxis().SetRangeUser(options.ymin if options.ymin else 0.0, 1.05)
              frame.GetYaxis().SetDecimals(True)
              leg = ROOT.TLegend(0.50,0.18,0.93,0.18+0.055*len(things))
              plotname = '%s%s-%s_eta%s_pt%d_%.0fkHz' % (options.var, plotkind, objset, options.eta, options.pt, plotparam)
          elif plotkind == "plateff":
              frame = ROOT.TH1D("",";Gen %s (%s); Eff" % (options.varlabel, qualif), 100, 0, options.xmax)
              frame.GetYaxis().SetRangeUser(options.ymin if options.ymin else 0.0, 1.05)
              frame.GetYaxis().SetDecimals(True)
              leg = ROOT.TLegend(0.30,0.18,0.93,0.18+0.055*len(things))
              plotname = '%s%s-%s_eta%s_pt%d_eff%.2fat%g' % (options.var, plotkind, objset, options.eta, options.pt, plotparam, options.genht)
          elif plotkind == "roc":
              plotter.SetLogy(True)
              gentext, genpost = "iciency", "" 
              if options.genht > 0:
                  gentext = " (Gen %s > %s)" % (options.varlabel, options.genht)
                  genpost = "_gen%.0f" % (options.genht)
              frame = ROOT.TH1D("",";Eff%s; Minbias rate @ PU200 [kHz]" % (gentext), 100, 0, 1)
              frame.GetYaxis().SetDecimals(True)
              frame.GetXaxis().SetNdivisions(505)
              frame.GetYaxis().SetRangeUser(options.ymin if options.ymin else 0.5, 100e3)
              leg = ROOT.TLegend(0.2,0.93,0.55,0.93-0.055*len(things))
              plotname = '%s%s-%s_eta%s_pt%d%s' % (options.var, plotkind, objset, options.eta, options.pt, genpost)
              plotstyle = "PLX"
          elif plotkind == "platroc":
              plotter.SetLogy(True)
              xtit = "Gen threshold at %g%% eff" % (plotparam*100)
              frame = ROOT.TH1D("",";%s; Minbias rate @ PU200 [kHz]" % (xtit), 100, 0, options.xmax)
              frame.GetYaxis().SetDecimals(True)
              #frame.GetXaxis().SetNdivisions(505)
              frame.GetYaxis().SetRangeUser(options.ymin if options.ymin else 0.5, 100e3)
              leg = ROOT.TLegend(0.6,0.93,0.93,0.93-0.055*len(things))
              plotname = '%s%s-%s_eta%s_pt%d_eff%.3f' % (options.var, plotkind, objset, options.eta, options.pt, plotparam)
          elif plotkind == "lepeff":
              gentext, genpost = "iciency", "" 
              if options.genprompt != 2: genpost = "_id%d" % options.genprompt
              xmin = options.etamin*0.8 if "eta" in options.xvar else 0
              frame = ROOT.TH1D("",";Gen %s; Eff%s" % (options.varlabel, gentext), 100, xmin, options.xmax)
              frame.GetYaxis().SetRangeUser(options.ymin if options.ymin else 0.0, 1.05)
              leg = ROOT.TLegend(0.57,0.18,0.93,0.18+0.055*len(things))
              streta = ("%s_%s" % options.etas) if options.etamin > 0 else str(options.eta)
              plotname = '%s%s-%s_eta%s_genpt%s%s_recpt%s_%s%s' % (options.xvar, plotkind, objset, streta, options.genleppt, genpost, options.pt, options.var, plotparam)
              if "eta" in options.xvar:
                  ROOT.gStyle.SetErrorX(0.5)
                  plotstyle = "PE"
          frame.Draw()
          if plotkind in ("rate","roc","platroc"):
              line = ROOT.TLine()
              line.SetLineStyle(7)
              for y in 40e3, 100, 10:
                  line.DrawLine(frame.GetXaxis().GetXmin(),y,frame.GetXaxis().GetXmax(),y)
          for n,p in plots: 
              if ("TH1" not in p.ClassName()): 
                  plotstyle = plotstyle.lstrip("R").replace("X","")
              p.Draw(plotstyle+" SAME")
          for n,p in plots: 
              leg.AddEntry(p, n, "L" if plotkind in ("rate","roc") else "LP")
          leg.Draw()
          plotter.decorations()
          #printprogress = _progress("  Printing plot %s%s (plot)" % (plotname, ("_"+options.label) if options.label else ""))
          plotter.Print('%s%s' % (plotname, ("_"+options.label) if options.label else ""))
          fout = ROOT.TFile.Open('%s/%s%s.root' % (odir, plotname, ("_"+options.label) if options.label else ""), "RECREATE")
          fout.WriteTObject(frame,"frame")
          for n,p in plots: 
              p.SetTitle(n)
              fout.WriteTObject(p)
          fout.Close()
          del frame
          #printprogress.done()
  progress.done("  Done plot %s" % plotkind)
  print("")



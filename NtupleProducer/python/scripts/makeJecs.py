import os, sys
import ROOT
from array import array
from timeit import default_timer
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kWarning
from FastPUPPI.NtupleProducer.scripts.respPlots import doRespPt

class _progress:
    def __init__(self,statement):
        self._t0 = default_timer()
        print(statement, end=' ')
        sys.stdout.flush()
    def done(self, statement="done"):
        print("%s (%.2f s)" % (statement, default_timer() - self._t0))

def _list2graph(datapoints):
    ret = ROOT.TGraph(len(datapoints))
    for i,(x,y) in enumerate(datapoints):
        ret.SetPoint(i,x,y)
    return ret

if __name__ == "__main__":
    from optparse import OptionParser
    parser = OptionParser("%(prog) infile [ src [ dst ] ]")
    parser.add_option("-m", dest="methods",  default="fitlin")
    parser.add_option("--etabins", dest="etabins", default="0,1.3,1.7,1.9,2.1,2.4,2.8,3.0,3.3,3.6,4.0,4.8")
    parser.add_option("--yr","--yrange", dest="yrange", nargs=2, default=(0.65,1.3), type="float")
    parser.add_option("--fr","--fitrange", dest="fitrange", nargs=2, default=(20,999), type="float")
    parser.add_option("--cp","--closure-plots", dest="closurePlots", default=None, type="string")
    parser.add_option("--raw-in-closure-plots", dest="rawInClosurePlots", default=False, action="store_true", help="Include raw in closure plots")
    parser.add_option("--ptmin", dest="ptmin", nargs=1, default=0, type="float")
    parser.add_option("--ptmax", dest="ptmax", nargs=1, default=0, type="float")
    parser.add_option("-A","--auto", dest="auto", default=False, action="store_true", help="make JECs for all inputs")
    parser.add_option("-g","--gauss", dest="gauss", default=False, action="store_true", help="make also gaussian estimates")
    parser.add_option("-o","--out", dest="rootfile", default="jecs.root")
    options, args = parser.parse_args()

    tree = ROOT.TChain("Events")
    filenames = []
    for a in args[:]:
        if a.endswith(".root"): 
            filenames.append(a)
            tree.Add(a)
            args.remove(a)
    whats = args
    if options.auto:
        for b in tree.GetListOfBranches():
            bname = str(b.GetName())
            if bname.endswith("Jets_genpt"):
                whats.append(bname.replace("_genpt",""))
                print("Got branch "+whats[-1])

    ROOT.gSystem.Load("libL1TriggerPhase2L1ParticleFlow")
    ROOT.gInterpreter.ProcessLine('#include "L1Trigger/Phase2L1ParticleFlow/interface/corrector.h"')

    if options.closurePlots:
        os.system("mkdir -p "+options.closurePlots)
        os.system("cp %s/src/FastPUPPI/NtupleProducer/python/display/index.php %s/" % (os.environ['CMSSW_BASE'], options.closurePlots));
        ROOT.gROOT.ProcessLine(".x %s/src/FastPUPPI/NtupleProducer/python/display/tdrstyle.cc" % os.environ['CMSSW_BASE']);
        ROOT.gStyle.SetOptStat(False)
        ROOT.gStyle.SetOptFit(False)
        ROOT.gStyle.SetErrorX(0.5)
        c1 = ROOT.TCanvas("c1","c1")

    etas = list(map(float,options.etabins.split(",")))
    indextemplate = ROOT.TH1F("indextemplate", "", len(etas)-1, array('f', etas))
    indextemplate.SetDirectory(ROOT.nullptr)

    methods = options.methods.split(",")
    ptCorrs = {}
    for w in whats:
        for m in methods:
            h = ROOT.l1tpf.corrector(indextemplate)
            ptCorrs[(m,w)] = h
    for w in whats:
        print("Processing "+w)
        expr = w+"_pt"; mcpt = w + "_genpt"
        for ieta in range(1,len(etas)):
            ptForResp = [20,25,30,35,40,45,50,55,60,70,80,90,100,120,140,160,200,250,320,400]
            if etas[ieta-1] > 3: ptForResp = [25,35,45,60,80,100,120,150,190,250]
            recoPtMin = 5 if etas[ieta-1] < 2.4 else 10
            cut = "{w}_gendr < 0.2 && {mcpt} > 20 && {etamin} < abs({w}_eta) && abs({w}_eta) < {etamax} && {expr} > {recoPtMin}".format(w=w, mcpt=mcpt, expr=expr, etamin=etas[ieta-1],etamax=etas[ieta], recoPtMin=recoPtMin)
            perf = _progress( "Making response for %s (%.3f < |eta| < %.3f)..." % (w,etas[ieta-1],etas[ieta]))
            prof, _ = doRespPt("jet",tree,w,expr,cut,mcpt=mcpt,xpt=mcpt,fitrange=options.fitrange,fitfunc="none",ptbins=ptForResp,keepGraph=options.closurePlots)
            perf.done()
            if prof and options.gauss: prof = getattr(prof,'gaus',None) 
            if not prof: 
                print("No response plot for "+cut)
                continue
            responses = [ prof.GetY()[i] for i in range(prof.GetN()) ]
            if len([r for r in responses if r > 0.2]) < 3:
                print("Response too low for this bin. skipping")
                continue
            for m in methods:
                h = ptCorrs[(m,w)]
                if m in ("fitlin",):
                    respformula = "1/x++1"
                    tf1 = ROOT.TF1(w+m+"_f1", respformula, 0, ptForResp[-1])
                    prof.Fit(tf1, "WNQ0C EX0", "", options.fitrange[0], options.fitrange[1])
                    (scale,offs) = tf1.GetParameter(1), tf1.GetParameter(0) 
                    xmin, xmax, npoints = recoPtMin, min(500.,(400.-offs)/scale), 200
                    corrgraph = ROOT.TGraph(npoints)
                    for i in range(npoints):
                        pt = xmin + i*float((xmax-xmin)/(npoints-1))
                        corrgraph.SetPoint(i, pt, (pt-offs)/scale)
                    h.setGraph(corrgraph, ieta-1);
            if options.closurePlots:
                cols = [ ROOT.kRed, ROOT.kGreen+2, ROOT.kBlue, ROOT.kMagenta+1, ROOT.kOrange+7, ROOT.kCyan+1, ROOT.kGray+2, ROOT.kViolet+5, ROOT.kSpring+5, ROOT.kAzure+1, ROOT.kPink+7, ROOT.kOrange+3, ROOT.kBlue+3, ROOT.kMagenta+3, ROOT.kRed+2, ] 
                perf = _progress( "   making closure plots...")
                plots = []; 
                ptForClosure = [10,20,25,30,35,40,45,50,55,60,70,80,90,100,120,140,160,200,250] 
                if etas[ieta-1] > 3: ptForClosure = [20,30,40,50,60,80,100,135,200]
                if options.rawInClosurePlots:
                    rresp, _ = doRespPt("jet",None,w+"_"+m, None, None, fitfunc="none", ptbins=ptForResp, fromGraph=prof.graph)
                    if rresp: plots.append(("raw",rresp,rresp.gaus))
                    cols = [ ROOT.kGray+2 ] + [ c for c in cols if c != ROOT.kGray+2 ]
                ptgen, ratio = prof.graph.GetX(), prof.graph.GetY()
                for m in methods:
                    corrgraph = ptCorrs[(m,w)].getGraph(ieta-1)
                    datapoints = []
                    for i in range(prof.graph.GetN()):
                        ptrec = ptgen[i]*ratio[i]
                        if ptrec < recoPtMin: continue
                        datapoints.append((ptgen[i], corrgraph.Eval(ptrec)/ptgen[i]))
                    cresp, _ = doRespPt("jet",None,w+"_"+m, None, None, fitfunc="none", ptbins=ptForResp, fromGraph=_list2graph(datapoints))
                    if cresp:
                        plots.append((m,cresp,cresp.gaus))
                if not plots: continue
                frame = ROOT.TH1F("%s_clos_eta_%.2f_%.2f"%(w,etas[ieta-1],etas[ieta]),"stk;p_{T}^{ref};<p_{T}^{corr}/p_{T}^{ref}>",100,0.0,ptForClosure[-1])
                frame.GetYaxis().SetDecimals(True)
                frame.GetYaxis().SetRangeUser(options.yrange[0], options.yrange[1])
                frame.Draw()
                line = ROOT.TLine()
                line.SetLineStyle(7)
                line.DrawLine(0.0,1.0,ptForClosure[-1],1.0)
                leg = ROOT.TLegend(0.3,0.99,0.95,0.99-0.042*(len(plots)+0.5))
                leg.SetTextSize(0.04);
                leg.SetNColumns(2)
                for isample,(pname,histm,histg) in enumerate(plots):
                    for hist,gauss in [(histg,True),(histm,False)]:
                        hist.SetLineColor(cols[isample])
                        hist.SetMarkerColor(cols[isample])
                        hist.SetMarkerSize(1.1 if gauss else 0.85)
                        hist.SetMarkerStyle(24 if gauss else 20)
                        hist.SetLineWidth(1 if gauss else 2)
                        hist.SetLineStyle(7 if gauss else 1)
                        hist.Draw("P SAME")
                    leg.AddEntry(histm, pname, "LP")
                    leg.AddEntry(histg, "(gauss)", "LP")
                leg.Draw()
                c1.Print("%s/%s.png" % (options.closurePlots, frame.GetName()))
                del frame
                perf.done()
    tfout = ROOT.TFile.Open(options.rootfile, "RECREATE");
    print("Saving to output file: ")
    for (m,w),h in ptCorrs.items():
        dirname = ("%s_%s" % (w,m)) if len(methods) > 1 else w 
        h.writeToFile(tfout.mkdir(dirname))
        print("  - ",dirname)
    tfout.Close()
    print("Saved responses in %s" % options.rootfile)


import os, re, sys
from collections import defaultdict
import ROOT
from FastPUPPI.NtupleProducer.plotTemplate import plotTemplate


V="v0"; VT="v0"
PLOTDIR="plots/106X/from104X/"+V+"/tdr/"+VT+"/"


PARTICLES=["pion","pizero"]
ETAS=[(0,1.3),(1.7,2.5)]
if False:
    NTUPLES = " ".join("respTupleNew_Particle"+G+"_PU0."+V+".root" for G in ("Gun","GunPt0p5To5","GunPt80To300"))
    for P in PARTICLES:
        for E in ETAS:
            os.system("python scripts/respPlots.py %s %s -w l1pf -p %s --eta %g %g -g --no-eta --ptmax 100 --corr-resol=none" % (NTUPLES,PLOTDIR,P,E[0],E[1]))
    sys.exit(0)

plotter = plotTemplate(PLOTDIR)
markers = [ 21, 20,  25, 24 ]
sizes   = [ 1.5, 1.35, 1.4, 1.2 ]
colors  = [ ROOT.kGreen, ROOT.kGreen+2, ROOT.kAzure+1, ROOT.kBlue+2, ]
for (R,G) in ("",""), ("_res","_gauss"):
    tfiles = {}
    graphs = []; labels = []
    i = 0;
    for P in PARTICLES:
        for E in ETAS:
            fname = "%s_eta_%02d_%02d%s%s-l1pf_pt02.root" % (P,E[0]*10,E[1]*10,R,G)
            tfiles[fname] = ROOT.TFile.Open(PLOTDIR+"/"+fname)
            g = tfiles[fname].Get("PF%s%s"%(G,R))
            g.SetMarkerStyle(markers[i])
            g.SetMarkerSize(sizes[i])
            g.SetMarkerColor(colors[i])
            graphs.append(g)
            labels.append( "#pi^{%s}, %.1f < |#eta| < %.1f" % ("0" if P == "pizero" else "#pm", E[0], E[1]) )
            i += 1
    frame = ROOT.TH1F("frame","frame",100,0,100)
    frame.GetXaxis().SetTitle("p_{T} (GeV)")
    if R == "_res":
        frame.GetYaxis().SetTitle("#sigma(p_{T})/p_{T}^{gen}")
        frame.GetYaxis().SetRangeUser(0, 0.35)
        frame.GetYaxis().SetNdivisions(1005)
    else:
        frame.GetYaxis().SetTitle("median p_{T}/p_{T}^{gen}")
        frame.GetYaxis().SetRangeUser(0.5, 1.1)
    frame.GetYaxis().SetDecimals(True)
    frame.Draw()

    legy = 0.91 if R == "_res" else 0.40
    leg = ROOT.TLegend(0.45,legy,0.92,legy-0.06*len(labels))
    #leg.SetNColumns(2)
    leg.SetTextSize(0.05)
    for n,p in zip(labels,graphs): 
        p.Draw("PX SAME")
        leg.AddEntry(p, n, "P")
    leg.Draw()
    
    plotter.decorations(pu=0)
    plotter.Print("pions_stack%s%s" % (R,G))

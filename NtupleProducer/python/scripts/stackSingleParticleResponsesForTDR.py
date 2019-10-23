import os, re, sys
import ROOT
from FastPUPPI.NtupleProducer.plotTemplate import plotTemplate


V="v0.1_TDR"; VT="v0"
PLOTDIR="plots/106X/from104X/"+V+"/tdr/"+VT+"/"

PARTICLES=["pion","pizero","klong"]
if "pions" in sys.argv: 
    PARTICLES=PARTICLES[:2]
    pions=True
else:
    pions=False

LABELS=dict(pion="#pi^{#pm}", pizero="#pi^{0}", klong="K^{0}_{L}")
ETAS=[(0,1.3),(1.7,2.5)]
if "replot" in sys.argv:
    NTUPLES = " ".join("respTupleNew_Particle"+G+"_PU0."+V+".root" for G in ("Gun","GunPt0p5To5","GunPt80To300"))
    for P in PARTICLES:
        for E in ETAS:
            os.system("python scripts/respPlots.py %s %s -w l1pf -p %s --eta %g %g -g --no-eta --ptmax 100 --corr-resol=none" % (NTUPLES,PLOTDIR,P,E[0],E[1]))
    sys.exit(0)

plotter = plotTemplate(PLOTDIR)
markers = [ 21, 25, 20, 24, 33, 27 ]
sizes   = [ 1.65, 1.40, 1.55, 1.25, 1.9, 1.6 ]
colors  = [ ROOT.kGreen, ROOT.kGreen+3, ROOT.kAzure+1, ROOT.kBlue+2, ROOT.kRed+1, ROOT.kRed+1 ]
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
            g.SetLineWidth(2)
            graphs.append(g)
            labels.append( "%s, %.1f < |#eta| < %.1f" % (LABELS[P], E[0], E[1]) )
            i += 1
    frame = ROOT.TH1F("frame","frame",100,0,100)
    frame.GetXaxis().SetTitle("p_{T}^{gen} (GeV)")
    if R == "_res":
        frame.GetYaxis().SetTitle("#sigma(p_{T})/p_{T}^{gen}")
        frame.GetYaxis().SetRangeUser(0, 0.35 if pions else 0.55)
        if pions: frame.GetYaxis().SetNdivisions(1005)
    else:
        frame.GetYaxis().SetTitle("median p_{T}/p_{T}^{gen}")
        frame.GetYaxis().SetRangeUser(0.5, 1.1)
    frame.GetYaxis().SetDecimals(True)
    frame.Draw()

    legysize = 0.065*len(labels)
    legy = 0.91 if R == "_res" else 0.16+legysize
    leg = ROOT.TLegend(0.45,legy,0.92,legy-legysize)
    #leg.SetNColumns(2)
    leg.SetTextSize(0.05)
    for n,p in zip(labels,graphs): 
        p.Draw("PX SAME")
        leg.AddEntry(p, n, "P")
    leg.Draw()
    
    plotter.decorations(pu=0)
    if pions: 
        plotter.Print("pions_stack%s%s" % (R,G))
        print "Made pions"
    else:
        plotter.Print("hadrons_stack%s%s" % (R,G))
        print "Made hadrons"

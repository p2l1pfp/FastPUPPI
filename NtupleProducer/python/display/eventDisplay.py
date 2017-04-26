#!/usr/bin/env python
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(True)
ROOT.gSystem.Load("libFWCoreFWLite.so");
ROOT.gSystem.Load("libDataFormatsFWLite.so");
ROOT.FWLiteEnabler.enable()

import os
from sys import argv
from math import *

from DataFormats.FWLite import Handle, Events

from FastPUPPI.NtupleProducer.display.drawers import *
from FastPUPPI.NtupleProducer.display.physobjlist import *

if len(argv) <= 1:
    print """
        usage: python eventDisplay.py file.root output_directory/    [event number]
           or: python eventDisplay.py file.root secondaryfile.root output_directory/
"""

if ".root" in argv[2]:
    print "Primary: %s, Secondary: %s" % (argv[1], argv[2])
    events = Events(options=Opts(files=[argv[1]], secFiles=[argv[2]]))
    del argv[2]
else:
    events = Events(argv[1])

out = argv[2]
if not os.path.isdir(out):
    os.system("mkdir -p "+out)
    os.system("cp index.php "+out+"/")

genj  = Handle("std::vector<reco::GenJet>")
genp  = Handle("std::vector<reco::GenParticle>")
tracks   = Handle("std::vector<reco::Track>")
vertices = Handle("std::vector<reco::Vertex>")
pfcands  = Handle("std::vector<reco::PFCandidate>")
pfclust  = Handle("std::vector<reco::PFCluster>")
l1pfp    = Handle("std::vector<l1tpf::Particle>")
puppiw   = Handle("edm::ValueMap<float>")

for iev,event in enumerate(events):
    if len(argv) > 3 and event.eventAuxiliary().event() != int(argv[3]): continue
    print "\nEvent %1d %5d %12d" % (
        event.eventAuxiliary().run(),
        event.eventAuxiliary().luminosityBlock(),
        event.eventAuxiliary().event())
 
    phystruth = []
    physobj  = []

    genps = read(event, "genParticles", genp, filter = lambda g : g.status() == 1 and abs(g.pdgId()) not in (12,14,16))
    phystruth.append(PhysObjList("genParticles", genps, drawGenCands, ["all"], printer = lambda p : "charge %+1d pdgId %+4d" % (p.charge(), p.pdgId())))

    genjs = read(event, "ak4GenJetsNoNu", genj, filter = lambda g : g.pt() > 20)
    phystruth.append(PhysObjList("genJets", genjs, drawGenJets, ["all"]))

    ## === RAW L1 objects (reverse-ordered by depth)===
    physobj.append(PhysObjList("L1HgcBH", read(event, "l1tPFHGCalBHProducerFromOfflineRechits:towers", l1pfp), drawHGCalB, ["l1","l1in","l1in+puppi"]))
    physobj.append(PhysObjList("L1Hcal", (read(event, "l1tPFHcalProducerFromOfflineRechits:towers", l1pfp) + 
                                          read(event, "l1tPFHFProducerFromOfflineRechits:towers", l1pfp)), drawHcal,       ["l1","l1in","l1in+puppi"]))
    physobj.append(PhysObjList("L1HgcFH", read(event, "l1tPFHGCalFHProducerFromOfflineRechits:towers", l1pfp), drawHGCalH, ["l1","l1in","l1in+puppi"]))
    physobj.append(PhysObjList("L1HgcEE", read(event, "l1tPFHGCalEEProducerFromOfflineRechits:towers", l1pfp), drawHGCalE, ["l1","l1in","l1in+puppi"]))
    physobj.append(PhysObjList("L1Ecal", read(event, "l1tPFEcalProducerFromOfflineRechits:towers", l1pfp), drawEcal,       ["l1","l1in","l1in+puppi"]))
    physobj.append(PhysObjList("L1Tk",    read(event, "l1tPFTkProducersFromOfflineTracksStrips", l1pfp), drawTracks,       ["l1","l1in","l1in+puppi"]))
    
    ## === CORRELATOR OUTPUTS ===
    physobj.append(PhysObjList("L1C_Calo",  read(event, "InfoOut:Calo",  pfcands), drawCaloCands,  ["l1","l1out"]))
    physobj.append(PhysObjList("L1C_TK",    read(event, "InfoOut:TK",    pfcands), drawTkCands,    ["l1","l1out"]))
    physobj.append(PhysObjList("L1C_PF",    read(event, "InfoOut:PF",    pfcands), drawPFCands,    ["l1","l1out"]))
    #physobj.append(PhysObjList("L1C_Puppi", read(event, "InfoOut:Puppi", pfcands), drawPuppiCands, ["l1","l1out","l1in+puppi"]))
    physobj.append(PhysObjList("L1C_Puppi", read(event, "InfoOut:Puppi", pfcands), drawPFCands,    ["l1in+puppi","l1puppi"], printer = lambda p : "pdgId %+4d" % p.pdgId()))

    ## === CORRELATOR OUTPUTS (INTEGER MATH) ===
    physobj.append(PhysObjList("L1IC_Calo",  read(event, "InfoOut:L1Calo",  pfcands), drawCaloCands,  ["l1out-int"]))
    physobj.append(PhysObjList("L1IC_TK",    read(event, "InfoOut:L1TK",    pfcands), drawTkCands,    ["l1out-int"]))
    physobj.append(PhysObjList("L1IC_PF",    read(event, "InfoOut:L1PF",    pfcands), drawPFCands,    ["l1out-int"]))
    physobj.append(PhysObjList("L1IC_Puppi", read(event, "InfoOut:L1Puppi", pfcands), drawPFCands,    ["l1puppi-int"], printer = lambda p : "pdgId %+4d" % p.pdgId()))

    evname = "r%d_ls%d_ev%d"  % (event.eventAuxiliary().run(), event.eventAuxiliary().luminosityBlock(), event.eventAuxiliary().event())
    zoomRadius = 0.4
    #for view in "l1", "l1in", "l1out","l1in+puppi":
    for view in "l1", "l1out", "l1in+puppi", "l1puppi", "l1out-int", "l1puppi-int":
        cMaker.clear()
        vlog = open(out+"/%s_%s.txt" % (evname, view), "w")
        theLegend = ROOT.TLegend(0.01, 0.15, 0.11, 0.75);
        theLegend.SetLineColor(0)
        for p in physobj + phystruth: 
            p.draw(view)
            p.write(view,vlog)
            p.addToLegend(view,theLegend)
        theLegend.Draw()
        cMaker.write(event, out+"/%s_%s.png" % (evname, view))
        vlog.close()
        os.system("mkdir -p "+out+"/%s_%s.dir" % (evname, view))
        os.system("cp index.php "+out+"/%s_%s.dir/" % (evname, view))
        for i,gj in enumerate(genjs):
            cMaker.clear()
            cMaker.zoom((gj.eta(), gj.phi()), zoomRadius*1.7)
            vlog = open(out+"/%s_%s.dir/jet%d.txt" % (evname,view, i+1), "w")
            for p in physobj + phystruth: 
                p.draw(view)
                p.writeZoom(view, (gj.eta(), gj.phi()), zoomRadius*1.7, zoomRadius, vlog)
            cMaker.write(event, out+"/%s_%s.dir/jet%d.png" % (evname,view, i+1) )
            vlog.close()


    if len(argv) > 3: break
    if iev >= 10: break


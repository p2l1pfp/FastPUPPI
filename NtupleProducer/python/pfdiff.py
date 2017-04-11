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
from PhysicsTools.HeppyCore.statistics.tree import Tree
from PhysicsTools.HeppyCore.utils.deltar import *

tree = Tree("t","t",defaultFloatType="F") 
events = Events(argv[1])

pf1   = Handle("std::vector<reco::PFCandidate>")
pf2   = Handle("std::vector<reco::PFCandidate>")
genj  = Handle("std::vector<reco::GenJet>")
genp  = Handle("std::vector<reco::GenParticle>")

for iev,event in enumerate(events):
    print "\nEvent %1d %5d %12d" % (
        event.eventAuxiliary().run(),
        event.eventAuxiliary().luminosityBlock(),
        event.eventAuxiliary().event())
 
    if len(argv) == 2:
        event.getByLabel("ak4GenJetsNoNu", genj)
        allGenJ = [ g for g in genj.product() ]
        event.getByLabel("genParticles", genp)
        genLeptons = [ p for p in genp.product() if p.status() == 1 and abs(p.pdgId()) in (11,13) and (p.isPromptFinalState() or p.isDirectPromptTauDecayProductFinalState()) and p.pt() > 5 ]
        genTau     = [ p for p in genp.product() if abs(p.pdgId()) == 15 and p.isLastCopy() and p.isPromptDecayed() and p.pt() > 5 ]
        genGamma   = [ p for p in genp.product() if p.status() == 1 and abs(p.pdgId()) == 22 and p.isPromptFinalState() and p.pt() > 10 ]
        genItems = genLeptons + genTau + genGamma
        for j in allGenJ:
            incone = [ i for i in genItems if deltaR(i,j) < 0.4 ]
            if len(incone) == 0 or sum(i.pt() for i in incone) < 0.25*j.pt():
                j.mcId = 1
            elif len(incone) == 1 and j.pt() < 2*incone[0].pt():
                j.mcId = abs(incone[0].pdgId())
            else:
                j.mcId = 99
        allGenJ.sort(key = lambda j : - j.pt())
        print "Gen jets in the event:"
        for j in allGenJ:
            if abs(j.eta()) > 5: continue
            if j.mcId in (0, 1, 99) and j.pt() < 5: continue
            print "   pt %7.2f eta %+5.2f phi %+5.2f   id %d" % (j.pt(), j.eta(), j.phi(), j.mcId) 
        print ""
    else:
        print "Comparing %s to %s" % (argv[2], argv[3])
        event.getByLabel(argv[2], pf1)
        event.getByLabel(argv[3], pf2)
        o1 = [ o for o in sorted(pf1.product(), key = lambda p : -p.pt()) ]
        o2 = [ o for o in sorted(pf2.product(), key = lambda p : -p.pt()) ]
        match = matchObjectCollection3(o1, o2, deltaRMax=0.3)
        for p2 in o2: p2.seen = False
        for p1 in o1:
            print "pdgId %+4d pt %7.2f eta %+5.2f phi %+5.2f " % (p1.pdgId(), p1.pt(), p1.eta(), p1.phi()), 
            m1 = match[p1]
            if m1: 
                print " -> pdgId %+4d pt %7.2f eta %+5.2f phi %+5.2f " % (m1.pdgId(), m1.pt(), m1.eta(), m1.phi()),
                if m1.pdgId() == p1.pdgId() and abs(p1.pt()-m1.pt())/max(m1.pt()+p1.pt(),2) < 0.005 and deltaR2(p1,m1) < 0.001:
                    print "  [ same ]"
                elif m1.pdgId() == p1.pdgId() and abs(p1.pt()-m1.pt())/max(m1.pt()+p1.pt(),2) < 0.025 and deltaR2(p1,m1) < 0.03:
                    print "  [ Similar ]"
                else:
                    print "  [ DIFFERENT ]"
                m1.seen = True
            else:
                print " -> ..."
        for p2 in o2:
            if p2.seen: continue
            print "                                        ... -> pdgId %+4d pt %7.2f eta %+5.2f phi %+5.2f " % (p2.pdgId(), p2.pt(), p2.eta(), p2.phi()) 
    if len(argv) > 4 and iev < int(argv[4]):
        print "\n\n"; continue
    exit()

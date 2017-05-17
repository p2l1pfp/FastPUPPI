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
from PhysicsTools.HeppyCore.utils.deltar import *

from optparse import OptionParser
parser = OptionParser("%(prog) infile [ src [ dst ] ]")
parser.add_option("-n", "--events", type=int, nargs=1, default=3, help="Number of events to consider")
parser.add_option("-E", "--select-events", type=str, nargs=1, action="append", default=[], help="Specific event to select (run:lumi:event)")
options, args = parser.parse_args()
while len(args) < 3: args.append("")

events = Events(args[0])

ihpf = Handle("std::vector<l1tpf::Particle>")
hpf  = Handle("std::vector<reco::PFCandidate>")
genj  = Handle("std::vector<reco::GenJet>")
genp  = Handle("std::vector<reco::GenParticle>")

for iev,event in enumerate(events):
    if iev >= options.events: break
    idev = "%d:%d:%d" % ( event.eventAuxiliary().run(), event.eventAuxiliary().luminosityBlock(), event.eventAuxiliary().event())
    if options.select_events:
       if idev not in options.select_events: continue
    print "Event %s" % idev
 
    event.getByLabel("ak4GenJetsNoNu", genj)
    allGenJ = [ g for g in genj.product() ]
    event.getByLabel("genParticles", genp)
    genLeptons = [ p for p in genp.product() if p.status() == 1 and abs(p.pdgId()) in (11,13,211) and (p.isPromptFinalState() or p.isDirectPromptTauDecayProductFinalState()) and p.pt() > 5 ]
    genTau     = [ p for p in genp.product() if abs(p.pdgId()) == 15 and p.isLastCopy() and p.isPromptDecayed() and p.pt() > 5 ]
    genGamma   = [ p for p in genp.product() if p.status() == 1 and abs(p.pdgId()) == 22 and p.isPromptFinalState() and p.pt() > 5 ]
    genItems = genLeptons + genTau + genGamma
    for j in allGenJ:
        incone = [ i for i in genItems if deltaR(i,j) < 0.4 ]
        if len(incone) == 0 or sum(i.pt() for i in incone) < 0.25*j.pt():
            j.mcId = 1
        elif len(incone) == 1 and j.pt() < 2*incone[0].pt():
            j.mcId = abs(incone[0].pdgId())
        else:
            j.mcId = 99
        if abs(j.eta()) > 5: continue
        #daus = [ j.daughter(i) for i in xrange(j.numberOfDaughters()) ]
        #chall = sum(d.pt() for d in daus if d.charge() != 0)
        #chacc = sum(d.pt() for d in daus if d.charge() != 0 and d.pt() > 2 and abs(d.eta()) < 2.5)
        #phall = sum(d.pt() for d in daus if d.pdgId() == 22)
        #phpt1 = sum(d.pt() for d in daus if d.pdgId() == 22 and d.pt() > 1)
        #nhall = sum(d.pt() for d in daus if d.charge() ==0 and d.pdgId() != 22)
        #nhpt1 = sum(d.pt() for d in daus if d.charge() ==0 and d.pdgId() != 22)
    allGenJ.sort(key = lambda j : - j.pt())
    

    print "Gen jets in the event:"
    for j in allGenJ:
        if abs(j.eta()) > 2.5: continue
        if j.mcId != 1 or j.pt() < 20: continue
        print "------------------------"     
        print "   pt %7.2f eta %+5.2f phi %+5.2f   id %d" % (j.pt(), j.eta(), j.phi(), j.mcId) 
        daus = [ j.daughter(i) for i in xrange(j.numberOfDaughters()) ]
        for d in daus: d.dr = deltaR(d,j)
        daus.sort(key = lambda p : p.dr)
        ptsum = 0
        for d in daus:
            ptsum += d.pt()
            print "        dau pt %7.2f eta %+5.2f phi %+5.2f dr %.2f   id % +5d  ptsum %7.2f" % (d.pt(), d.eta(), d.phi(), d.dr, d.pdgId(), ptsum) 
        print "   total pt from charged particles in acceptance: %7.2f" % sum(d.pt() for d in daus if d.charge() != 0 and d.pt() > 2 and abs(d.eta()) < 2.5)
        print ""
        
        for a in args[1:]:
            h = ihpf if a.startswith("l1t") else hpf
            event.getByLabel(a, h)
            print "   --> ",a
            matches = [ p for p in h.product() if deltaR(p,j) < 0.5 ]
            for d in matches: d.dr = deltaR(d,j)
            matches.sort(key = lambda p : p.dr)
            ptsum = 0
            for d in matches:
                if d.dr < 0.4: ptsum += d.pt()
                print "        cand  pt %7.2f eta %+5.2f phi %+5.2f dr %.2f   id % +5d  ptsum %7.2f" % (d.pt(), d.eta(), d.phi(), d.dr, d.pdgId(), ptsum),
                #if "Track" in a or "TK" in a:
                #    tomatch = [g for g in daus if g.charge() != 0 and g.pt() > 1.5 and abs(d.eta()) < 2.7)
                #    matched, dr2 = bestMatch
                print ""
            print "   total pt within 0.4: %7.2f" %ptsum
            print ""
    print "\n========================\n"     

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
        print "   total pt from neutral particles in acceptance: %7.2f" % sum(d.pt() for d in daus if d.charge() == 0 and d.pt() > 1)
        print "   total pt from any     particles in acceptance: %7.2f" % sum(d.pt() for d in daus if d.pt() > (2 if d.charge() != 0 else 1))
        print ""
        
        c_tomatch = [g for g in daus if g.charge() != 0 and g.pt() > 1.5 and abs(d.eta()) < 2.7]
        n_tomatch = [g for g in daus if g.charge() == 0 and g.pt() > 1]
        a_tomatch = c_tomatch + n_tomatch

        for a in args[1:]:
            h = ihpf if a.startswith("l1t") else hpf
            event.getByLabel(a, h)
            print "   --> ",a
            matches = [ p for p in h.product() if deltaR(p,j) < 0.5 ]
            for d in matches: d.dr = deltaR(d,j)
            matches.sort(key = lambda p : p.dr)
            ptsum = 0
            # OK, try some MC matching for them
            if "PF" in a: 
                mcmatch =  dict((d,[]) for d in matches )
                # first, match and remove the charged
                charged = [d for d in matches if d.charge() != 0]
                neutral = [d for d in matches if d.charge() == 0]
                ch_match = matchObjectCollection3( charged, c_tomatch, 0.07, filter = lambda d,g: abs(d.pt()-g.pt())/(d.pt()+g.pt()) < 0.2 )
                gen_notused = a_tomatch[:]
                reco_notused = matches[:]
                for c in charged:
                    g = ch_match[c]
                    if g == None: continue
                    mcmatch[c] = [g]
                    gen_notused.remove(g)
                    reco_notused.remove(c)
                # then assign each gen to the nearest unmatched reco that has reco pt > 0.5 * gen pt
                for g in gen_notused[:]:
                    d, dr2 = bestMatch(g, reco_notused)
                    if dr2 < 0.01 and d.pt() > 0.5*g.pt(): 
                        gen_notused.remove(g)
                        mcmatch[d].append(g)
                # then take any remaining gen and just attach it to the nearest reco, except well-matched tracks
                rematch = neutral + [ c for c in charged if c in reco_notused ]
                for g in gen_notused:
                    d, dr2 = bestMatch(g, rematch)
                    if dr2 < 0.01: mcmatch[d].append(g) 
            for d in matches:
                if d.dr < 0.4: ptsum += d.pt()
                print "        cand  pt %7.2f eta %+5.2f phi %+5.2f dr %.2f   id % +5d  ptsum %7.2f" % (d.pt(), d.eta(), d.phi(), d.dr, d.pdgId(), ptsum),
                if "Track" in a or "TK" in a:
                    g, dr2 = bestMatch(d, c_tomatch)
                    if (dr2 < .01): print " --> match with  pt %7.2f eta %+5.2f phi %+5.2f dr %.3f   id % +5d " % (g.pt(), g.eta(), g.phi(), sqrt(dr2), g.pdgId()),
                    else:           print " --> unmatched",
                elif "PF" in a:
                        incone = mcmatch[d]
                        if len(incone) == 0:
                            print " --> unmatched",
                        elif len(incone) == 1:
                            g = incone[0]
                            print " --> match with  pt %7.2f eta %+5.2f phi %+5.2f dr %.3f   id % +5d " % (g.pt(), g.eta(), g.phi(), deltaR(g,d), g.pdgId()),
                        else:
                            print " --> match with  pt %7.2f from %2d particles: " % (sum(g.pt() for g in incone),len(incone)),
                            for g in incone:
                                print "%+d[pt %.1f, dr %.2f]" % (g.pdgId(), g.pt(), deltaR(g,d)),
                print ""
            print "   total pt within 0.4: %7.2f" %ptsum
            if "PF" in a:
                print "   total charged pt within 0.4: %7.2f" % sum(d.pt() for d in charged if d.dr < 0.4)
                print "   total neutral pt within 0.4: %7.2f" % sum(d.pt() for d in neutral if d.dr < 0.4)

            print ""
    print "\n========================\n"     

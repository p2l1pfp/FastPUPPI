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
parser.add_option("-g", "--global", dest="ged", action="store_true", default=False, help="Look at the global event description instead of individual PFJets")
parser.add_option("--min-pt", dest="minPt", type=float, default=20, help="Number of events to consider")
parser.add_option("--max-eta", dest="maxEta", type=float, default=5, help="Number of events to consider")
parser.add_option("--ld", "--link-debug", dest="linkDebug", default=None, help="show linking information")
parser.add_option("-J", "--jet", dest="jets", type=str, action="append", default=[], help="Specific jet to select eta,phi")
options, args = parser.parse_args()

events = Events(args[0])

handles_ = {}; loadedCaloTools_ = False
def getHandle(coll):
    global handles_, loadedCaloTools_
    if coll not in handles_: 
        handles_[coll] = Handle(coll)
        if coll in ("BXVector<l1t::CaloTower>","BXVector<l1t::CaloCluster>") and not loadedCaloTools_:
            ROOT.gSystem.Load("libL1TriggerL1TCalorimeter.so");
            ROOT.gInterpreter.ProcessLine("#include <L1Trigger/L1TCalorimeter/interface/CaloTools.h>")
            loadedCaloTools_ = True
    return handles_[coll]
hopf = getHandle("std::vector<reco::PFCandidate>")
hpf1 = getHandle("std::vector<l1t::PFCandidate>")
hpfc = getHandle("std::vector<l1t::PFCluster>")
hpft = getHandle("std::vector<l1t::PFTrack>")
htmu = getHandle("BXVector<l1t::Muon>")
hop  = getHandle("edm::OwnVector<reco::Candidate>")
hhgc = getHandle("BXVector<l1t::HGCalMulticluster>")
hhgt = getHandle("BXVector<l1t::HGCalTower>")

genj  = Handle("std::vector<reco::GenJet>")
genp  = Handle("std::vector<reco::GenParticle>")
geno  = Handle("math::XYZPointF")
recf  = Handle("float")

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
        if abs(j.eta()) > options.maxEta: continue
        #daus = [ j.daughter(i) for i in xrange(j.numberOfDaughters()) ]
        #chall = sum(d.pt() for d in daus if d.charge() != 0)
        #chacc = sum(d.pt() for d in daus if d.charge() != 0 and d.pt() > 2 and abs(d.eta()) < 2.5)
        #phall = sum(d.pt() for d in daus if d.pdgId() == 22)
        #phpt1 = sum(d.pt() for d in daus if d.pdgId() == 22 and d.pt() > 1)
        #nhall = sum(d.pt() for d in daus if d.charge() ==0 and d.pdgId() != 22)
        #nhpt1 = sum(d.pt() for d in daus if d.charge() ==0 and d.pdgId() != 22)
    allGenJ.sort(key = lambda j : - j.pt())
    event.getByLabel("genParticles:xyz0", geno)
    gen_Z0 = geno.product().Z()
   
    rec_Z0 = 0;
    event.getByLabel("InfoOut:z0", recf)
    if recf.isValid():
        rec_Z0 = recf.product()[0];
        print "PV  generated %+7.3f  reconstructed %+7.3f  diff %+5.3f " % (gen_Z0, rec_Z0, rec_Z0-gen_Z0)

    if options.ged:
        allGenJ = [None]
        print "Full event:"
    elif options.jets:
        centers = [ map(float,c.split(",")) for c in options.jets ]
        print "Selected gen jets in the event:"
        allGenJ = [ j for j in allGenJ if min(deltaR(j.eta(), j.phi(), c[0], c[1]) for c in centers) < 0.2 ]
    else:
        print "Gen jets in the event:"
    for j in allGenJ:
        if j:
            if abs(j.eta()) > 2.5: continue
            #if j.mcId != 1 or j.pt() < 20: continue
            if j.pt() < options.minPt: continue
            print "------------------------"     
            print "   pt %7.2f eta %+5.2f phi %+5.2f   id %d" % (j.pt(), j.eta(), j.phi(), j.mcId) 
            daus = [ j.daughter(i) for i in xrange(j.numberOfDaughters()) ]
            for d in daus: d.dr = deltaR(d,j)
            daus.sort(key = lambda p : p.dr)
        else:
            daus = [ p for p in genp.product() if p.status() == 1 and abs(p.pdgId()) not in (12,14,16) and p.pt() > 0.5 and abs(p.eta()) < options.maxEta ]
            daus.sort(key = lambda p : -p.pt())
            for d in daus: d.dr = 0
        ptsum = 0
        for d in daus:
            ptsum += d.pt()
            print "        dau pt %7.2f eta %+5.2f phi %+5.2f dr %.2f   id % +5d charge %+1d  ptsum %7.2f" % (d.pt(), d.eta(), d.phi(), d.dr, d.pdgId(), d.charge(), ptsum) 
        print "   total pt from charged particles in acceptance: %7.2f" % sum(d.pt() for d in daus if d.charge() != 0 and d.pt() > 2 and abs(d.eta()) < 2.5)
        print "   total pt from neutral particles in acceptance: %7.2f" % sum(d.pt() for d in daus if d.charge() == 0 and d.pt() > 1)
        print "   total pt from photons           in acceptance: %7.2f" % sum(d.pt() for d in daus if d.charge() == 0 and d.pdgId() == 22 and d.pt() > 1)
        print "   total pt from other neutrals    in acceptance: %7.2f" % sum(d.pt() for d in daus if d.charge() == 0 and d.pdgId() != 22 and d.pt() > 1)
        print "   total pt from any     particles in acceptance: %7.2f" % sum(d.pt() for d in daus if d.pt() > (2 if d.charge() != 0 else 1))
        if not j:
            print "   total pt from charged particles (pt > 0.5)   : %7.2f" % sum(d.pt() for d in daus if d.charge() != 0)
            print "   total pt from neutral particles (pt > 0.5)   : %7.2f" % sum(d.pt() for d in daus if d.charge() == 0 )
            print "   total pt from photons           (pt > 0.5)   : %7.2f" % sum(d.pt() for d in daus if d.charge() == 0 and d.pdgId() == 22)
            print "   total pt from other neutrals    (pt > 0.5)   : %7.2f" % sum(d.pt() for d in daus if d.charge() == 0 and d.pdgId() != 22)
            print "   total pt from any     particles (pt > 0.5)   : %7.2f" % sum(d.pt() for d in daus)
        print ""
        
        c_tomatch = [g for g in daus if (g.charge() != 0 and g.pt() > 1.5 and abs(g.eta()) < 2.7)]
        n_tomatch = [g for g in daus if g.charge() == 0 and g.pt() > 1]
        a_tomatch = c_tomatch + n_tomatch

        layers = {}
        for ia,a in enumerate(args[1:]):
            layerlabel = chr(ia+ord("A"))
            layers[layerlabel] = []
            alabel = a
            if "@" in a:
                (atype,alabel) = a.split("@")
                h = getHandle(atype)
            elif a.startswith("l1tPF"): h = ihpf
            elif a.startswith("pfClu"): h = hpfc
            elif a.startswith("pfTra"): h = hpft
            elif a.startswith("l1pf") : h = hpf1
            elif a.startswith("hgcalBackEndLayer2Producer"): 
                h = hhgc; alabel += ":HGCalBackendLayer2Processor3DClustering"
            elif a.startswith("hgcalTowerProducer"): 
                h = hhgt; alabel += ":HGCalTowerProcessor"
            elif a.startswith("simGmt"): h = htmu
            else: h = hop
            event.getByLabel(alabel, h)
            objs = h.product()
            print "    %s   --> %s (%d)" %(layerlabel,a,objs.size())
            if ("l1t::CaloTower" in a) or ("l1t::CaloCluster" in a):
                objs = [ o for o in objs if o.hwPt() > 0 ]
                p4unpack = ROOT.l1t.CaloTools.p4MP if ":MP" in a else ROOT.l1t.CaloTools.p4Demux
                for o in objs: 
                    o.setP4(p4unpack(o))
            if j:
                matches = [ p for p in objs if deltaR(p,j) < 0.5 ]
                for d in matches: d.dr = deltaR(d,j)
                matches.sort(key = lambda p : p.dr)
            else:
                matches = [ p for p in objs if abs(p.eta()) < options.maxEta ]
                for d in matches: d.dr = 0
                matches.sort(key = lambda p : -p.pt())

            ptsum = 0
            # OK, try some MC matching for them
            if "PF" in a or "Puppi" in a: 
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
            elif "Calo" in a or "Cluster" in a: 
                mcmatch =  dict((d,[]) for d in matches )
                for g in a_tomatch:
                    d, dr2 = bestMatch(g, matches)
                    if dr2 < 0.02: mcmatch[d].append(g)
            for im,d in enumerate(matches):
                if d.dr < 0.4: ptsum += d.pt()
                print "  %s%02d   cand  pt %7.2f eta %+5.2f phi %+5.2f dr %.2f   id % +5d  ptsum %7.2f" % (layerlabel,im, d.pt(), d.eta(), d.phi(), d.dr, d.pdgId(), ptsum),
                layers[layerlabel].append(("%s%02d"%(layerlabel,im), d, (d.caloEta(), d.caloPhi()) if "pfTrack" in a else (d.eta(), d.phi())))
                if a.startswith("pfCl"):
                    print "  emEt %7.2f isEM %1d  " % (d.emEt(), d.isEM()),
                if "Track" in a or "TK" in a:
                    if "pfTrack" in a:
                        print "  calo eta %+5.2f phi %+5.2f " % (d.caloEta(), d.caloPhi()),
                    print "  charge %+1d  dz %+7.2f "% (d.charge(), d.vertex().Z()-rec_Z0),
                    g, dr2 = bestMatch(d, c_tomatch)
                    if (dr2 < .01): 
                        print " --> match with  pt %7.2f eta %+5.2f phi %+5.2f dr %.3f   id % +5d " % (g.pt(), g.eta(), g.phi(), sqrt(dr2), g.pdgId()),
                        # check match back
                        db, dr2b = bestMatch(g, matches)
                        if db != d: print " <-- NOT matched back (gen prefers reco pt %7.2f eta %+5.2f phi %+5.2f dr %.3f   id % +5d)"  % (db.pt(), db.eta(), db.phi(), sqrt(dr2b), db.pdgId()),
                    else:           
                        print " --> unmatched",
                elif "PF" in a or "Puppi" in a or "Calo" in a:
                        if "Calo" not in a:
                            print "  charge %+1d"% (d.charge()),
                            if d.charge():
                                print "  dz %+7.2f "% (d.vertex().Z()-rec_Z0),
                            else:
                                print "     "+" "*7+" ",
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
                elif "simGmt" in a:
                    print "  charge %s  quality %d "% (("%+d" % d.hwCharge()) if d.hwChargeValid() else "n/a", d.hwQual()),
                    g, dr2 = bestMatch(d, [c for c in c_tomatch if abs(c.pdgId()) == 13])
                    if (dr2 < .1): 
                        print " --> match with  pt %7.2f eta %+5.2f phi %+5.2f dr %.3f   id % +5d " % (g.pt(), g.eta(), g.phi(), sqrt(dr2), g.pdgId()),
                    else:           
                        print " --> unmatched",
                print ""
            if len(matches) > 1: print "   total pt within 0.4: %7.2f" %ptsum
            if "PF" in a or "Puppi" in a:
                print "   total charged pt within 0.4: %7.2f" % sum(d.pt() for d in charged if d.dr < 0.4)
                print "   total neutral pt within 0.4: %7.2f" % sum(d.pt() for d in neutral if d.dr < 0.4)
                print "   total photon  pt within 0.4: %7.2f" % sum(d.pt() for d in neutral if d.dr < 0.4 and d.pdgId() == 22)
                print "   total neu had pt within 0.4: %7.2f" % sum(d.pt() for d in neutral if d.dr < 0.4 and d.pdgId() != 22)

            print ""
        if options.linkDebug:
            for link in options.linkDebug.split(","):
                layernames, drmax = link, 0.2
                if ":" in link:
                    layernames = link.split(":")[0]
                    drmax = float(link.split(":")[1])
                for (ifrom, ofrom, (etafrom,phifrom)) in layers[layernames[0]]:
                    print "   %s pt %6.2f -> " % (ifrom, ofrom.pt()),
                    ptsum = 0; matches = []
                    for (ito, oto, (etato,phito)) in layers[layernames[1]]:
                        dr = deltaR(etafrom,phifrom,etato,phito)
                        if dr < drmax: 
                            matches.append( (ito,oto.pt(),dr,oto.pt()-ofrom.pt()) );
                            ptsum += oto.pt()
                    matches.sort(key = lambda (i,pt,dr,dpt) : dr)
                    for m in matches:
                        print " %s pt %6.2f (dr %.3f, dpt %5.2f)   " % m,
                    if len(matches) > 1:
                        print " SUM pt %6.2f (dpt %5.2f)   " % (ptsum, ptsum-ofrom.pt()),
                    print ""
                print ""
            print ""
    print "\n========================\n"     

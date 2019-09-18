#!/usr/bin/env python
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(True)

import os
from sys import argv

from collections import defaultdict

from PhysicsTools.HeppyCore.utils.deltar import *
from PhysicsTools.HeppyCore.statistics.tree import *


from optparse import OptionParser
parser = OptionParser("%(prog) -w what infile1 infile2 ")
parser.add_option("--id", "--mc_id", dest="mc_id", type=int, nargs=1, default=0, help="MC ")
parser.add_option("-w", "--what", dest="what", action="append", default=[])
parser.add_option("-o", "--out", dest="out", default="correlatedRespTuple.root")
options, args = parser.parse_args()
nfiles = len(args)

fout = ROOT.TFile(options.out, "RECREATE")
otree = Tree("t","t","F")
for x in "run","lumi","event","mc_id": otree.var(x,int);
for x in "pt","eta","phi": otree.var("mc_"+x);
for j in range(nfiles):
    for x in options.what:  otree.var("%s_%d" % (x,j))

events = {}

for ia, arg in enumerate(args):
    tfile = ROOT.TFile(arg)
    tree = tfile.Get("ntuple/tree")
    for e in tree:
        if e.mc_id != options.mc_id: continue
        evid = e.run, e.lumi, e.event
        if evid not in events: 
            if ia != 0: continue
            events[evid] = {}
        pid = tuple(map(int, (e.mc_pt*100, e.mc_eta*100, e.mc_phi*100)))
        if ia == 0: 
            events[evid][pid] = [ [ getattr(e,x) for x in options.what ] ]
        else:
            if pid not in events[evid]: continue
            events[evid][pid].append( [ getattr(e,x) for x in options.what ] )

print "Found %d events" % len(events)

otree.fill("mc_id", options.mc_id)
for evid, pids in sorted(events.iteritems()):
    otree.fill("run", evid[0])
    otree.fill("lumi", evid[1])
    otree.fill("event", evid[2])
    for pid, data in pids.iteritems():
        otree.fill("mc_pt", 0.01*pid[0])
        otree.fill("mc_eta", 0.01*pid[1])
        otree.fill("mc_phi", 0.01*pid[2])
        for j, row in enumerate(data):
            for (x,v) in zip(options.what,row):
                otree.fill("%s_%d" % (x,j), v)
        otree.tree.Fill()
print ""
print "Output entries: %d" % otree.tree.GetEntries()
fout.cd()
fout.WriteTObject(otree.tree)
fout.Close()
        



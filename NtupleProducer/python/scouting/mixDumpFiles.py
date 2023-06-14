#!/usr/bin/env python3
from optparse import OptionParser
parser = OptionParser("%(prog) infile [ src [ dst ] ]")
parser.add_option("--mux-bx", dest="bxMux", type=int, default=1)
parser.add_option("--mux-orbits", dest="orbitMux", type=int, default=1)
parser.add_option("--type", dest="dumpType", type=str, default="puppi")
parser.add_option("-n","--norbits", dest="orbits", type=int, default=100)
parser.add_option("--bx", dest="nbx", type=int, default=3564, help="BX per orbit")
parser.add_option("-o","--out", dest="out", type=str, default="dump")
parser.add_option("--orbit-header","--OH", dest="orbitHeader", default=False, action="store_true")
options, args = parser.parse_args()

import sys, random
rnd = random.Random(37)

class PuppiFile:
    def __init__(self, fname):
        self._fname = fname
        if fname != None:
            fin = open(fname,'rb')
            index = []
            sumn = 0; nent = 0
            while fin.readable():
                w64 = fin.read(8)
                if not w64: break
                npuppi = (int.from_bytes(w64, sys.byteorder) & 0xFFF)
                puppidata = fin.read(8*npuppi)
                index.append((npuppi,puppidata))
                sumn += npuppi
                nent += 1
                if nent % 10000 == 0:
                    print("Processed %7d entries" % nent)
            self._index = index
            print("File %s with %d events, avg multiplicity %.2f" % (fname,len(index),sumn/float(len(index))));
    def push(self,lout,rnd,header):
        if self._fname != None:
            (n,data) = rnd.choice(self._index)
            if (n >= 256):
                print("Warning, cropping event with %d candidates" % n)
                n = 255
                data = data[:8*255]
            header = header | (n & 0xFFF)
            lout.append(header.to_bytes(8, sys.byteorder))
            lout.append(data)
        else:
            lout.append(header.to_bytes(8, sys.byteorder))
    def write(self,fout,rnd,header):
        if self._fname != None:
            (n,data) = rnd.choice(self._index)
            if (n >= 256):
                print("Warning, cropping event with %d candidates" % n)
                n = 255
                data = data[:8*255]
            header = header | (n & 0xFFF)
            fout.write(header.to_bytes(8, sys.byteorder))
            fout.write(data)
        else:
            fout.write(header.to_bytes(8, sys.byteorder))


makers = dict(puppi=PuppiFile)

samples = []; tot = 0.
for a in args:
    fname, prob = a.split(":")
    tot += float(prob)
    if tot > 1: raise RuntimeError()
    samples.append((makers[options.dumpType](fname), tot))
if tot < 1:
    samples.append((makers[options.dumpType](None), 1))


outs = []
fname = options.out
for iorb in range(options.orbitMux):
    orbitname = ".orbit_%d_of_%d" % (iorb,options.orbitMux) if options.orbitMux > 1 else ""
    for ibx in range(options.bxMux):
        bxname = ".bx_%d_of_%d" % (ibx,options.bxMux) if options.bxMux > 1 else ""
        fname = options.out + orbitname + bxname + ".dump"
        print("Opening "+fname)
        outs.append(open(fname, "wb"))

run = 0
for orbit in range(options.orbits):
    if options.orbitHeader:
        orbitLists = [[] for o in range(options.bxMux)]
        for bx in range(options.nbx):
            ibx = bx % options.bxMux
            header = ((((run << 30) | orbit) << 12) | bx) << 12 
            u = rnd.random()
            for (f,t) in samples:
                if u < t:
                    f.push(orbitLists[ibx],rnd,header)
                    break
        for ibx,data in enumerate(orbitLists):
            imux = (orbit % options.orbitMux) * options.bxMux + ibx
            #print("Orbit %d ibx %d of length %d" % (orbit, ibx, sum(len(d)/8 for d in data)))
            orbitHeader = ((run << 30) | orbit) << 24
            for d in data: 
                orbitHeader += len(d)>>3
            outs[imux].write(orbitHeader.to_bytes(8, sys.byteorder))
            for d in data: 
                outs[imux].write(d)
    else:
        for bx in range(options.nbx):
            out = outs[(orbit % options.orbitMux) * options.bxMux + (bx % options.bxMux)]
            header = ((((run << 30) | orbit) << 12) | bx) << 12 
            u = rnd.random()
            for (f,t) in samples:
                if u < t:
                    f.write(out,rnd,header)
                    break
    if orbit % 10 == 9: 
        print("Processed %7d orbits" % (orbit+1))

for o in outs:
    o.close()

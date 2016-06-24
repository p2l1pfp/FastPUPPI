import ROOT
import copy
import sys
from math import *
from array import array
from operator import itemgetter, attrgetter

ROOT.gSystem.Load("libFWCoreFWLite.so")
ROOT.AutoLibraryLoader.enable()

from DataFormats.FWLite import Events, Handle

ROOT.gSystem.Load("libDataFormatsFWLite.so");
ROOT.gROOT.ProcessLine("namespace edm {typedef edm::Wrapper<vector<float> > Wrapper<vector<float,allocator<float> > >; }");
ROOT.gROOT.ProcessLine("namespace edm {typedef edm::Wrapper<vector<double> > Wrapper<vector<double,allocator<double> > >; }");
# ROOT.gROOT.ProcessLine("namespace edm {typedef edm::Wrapper< vector<TTTrack<edm::Ref<edm::DetSetVector<PixelDigi>,PixelDigi,edm::refhelper::FindForDetSetVector<PixelDigi> > > >; }");

############################################################
############################################
#            Job steering                  #
############################################
from optparse import OptionParser

parser = OptionParser()

parser.add_option('-b', action='store_true', dest='noX', default=False,
                  help='no X11 windows')

#parser.add_option('--makeCards', action='store_true', dest='makeCards', default=False,
#                  help='no X11 windows')
#
#parser.add_option('--channel',action="store",type="string",dest="channel",default="mu")
parser.add_option('--nPU',action="store",type="int",dest="nPU",default=50)
parser.add_option('--nBx',action="store",type="int",dest="nBx",default=50)
parser.add_option('--nEv',action="store",type="int",dest="nEv",default=500)
parser.add_option('--index',action="store",type="int",dest="index",default=1)



(options, args) = parser.parse_args()
############################################################
############################################################


def main():

    nEvents = options.nEv;
    nPU = options.nPU; 
    nBx = options.nBx;
    index = options.index;
    ofilename = "hists/fileOut_PU"+str(nPU)+"bx"+str(nBx)+"_"+str(options.index)+".root";

    dir = "/afs/cern.ch/work/n/ntran/private/Correlator/eos/cms/store/mc/TTI2023Upg14D/Zmumu_TuneZ2star_14TeV_Eta4-pythia6/GEN-SIM-DIGI-RAW/PU140bx25_PH2_1K_FB_V3-v2/00000/";
    files = []
    dirEos = "file:%s" % dir;
    files.append(dirEos+"08B6D2DB-D5E6-E311-9356-002618943975.root");
        
    events = Events( files );

    # loop over events
    count = 0
    ntotal = events.size()
    print "Nevents = "+str(ntotal)
    print "Start looping"

    # https://github.com/cms-sw/cmssw/blob/CMSSW_6_2_0_SLHC12_patch/DataFormats/L1TrackTrigger/interface/TTTrack.h
    l1trackHandle = Handle("vector<TTTrack<edm::Ref<edm::DetSetVector<PixelDigi>,PixelDigi,edm::refhelper::FindForDetSetVector<PixelDigi> > > >");
    l1trackLabel  = ("TTTracksFromPixelDigis","Level1TTTracks","RECO")

    # https://github.com/cms-sw/cmssw/blob/CMSSW_6_2_0_SLHC12_patch/DataFormats/EcalDigi/interface/EcalTriggerPrimitiveDigi.h
    ECalPrimHandle = Handle("edm::SortedCollection<EcalTriggerPrimitiveDigi,edm::StrictWeakOrdering<EcalTriggerPrimitiveDigi> >");
    ECalPrimLabel = ("simEcalTriggerPrimitiveDigis","","RECO" );

    # https://github.com/cms-sw/cmssw/blob/CMSSW_6_2_0_SLHC12_patch/DataFormats/HcalDigi/interface/HcalTriggerPrimitiveDigi.h
    HCalPrimHandle = Handle("edm::SortedCollection<HcalTriggerPrimitiveDigi,edm::StrictWeakOrdering<HcalTriggerPrimitiveDigi> >");
    HCalPrimLabel = ("simHcalTriggerPrimitiveDigis","","RECO" );

    #for event in events:
    # ctr = 0;
    for i,ev in enumerate(events):

        if i == 10: break;
        if i % 1 == 0: print "i = ", i;

        ev.getByLabel(l1trackLabel,l1trackHandle);
        tracks = l1trackHandle.product();
        print len(tracks);
        for t in tracks:
            print "t: ", t.getRInv(), t.getChi2();

        ev.getByLabel(ECalPrimLabel,ECalPrimHandle);
        ecaldigis = ECalPrimHandle.product();
        print len(ecaldigis);
        for e in ecaldigis:
            print "e: ", e.id().ieta(),e.id().iphi(),e.compressedEt();


        ev.getByLabel(HCalPrimLabel,HCalPrimHandle);
        hcaldigis = HCalPrimHandle.product();
        print len(hcaldigis);
        for h in hcaldigis:
            print "h: ", h.id().ieta(),h.id().iphi(),h.SOI_compressedEt();
            hscale = L1CaloHcalScale(0.5);
            # et = L1CaloHcalScale

#--- --- --- --- --- --- --- --- --- ---
if __name__ == '__main__':
    main();
#--- --- --- --- --- --- --- --- --- ---




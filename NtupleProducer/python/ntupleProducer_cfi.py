import FWCore.ParameterSet.Config as cms

InfoOut = cms.EDProducer('NtupleProducer',
         #L1TrackTag  = cms.InputTag('l1tPFTkProducersFromOfflineTracksStrips'),
         L1TrackTag  = cms.InputTag('l1tPFTkProducersFromL1Tracks'),
         CaloClusterTags  = cms.VInputTag( cms.InputTag('CaloInfoOut','calibrated') ),
         EmClusterTags  = cms.VInputTag( ),
         correctCaloEnergies = cms.bool(False),
         MuonTPTag   = cms.InputTag('l1tPFMuProducerFromL1Mu'), 
         genParTag   = cms.InputTag('genParticles'),
         corrector   = cms.string("FastPUPPI/NtupleProducer/data/pion_eta_phi.root"),
         ecorrector  = cms.string("FastPUPPI/NtupleProducer/data/ecorr.root"),
         trackres    = cms.string("FastPUPPI/NtupleProducer/data/tkres.root"),
         eleres      = cms.string("FastPUPPI/NtupleProducer/data/eres.root"),
         pionres     = cms.string("FastPUPPI/NtupleProducer/data/pionres.root"),
         trkPtCut    = cms.double(2.0),
         metRate     = cms.bool(False),
         etaCharged  = cms.double(2.5),
         puppiPtCut  = cms.double(4.0),
         puppiDr     = cms.double(0.5),
         vtxRes      = cms.double(0.333),
         linking = cms.PSet(
                        trackCaloDR = cms.double(0.15),
                        trackCaloNSigmaLow = cms.double(2.0),
                        trackCaloNSigmaHigh = cms.double(2.0),
                        useTrackCaloSigma = cms.bool(True), # take the uncertainty on the calo cluster from the track, for linking purposes
                        rescaleUnmatchedTrack = cms.bool(False),
                        maxInvisiblePt = cms.double(10.0),
                        ),
         outputName  = cms.untracked.string("ntuple.root"),
         debug       = cms.untracked.int32(0),
)
# these below are temporary resolutions, to be improved later
if True:
    InfoOut.simpleResolHad = cms.PSet(
            etaBins = cms.vdouble( 1.300,  1.700,  2.800,  3.200,  4.000,  5.000),
            offset  = cms.vdouble( 3.522,  0.078,  2.071,  1.708,  1.148, -0.265),
            scale   = cms.vdouble( 0.124,  0.494,  0.183,  0.257,  0.162,  0.428),
            kind    = cms.string('calo'),
            )
    InfoOut.simpleResolEm = cms.PSet(
            etaBins = cms.vdouble( 1.300,  1.700,  2.800,  3.200,  4.000,  5.000),
            offset  = cms.vdouble( 0.849,  0.626,  0.157, -1.305,  0.607, -3.985),
            scale   = cms.vdouble( 0.016,  0.097,  0.043,  0.305,  0.142,  0.626),
            kind    = cms.string('calo'),
            )
    InfoOut.simpleResolTrk  = cms.PSet(
            etaBins = cms.vdouble( 0.800,  1.200,  1.500,  2.000,  2.500),
            offset  = cms.vdouble( 0.006,  0.010,  0.010,  0.019,  0.027),
            scale   = cms.vdouble( 0.303,  0.465,  1.003,  1.219,  1.518),
            kind    = cms.string('track'),
            )


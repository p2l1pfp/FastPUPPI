import FWCore.ParameterSet.Config as cms

InfoOut = cms.EDProducer('NtupleProducer',
         #L1TrackTag  = cms.InputTag('l1tPFTkProducersFromOfflineTracksStrips'),
         L1TrackTag  = cms.InputTag('l1tPFTkProducersFromL1Tracks'),
         CaloClusterTags  = cms.VInputTag( cms.InputTag('CaloInfoOut','uncalibrated') ),
         EmClusterTags  = cms.VInputTag( cms.InputTag('l1tPFHGCalProducerFrom3DTPsEM'),cms.InputTag('l1tPFEcalProducerFromL1EGCrystalClusters'), ),
         correctCaloEnergies = cms.bool(True),
         rawEmCaloPtMin = cms.double(0.5),
         rawCaloPtMin   = cms.double(1.0),
         MuonTPTag   = cms.InputTag('l1tPFMuProducerFromL1Mu'), 
         genParTag   = cms.InputTag('genParticles'),
         genOriginTag = cms.InputTag('genParticles','xyz0'),
         corrector   = cms.string("FastPUPPI/NtupleProducer/data/pion_eta_phi.root"),
         ecorrector  = cms.string("FastPUPPI/NtupleProducer/data/ecorr.root"),
         trackres    = cms.string("FastPUPPI/NtupleProducer/data/tkres.root"),
         eleres      = cms.string("FastPUPPI/NtupleProducer/data/eres.root"),
         pionres     = cms.string("FastPUPPI/NtupleProducer/data/pionres.root"),
         correctorEmfBins = cms.uint32(11),
         correctorEmfMax  = cms.double(1.0),
         trkPtCut    = cms.double(2.0),
         metRate     = cms.bool(False),
         etaCharged  = cms.double(2.5),
         puppiDr     = cms.double(0.3),
         puppiPtCut  = cms.double(4.0), # for old algo
         puppiEtaCuts       = cms.vdouble(1.5, 2.5, 3.0, 5.5), 
         puppiPtCuts        = cms.vdouble(0.0, 3.0, 6.0, 8.0),
         puppiPtCutsPhotons = cms.vdouble(0.0, 3.0, 6.0, 8.0),
         vtxRes      = cms.double(0.333),
         vtxAlgo     = cms.string("TP"),
         vtxAdaptiveCut = cms.bool(True),
         linking = cms.PSet(
                        algo = cms.string("PFAlgo3"),
                        # track -> em linking configurables
                        trackEmDR   = cms.double(0.04), # 1 Ecal crystal size is 0.02, and ~2 cm in HGCal is ~0.007
                        # em -> calo linking configurables
                        emCaloDR    = cms.double(0.10),    # 1 Hcal tower size is ~0.09
                        caloEmPtMinFrac = cms.double(0.5), # Calo object must have an EM Et at least half of that of the EM cluster to allow linking
                        # track -> calo linking configurables
                        trackCaloDR = cms.double(0.12),
                        trackCaloNSigmaLow  = cms.double(2.0),
                        trackCaloNSigmaHigh = cms.double(1.2),
                        useTrackCaloSigma = cms.bool(True), # take the uncertainty on the calo cluster from the track, for linking purposes
                        sumTkCaloErr2 = cms.bool(True), # add up track calo errors in quadrature instead of linearly
                        # how to deal with unlinked tracks
                        maxInvisiblePt = cms.double(10.0), # max allowed pt of a track with no calo energy
                        tightTrackMinStubs = cms.uint32(6),
                        tightTrackMaxChi2  = cms.double(50),
                        tightTrackMaxInvisiblePt = cms.double(20),
                        # how to deal with neutrals
                        ecalPriority  = cms.bool(True), # take first ecal energy when making neutrals
                        # other features not turned on: reliniking of neutrals to track-matched calo clusters with track excess
                        caloReLink  = cms.bool(False),
                        caloReLinkDR = cms.double(0.3),
                        caloReLinkThreshold = cms.double(0.5),
                        # other features not turned on: matching too high pt tracks to calo but rescaling track pt (not implemented in PFAlgo3)
                        rescaleUnmatchedTrack = cms.bool(False),
         ),
         outputName  = cms.untracked.string("ntuple.root"),
         debug       = cms.untracked.int32(0),
)
# these below are temporary resolutions, to be improved later
if True:
    InfoOut.simpleResolHad = cms.PSet(
            etaBins = cms.vdouble( 1.300,  1.700,  2.800,  3.200,  4.000,  5.000),
            offset  = cms.vdouble( 3.749,  3.026,  2.611,  1.638,  1.231,  1.817),
            scale   = cms.vdouble( 0.122,  0.269,  0.138,  0.246,  0.213,  0.219),
            ptMin   = cms.vdouble( 10.00,  10.00,  10.00,  10.00,  10.00,  10.00),
            ptMax   = cms.vdouble(999999, 999999, 999999, 999999, 999999, 999999),
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


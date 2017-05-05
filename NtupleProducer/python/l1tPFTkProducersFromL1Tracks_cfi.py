import FWCore.ParameterSet.Config as cms

l1tPFTkProducersFromL1Tracks = cms.EDProducer("TkProducerFromL1Tracks",
    L1TrackTag = cms.InputTag("TTTracksFromTracklet", "Level1TTTracks"),
    MuonGunVeto = cms.bool(False),
)

import FWCore.ParameterSet.Config as cms

l1tPFTkProducersFromOfflineTracksStrips = cms.EDProducer("TkProducerFromOfflineTracks",
    TrackTag = cms.InputTag("generalTracks"),
    TrackSel = cms.string("hitPattern.stripLayersWithMeasurement >= 2"),
    MuonGunVeto = cms.bool(False),
)

l1tPFTkProducersFromOfflineTracksAll = cms.EDProducer("TkProducerFromOfflineTracks",
    TrackTag = cms.InputTag("generalTracks"),
    TrackSel = cms.string(""),
    MuonGunVeto = cms.bool(False),
)

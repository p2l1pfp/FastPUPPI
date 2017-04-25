import FWCore.ParameterSet.Config as cms

l1tPFTkProducersFromOfflineTracksStrips = cms.EDProducer("TkProducerFromOfflineTracks",
    TrackTag = cms.InputTag("generalTracks"),
    TrackSel = cms.string("hitPattern.stripLayersWithMeasurement >= 2"),
)

l1tPFTkProducersFromOfflineTracksAll = cms.EDProducer("TkProducerFromOfflineTracks",
    TrackTag = cms.InputTag("generalTracks"),
    TrackSel = cms.string(""),
)

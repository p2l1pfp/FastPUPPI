import FWCore.ParameterSet.Config as cms

l1tPFHGCalEEProducerFromOfflineRechits = cms.EDProducer("HGCalProducerFromOfflineRechits",
    src = cms.InputTag("HGCalRecHit:HGCEERecHits"),
    eMin = cms.double(0.0),
    etMin = cms.double(0.0), 
)

l1tPFHGCalFHProducerFromOfflineRechits = cms.EDProducer("HGCalProducerFromOfflineRechits",
    src = cms.InputTag("HGCalRecHit:HGCHEFRecHits"),
    eMin = cms.double(0.0),
    etMin = cms.double(0.0), 
)

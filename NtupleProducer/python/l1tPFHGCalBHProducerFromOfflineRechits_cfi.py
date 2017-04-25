import FWCore.ParameterSet.Config as cms

l1tPFHGCalBHProducerFromOfflineRechits = cms.EDProducer("HGCalBHProducerFromOfflineRechits",
    src = cms.InputTag("HGCalRecHit:HGCHEBRecHits"),
    eMin = cms.double(0.0),
    etMin = cms.double(0.0), 
)

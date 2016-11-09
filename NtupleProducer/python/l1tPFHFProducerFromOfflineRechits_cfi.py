import FWCore.ParameterSet.Config as cms

l1tPFHFProducerFromOfflineRechits = cms.EDProducer("HFProducerFromOfflineRechits",
    src = cms.InputTag("hfreco"),
    eMin = cms.double(0.0),
    etMin = cms.double(0.0), # a 0.2 GeV cut appeared to be fine
)

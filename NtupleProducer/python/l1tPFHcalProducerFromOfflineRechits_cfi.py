import FWCore.ParameterSet.Config as cms

l1tPFHcalProducerFromOfflineRechits = cms.EDProducer("HcalProducerFromOfflineRechits",
    src = cms.InputTag("hbhereco"),
    subdets = cms.vuint32(1), # 1=HB
    eMin = cms.double(0.0),
    etMin = cms.double(0.0), 
)

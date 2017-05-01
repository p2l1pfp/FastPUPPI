import FWCore.ParameterSet.Config as cms

l1tPFEcalProducerFromTPDigis = cms.EDProducer("EcalProducerFromTPDigi",
    EcalTPTag = cms.InputTag("simEcalEBTriggerPrimitiveDigis"),
    etMin = cms.double(0.0),
)

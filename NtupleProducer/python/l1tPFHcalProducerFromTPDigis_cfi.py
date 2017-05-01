import FWCore.ParameterSet.Config as cms

l1tPFHcalProducerFromTPDigis = cms.EDProducer("HcalProducerFromTPDigi",
    HcalTPTag = cms.InputTag("simHcalTriggerPrimitiveDigis"),
    #subdets = cms.vuint32(1), # 1=HB
)

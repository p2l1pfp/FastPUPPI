import FWCore.ParameterSet.Config as cms

l1tPFHGCalProducerFromTriggerCells = cms.EDProducer("HGCalProducerFromTriggerCells",
    src = cms.InputTag("hgcalTriggerPrimitiveDigiProducer","calibratedTriggerCells"),
    etMin = cms.double(0.0), 
)

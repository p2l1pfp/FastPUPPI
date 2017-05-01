import FWCore.ParameterSet.Config as cms

l1tPFHGCalProducerFrom3DTPs = cms.EDProducer("HGCalProducerFrom3DTPs",
    src = cms.InputTag("hgcalTriggerPrimitiveDigiProducer","cluster3D"),
    etMin = cms.double(0.0), 
)



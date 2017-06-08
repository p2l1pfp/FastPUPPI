import FWCore.ParameterSet.Config as cms

l1tPFHGCalProducerFrom3DTPsEM = cms.EDProducer('GetEMPart',
        src = cms.InputTag('l1tPFHGCalProducerFrom3DTPs'),
)

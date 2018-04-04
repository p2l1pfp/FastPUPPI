import FWCore.ParameterSet.Config as cms

from L1Trigger.L1CaloTrigger.l1EGammaCrystalsProducer_cfi import l1EGammaCrystalsProducer as L1EGammaCrystalsProducer

l1tPFEcalProducerFromL1EGCrystalClusters = cms.EDProducer("EcalProducerFromL1EGCrystalCluster",
    src = cms.InputTag("L1EGammaCrystalsProducer","L1EGXtalClusterNoCuts"),
    etMin = cms.double(0.0), 
)



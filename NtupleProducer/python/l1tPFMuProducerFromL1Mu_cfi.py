import FWCore.ParameterSet.Config as cms

l1tPFMuProducerFromL1Mu = cms.EDProducer("MuProducerFromL1Mu",
    MuonTag   = cms.InputTag('simGmtStage2Digis'),
    MuonGunVeto = cms.bool(False),
)

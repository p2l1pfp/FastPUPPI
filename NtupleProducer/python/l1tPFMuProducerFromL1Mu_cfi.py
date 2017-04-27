import FWCore.ParameterSet.Config as cms

l1tPFMuProducerFromL1Mu = cms.EDProducer("MuProducerFromL1Mu",
    MuonTag   = cms.InputTag('gmtStage2Digis','Muon'), 
    MuonGunVeto = cms.bool(False),
)

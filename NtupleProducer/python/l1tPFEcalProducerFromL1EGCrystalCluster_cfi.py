import FWCore.ParameterSet.Config as cms

L1EGammaCrystalsProducer = cms.EDProducer("L1EGCrystalClusterProducer",
   EtminForStore = cms.double(0.),
   EcalTpEtMin = cms.untracked.double(0.5), # 500 MeV default per each Ecal TP
   EtMinForSeedHit = cms.untracked.double(1.0), # 1 GeV decault for seed hit
   debug = cms.untracked.bool(False),
   useRecHits = cms.bool(False),
   doBremClustering = cms.untracked.bool(True), # Should always be True when using for E/Gamma
   ecalTPEB = cms.InputTag("simEcalEBTriggerPrimitiveDigis","","HLT"),
   ecalRecHitEB = cms.InputTag("ecalRecHit","EcalRecHitsEB","RECO"),
   hcalRecHit = cms.InputTag("hbhereco"),
   hcalTP = cms.InputTag("simHcalTriggerPrimitiveDigis","","HLT"),
   useTowerMap = cms.untracked.bool(False)
)

l1tPFEcalProducerFromL1EGCrystalClusters = cms.EDProducer("EcalProducerFromL1EGCrystalCluster",
    src = cms.InputTag("L1EGammaCrystalsProducer","L1EGXtalClusterNoCuts"),
    etMin = cms.double(0.0), 
)



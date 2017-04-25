import FWCore.ParameterSet.Config as cms

l1tPFEcalProducerFromOfflineRechits = cms.EDProducer("EcalProducerFromOfflineRechits",
    src = cms.InputTag("ecalRecHit:EcalRecHitsEB"),
    eMin = cms.double(0.0),
    etMin = cms.double(0.0),
    ## Note: Wisc people say that a 0.5 GeV threshold is good against PU
)

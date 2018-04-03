import FWCore.ParameterSet.Config as cms

CaloInfoOut = cms.EDProducer('CaloNtupleProducer',
         #EcalTPTags  = cms.VInputTag( cms.InputTag('l1tPFEcalProducerFromOfflineRechits','towers'),
         #                             cms.InputTag('l1tPFHGCalEEProducerFromOfflineRechits','towers') ),
         #HcalTPTags  = cms.VInputTag( cms.InputTag('l1tPFHcalProducerFromOfflineRechits','towers'),
         #                             cms.InputTag('l1tPFHGCalFHProducerFromOfflineRechits','towers'),
         #                             cms.InputTag('l1tPFHGCalBHProducerFromOfflineRechits','towers'),
         #                             cms.InputTag('l1tPFHFProducerFromOfflineRechits','towers') ),
         #EcalTPTags = cms.VInputTag( cms.InputTag('l1tPFEcalProducerFromTPDigis','towers'),
         #                            cms.InputTag('l1tPFHGCalProducerFromTriggerCells','towersEE')),
         EcalTPTags = cms.VInputTag( cms.InputTag('l1tPFHGCalProducerFrom3DTPsEM'),
                                     cms.InputTag('l1tPFEcalProducerFromL1EGCrystalClusters')),
         HcalTPTags =  cms.VInputTag( cms.InputTag('l1tPFHcalProducerFromTPDigis'),
                                      cms.InputTag('l1tPFHGCalProducerFromTriggerCells','towersFHBH')),
                                      #cms.InputTag('l1tPFHGCalBHProducerFromOfflineRechits','towers')),
         genParTag   = cms.InputTag('genParticles'),
         zeroSuppress = cms.bool(False),
         corrector   = cms.string("FastPUPPI/NtupleProducer/data/hadcorr.root"),
         ecorrector  = cms.string("FastPUPPI/NtupleProducer/data/ecorr.root"),
         correctIndividualEms = cms.bool(False),
         correctorEmfBins = cms.uint32(8),
         correctorEmfMax  = cms.double(1.0),
         caloClusterer = cms.PSet(
            ecal = cms.PSet(
                zsEt = cms.double(0.4),
                seedEt = cms.double(0.5),
                minClusterEt = cms.double(0.5),
                energyWeightedPosition = cms.bool(True),
                energyShareAlgo = cms.string("fractions"),
            ), 
            hcal = cms.PSet(
                zsEt = cms.double(0.4),
                seedEt = cms.double(0.5),
                minClusterEt = cms.double(0.8),
                energyWeightedPosition = cms.bool(True),
                energyShareAlgo = cms.string("fractions"),
            ),
            linker = cms.PSet(
                hoeCut = cms.double(0.1),
                minPhotonEt = cms.double(1.0),
                minHadronEt = cms.double(1.0),
                useCorrectedEcal = cms.bool(True), # use corrected ecal enery in linking
            ),
         ),
         outputName  = cms.untracked.string("calontuple.root"),
         debug       = cms.untracked.int32(0),
)



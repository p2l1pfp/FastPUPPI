import FWCore.ParameterSet.Config as cms

CaloInfoOut = cms.EDProducer('CaloNtupleProducer',
         #EcalTPTags  = cms.VInputTag( cms.InputTag('l1tPFEcalProducerFromOfflineRechits','towers'),
         #                             cms.InputTag('l1tPFHGCalEEProducerFromOfflineRechits','towers') ),
         #HcalTPTags  = cms.VInputTag( cms.InputTag('l1tPFHcalProducerFromOfflineRechits','towers'),
         #                             cms.InputTag('l1tPFHGCalFHProducerFromOfflineRechits','towers'),
         #                             cms.InputTag('l1tPFHGCalBHProducerFromOfflineRechits','towers'),
         #                             cms.InputTag('l1tPFHFProducerFromOfflineRechits','towers') ),
         EcalTPTags = cms.VInputTag( cms.InputTag('l1tPFEcalProducerFromTPDigis','towers'),
                                     cms.InputTag('l1tPFHGCalProducerFromTriggerCells','towersEE')),
         HcalTPTags =  cms.VInputTag( cms.InputTag('l1tPFHcalProducerFromTPDigis'),
                                      cms.InputTag('l1tPFHGCalProducerFromTriggerCells','towersFHBH') ),
                                    # could maybe add cms.InputTag('l1tPFHGCalBHProducerFromOfflineRechits','towers')
         genParTag   = cms.InputTag('genParticles'),
         zeroSuppress = cms.bool(False),
         corrector   = cms.string("FastPUPPI/NtupleProducer/data/pion_eta_phi.root"),
         ecorrector  = cms.string("FastPUPPI/NtupleProducer/data/ecorr.root"),
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

# these below are temporary calibrations, to be improved later
if True:
    CaloInfoOut.simpleCorrEm = cms.PSet(
                etaBins = cms.vdouble( 0.500,  1.000,  1.500,  2.000,  2.500,  3.000),
                offset  = cms.vdouble(-1.402, -1.733, -2.007, -0.983, -1.140, -1.362),
                scale   = cms.vdouble( 0.977,  0.976,  0.960,  0.915,  0.949,  0.986),
                )
    CaloInfoOut.simpleCorrHad = cms.PSet(
            etaBins = cms.vdouble( 0.500,  0.500,  0.500,  1.000,  1.000,  1.000,  1.500,  1.500,  1.500,  2.000,  2.000,  2.000,  2.500,  2.500,  2.500,  3.000,  3.000,  3.000,  3.500,  4.000,  4.500,  5.000),
            emfBins = cms.vdouble( 0.125,  0.500,  0.875,  0.125,  0.500,  0.875,  0.125,  0.500,  0.875,  0.125,  0.500,  0.875,  0.125,  0.500,  0.875,  0.125,  0.500,  0.875,  1.100,  1.100,  1.100,  1.100),
            offset  = cms.vdouble(-2.765, -1.101, -3.387, -3.069, -0.775, -2.698, -5.154,  0.823, -1.693, -2.871, -0.408, -1.276, -2.205, -0.923, -2.050, -3.338,  0.284, -1.839, -3.910, -3.679, -3.361, -4.131),
            scale   = cms.vdouble( 0.951,  0.719,  0.721,  0.977,  0.702,  0.722,  0.915,  0.651,  0.647,  0.586,  0.671,  0.722,  0.608,  0.670,  0.732,  0.544,  0.578,  0.674,  1.157,  1.154,  1.060,  0.744),
            )


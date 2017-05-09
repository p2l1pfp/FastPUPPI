import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras

process = cms.Process("IN", eras.Phase2C2_timing)
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.Geometry.GeometryExtended2023D4Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

#process.load('Configuration.StandardSequences.SimL1Emulator_cff')
#process.load('Configuration.StandardSequences.L1TrackTrigger_cff')
#process.load('SimCalorimetry.HcalTrigPrimProducers.hcalTTPDigis_cff')
process.load('FastPUPPI.NtupleProducer.reprocess_L1Phase2_MC_cff')

import EventFilter.EcalRawToDigi.EcalUnpackerData_cfi
process.ecalDigis = EventFilter.EcalRawToDigi.EcalUnpackerData_cfi.ecalEBunpacker.clone()
import EventFilter.ESRawToDigi.esRawToDigi_cfi
process.ecalPreshowerDigis = EventFilter.ESRawToDigi.esRawToDigi_cfi.esRawToDigi.clone()
import EventFilter.HcalRawToDigi.HcalRawToDigi_cfi
process.hcalDigis = EventFilter.HcalRawToDigi.HcalRawToDigi_cfi.hcalDigis.clone()
process.load('RecoLuminosity.LumiProducer.bunchSpacingProducer_cfi')
process.load('RecoLocalCalo.Configuration.RecoLocalCalo_cff')

process.load('FastPUPPI.NtupleProducer.l1tPFCaloProducersFromOfflineRechits_cff')

process.load('FastPUPPI.NtupleProducer.l1tPFEcalProducerFromTPDigis_cfi')
process.load('FastPUPPI.NtupleProducer.l1tPFHcalProducerFromTPDigis_cfi')
process.load('FastPUPPI.NtupleProducer.l1tPFHGCalProducerFrom3DTPs_cfi')
process.load('FastPUPPI.NtupleProducer.l1tPFHGCalProducerFromTriggerCells_cfi')
process.load('FastPUPPI.NtupleProducer.l1tPFTkProducersFromL1Tracks_cfi')
process.load('FastPUPPI.NtupleProducer.l1tPFMuProducerFromL1Mu_cfi')
process.l1tPFMuProducerFromL1Mu.MuonTag = "simGmtStage2Digis"
process.load('FastPUPPI.NtupleProducer.l1tPFEcalProducerFromL1EGCrystalCluster_cfi')

#process.load('FastPUPPI.NtupleProducer._cfi')

process.source = cms.Source("PoolSource",
    #fileNames = cms.untracked.vstring('/store/mc/PhaseIISpring17D/SinglePhoton_FlatPt-8to150/GEN-SIM-DIGI-RAW/NoPU_90X_upgrade2023_realistic_v9-v1/70000/008FC219-132A-E711-893F-A0000420FE80.root')
    fileNames = cms.untracked.vstring('file:SinglePion_17D.root',)
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(40))

process.s = cms.Sequence(
    process.reprocess_L1Phase2_MC +
    process.ecalDigis + process.ecalPreshowerDigis + process.hcalDigis + process.bunchSpacingProducer + process.ecalLocalRecoSequence + process.hcalLocalRecoSequence + process.hcalGlobalRecoSequence + process.hgcalLocalRecoSequence +
    process.l1tPFEcalProducerFromOfflineRechits + process.l1tPFHcalProducerFromOfflineRechits + process.l1tPFHFProducerFromOfflineRechits + process.l1tPFHGCalEEProducerFromOfflineRechits + process.l1tPFHGCalFHProducerFromOfflineRechits + process.l1tPFHGCalBHProducerFromOfflineRechits +
    process.l1tPFEcalProducerFromTPDigis +
    process.L1EGammaCrystalsProducer + process.l1tPFEcalProducerFromL1EGCrystalClusters +
    process.l1tPFHcalProducerFromTPDigis +
    process.l1tPFHGCalProducerFrom3DTPs +
    process.l1tPFHGCalProducerFromTriggerCells +
    process.l1tPFTkProducersFromL1Tracks +
    process.l1tPFMuProducerFromL1Mu
)

process.p = cms.Path(process.s)

process.out = cms.OutputModule("PoolOutputModule",
        fileName = cms.untracked.string("inputs_17D.root"),
        outputCommands = cms.untracked.vstring("drop *",
            "keep *_genParticles_*_*",
            "keep *_ak4GenJetsNoNu_*_*",
            "keep *_genMetTrue_*_*",
            "keep *_l1tPF*_*_IN",
        ),
        compressionAlgorithm = cms.untracked.string('LZMA'),
        compressionLevel = cms.untracked.int32(4),
        dropMetaData = cms.untracked.string('ALL'),
        fastCloning = cms.untracked.bool(False),
        overrideInputFileSplitLevels = cms.untracked.bool(True),
        eventAutoFlushCompressedSize = cms.untracked.int32(15728640),
)
process.e = cms.EndPath(process.out)



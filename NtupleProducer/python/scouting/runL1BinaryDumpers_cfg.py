import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras

process = cms.Process("L1Dump", eras.Phase2C9)

process.load('Configuration.StandardSequences.Services_cff')
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.MessageLogger.cerr.FwkReport.reportEvery = 10

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:/tmp/gpetrucc/005E74D6-B50E-674E-89E6-EAA9A617B476.root')
)

process.load('Configuration.Geometry.GeometryExtended2026D49Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2026D49_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff') # needed to read HCal TPs
process.load('SimCalorimetry.HGCalSimProducers.hgcalDigitizer_cfi') # needed for HGCAL_noise_fC
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '123X_mcRun4_realistic_v3', '')

process.pfClustersFromCombinedCaloHCal.phase2barrelCaloTowers = [cms.InputTag("L1EGammaClusterEmuProducer",)]


process.puppiDump = cms.EDAnalyzer("L1PuppiBinaryDumper",
                src = cms.InputTag("l1ctLayer1","Puppi"),
                outName = cms.string("puppi.dump"))

process.jetDump = cms.EDAnalyzer("L1JetBinaryDumper",
                src = cms.InputTag("scPFL1PuppiCorrectedEmulator"),
                ptMin = cms.double(15),
                outName = cms.string("puppiJets.dump"))

process.deps = cms.Task(
    process.L1SAMuonsGmt,
    process.L1GTTInputProducer,
    process.L1VertexFinderEmulator,
    process.l1ctLayer1TaskInputsTask,
    process.l1ctLayer1Task,
    process.l1ctLayer2Deregionizer, 
    process.scPFL1PuppiEmulator, 
    process.scPFL1PuppiCorrectedEmulator,
)

process.run = cms.EndPath(
    process.puppiDump +
    process.jetDump
)
process.run.associate(process.deps)

def noPU():
    for X in "", "EM", "Raw", "EMRaw":
        pfc = getattr(process, 'pfClustersFromHGC3DClusters'+X, None)
        if not pfc: continue
        pfc.emVsPUID.wp = "-1.0"

process.schedule = cms.Schedule(process.run)

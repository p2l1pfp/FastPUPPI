import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras

process = cms.Process("PF", eras.phase2_trigger)
process.load('Configuration.StandardSequences.Services_cff')
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True), allowUnscheduled = cms.untracked.bool(False) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1))
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

runCalo = True

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:test.root'),
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
)

process.load("L1Trigger.Phase2L1ParticleFlow.l1ParticleFlow_cff")

if runCalo:
    process.load('Configuration.Geometry.GeometryExtended2023D17Reco_cff')
    process.load('Configuration.StandardSequences.MagneticField_cff')
    process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff') # needed to read HCal TPs
    process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
    from Configuration.AlCa.GlobalTag import GlobalTag
    process.GlobalTag = GlobalTag(process.GlobalTag, '100X_upgrade2023_realistic_v1', '')
    process.l1ParticleFlow.remove(process.l1EGammaCrystalsProducer)
    process.runPF = cms.Sequence( process.l1ParticleFlow )

else:
    process.runPF = cms.Sequence( process.l1pfProducer )

process.ntuple = cms.EDAnalyzer("ResponseNTuplizer",
    genJets = cms.InputTag("ak4GenJetsNoNu"),
    genParticles = cms.InputTag("genParticles"),
    isParticleGun = cms.bool(False),
    doRandom = cms.bool(False),
    objects = cms.PSet(
       #TPEcal = cms.VInputTag('l1tPFHGCalProducerFrom3DTPsEM', 'l1tPFEcalProducerFromL1EGCrystalClusters', ),
       #TPHcal = cms.VInputTag('l1tPFHcalProducerFromTPDigis', 'l1tPFHGCalProducerFromTriggerCells:towersFHBH', ),
       #TPCalo = cms.VInputTag('l1tPFHGCalProducerFrom3DTPsEM', 'l1tPFEcalProducerFromL1EGCrystalClusters', 'l1tPFHcalProducerFromTPDigis', 'l1tPFHGCalProducerFromTriggerCells:towersFHBH',), 
       #TPTK   = cms.VInputTag('l1tPFTkProducersFromL1Tracks',),
       ## -- processed --
       #L1RawEcal = cms.VInputTag(cms.InputTag('InfoOut','RawEmCalo')),
       #L1RawCalo = cms.VInputTag(cms.InputTag('InfoOut','RawCalo')),
       #L1Ecal = cms.VInputTag(cms.InputTag('InfoOut','EmCalo')),
       #L1Calo = cms.VInputTag("InfoOut:Calo",),
       #L1TK = cms.VInputTag("InfoOut:TK",),
       #L1TKV = cms.VInputTag("InfoOut:TKVtx",),
       #L1PF = cms.VInputTag("InfoOut:PF",),
       #L1Puppi = cms.VInputTag("InfoOut:Puppi",),
        #---
        NewPFTK = cms.VInputTag('pfTracksFromL1Tracks',),
        NewEcal = cms.VInputTag('pfClustersFromL1EGClusters', 'pfClustersFromHGC3DClustersEM'),
        NewL1CRawEcal = cms.VInputTag('pfClustersFromCombinedCalo:emUncalibrated'),
        NewL1CEcal    = cms.VInputTag('pfClustersFromCombinedCalo:emCalibrated'),
        NewL1CRawCalo = cms.VInputTag('pfClustersFromCombinedCalo:uncalibrated'),
        NewL1CCalo    = cms.VInputTag('pfClustersFromCombinedCalo:calibrated'),
        NewL1Ecal = cms.VInputTag('l1pfProducer:EmCalo'),
        NewL1Calo = cms.VInputTag('l1pfProducer:Calo'),
        NewL1TK = cms.VInputTag("l1pfProducer:TK",),
        NewL1TKV = cms.VInputTag("l1pfProducer:TKVtx",),
        NewL1PF = cms.VInputTag("l1pfProducer:PF",),
        NewL1Puppi = cms.VInputTag("l1pfProducer:Puppi",),
       #OldL1CRawEcal = cms.VInputTag('CaloInfoOut:emUncalibrated'),
       #OldL1CEcal    = cms.VInputTag('CaloInfoOut:emCalibrated'),
       #OldL1CRawCalo = cms.VInputTag('CaloInfoOut:uncalibrated'),
       #OldL1CCalo    = cms.VInputTag('CaloInfoOut:calibrated'),
    ),
    copyUInts = cms.VInputTag(),
)

process.p = cms.Path(
    process.runPF +
    process.ntuple
)

process.TFileService = cms.Service("TFileService", fileName = cms.string("testTupleNew.root"))

def goGun():
    process.ntuple.isParticleGun = True

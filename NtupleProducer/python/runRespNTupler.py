import FWCore.ParameterSet.Config as cms

process = cms.Process("RESP")

process.load('Configuration.StandardSequences.Services_cff')
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True), allowUnscheduled = cms.untracked.bool(False) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1))
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:l1pf_out.root')
)

process.ntuple = cms.EDAnalyzer("ResponseNTuplizer",
    genJets = cms.InputTag("ak4GenJetsNoNu"),
    genParticles = cms.InputTag("genParticles"),
    # -- inputs --
    Ecal = cms.VInputTag('l1tPFEcalProducerFromOfflineRechits:towers','l1tPFHGCalEEProducerFromOfflineRechits:towers', 'l1tPFHFProducerFromOfflineRechits:towers'),
    Hcal = cms.VInputTag('l1tPFHcalProducerFromOfflineRechits:towers','l1tPFHGCalFHProducerFromOfflineRechits:towers', 'l1tPFHGCalBHProducerFromOfflineRechits:towers', 'l1tPFHFProducerFromOfflineRechits:towers'),
    Calo = cms.VInputTag('l1tPFEcalProducerFromOfflineRechits:towers','l1tPFHGCalEEProducerFromOfflineRechits:towers', 'l1tPFHcalProducerFromOfflineRechits:towers', 'l1tPFHGCalFHProducerFromOfflineRechits:towers', 'l1tPFHGCalBHProducerFromOfflineRechits:towers', 'l1tPFHFProducerFromOfflineRechits:towers'),
    TK   = cms.VInputTag('l1tPFTkProducersFromOfflineTracksStrips'),
    # -- processed --
    L1RawCalo = cms.VInputTag("InfoOut:RawCalo",),
    L1NewCalo = cms.VInputTag("InfoOut:NewCalo",),
    L1Calo = cms.VInputTag("InfoOut:Calo",),
    L1TK = cms.VInputTag("InfoOut:TK",),
    L1PF = cms.VInputTag("InfoOut:PF",),
    L1Puppi = cms.VInputTag("InfoOut:Puppi",),
   ## -- clustered --
   #L1ak4RawCalo = cms.VInputTag("ak4L1RawCalo",),
   #L1ak4Calo = cms.VInputTag("ak4L1Calo",),
   #L1ak4TK = cms.VInputTag("ak4L1TK",),
   #L1ak4PF = cms.VInputTag("ak4L1PF",),
   #L1ak4Puppi = cms.VInputTag("ak4L1Puppi",),
)
process.p = cms.Path(process.ntuple)
process.TFileService = cms.Service("TFileService", fileName = cms.string("respTupleNew.root"))
if True:
    process.load('FastPUPPI.NtupleProducer.ntupleProducer_cfi')
    process.InfoOut.outputName = ""; # turn off Ntuples
    process.p = cms.Path(process.InfoOut + process.ntuple)
    

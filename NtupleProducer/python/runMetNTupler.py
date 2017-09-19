import FWCore.ParameterSet.Config as cms

process = cms.Process("RESP")

process.load('Configuration.StandardSequences.Services_cff')
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True), allowUnscheduled = cms.untracked.bool(False) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1))
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:/eos/cms/store/cmst3/user/jngadiub/L1PFInputs/TTbar_PU140/inputs_17D_TTbar_PU140_job1.root'),
)
process.source.duplicateCheckMode = cms.untracked.string("noDuplicateCheck")

process.load('FastPUPPI.NtupleProducer.l1tPFHGCalProducerFrom3DTPsEM_cfi')
process.load('FastPUPPI.NtupleProducer.caloNtupleProducer_cfi')
process.load('FastPUPPI.NtupleProducer.ntupleProducer_cfi')
process.CaloInfoOut.outputName = ""; # turn off Ntuples
process.InfoOut.outputName = ""; # turn off Ntuples

from RecoMET.METProducers.PFMET_cfi import pfMet
pfMet.calculateSignificance = False
process.l1MetRawCalo = pfMet.clone(src = "InfoOut:RawCalo")
process.l1MetCalo    = pfMet.clone(src = "InfoOut:Calo")
process.l1MetTK      = pfMet.clone(src = "InfoOut:TK")
process.l1MetTKV      = pfMet.clone(src = "InfoOut:TKVtx")
process.l1MetPF      = pfMet.clone(src = "InfoOut:PF")
process.l1MetPuppi   = pfMet.clone(src = "InfoOut:Puppi")

process.ntuple = cms.EDAnalyzer("MetNTuplizer",
    mets = cms.PSet(
        Gen = cms.InputTag("genMetTrue"),
        RawCalo = cms.InputTag("l1MetRawCalo"),
        Calo = cms.InputTag("l1MetCalo"),
        TK = cms.InputTag("l1MetTK"),
        TKV = cms.InputTag("l1MetTKV"),
        PF = cms.InputTag("l1MetPF"),
        Puppi = cms.InputTag("l1MetPuppi"),
    ),
)

process.mets = cms.Sequence( process.l1MetRawCalo + process.l1MetCalo + process.l1MetTK + process.l1MetTKV + process.l1MetPF + process.l1MetPuppi)
process.p = cms.Path(process.l1tPFHGCalProducerFrom3DTPsEM + process.CaloInfoOut + process.InfoOut + process.mets + process.ntuple)
process.TFileService = cms.Service("TFileService", fileName = cms.string("metTuple.root"))

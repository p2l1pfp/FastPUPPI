import FWCore.ParameterSet.Config as cms

process = cms.Process("RESP")

process.load('Configuration.StandardSequences.Services_cff')
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True), allowUnscheduled = cms.untracked.bool(False) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1))
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:l1pf_out_TTbar_PU140_job1.root')
)

process.load('FastPUPPI.NtupleProducer.ntupleProducer_cfi')
process.load('FastPUPPI.NtupleProducer.l1tPFMuProducerFromL1Mu_cfi')
process.InfoOut.outputName = ""; # turn off Ntuples

from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets
process.ak4L1RawCalo = ak4PFJets.clone(src = 'InfoOut:RawCalo')
process.ak4L1Calo    = ak4PFJets.clone(src = 'InfoOut:Calo')
process.ak4L1TK      = ak4PFJets.clone(src = 'InfoOut:TK')
process.ak4L1PF      = ak4PFJets.clone(src = 'InfoOut:PF')
process.ak4L1Puppi   = ak4PFJets.clone(src = 'InfoOut:Puppi')
process.ak4L1ICalo   = ak4PFJets.clone(src = 'InfoOut:L1Calo')
process.ak4L1ITK     = ak4PFJets.clone(src = 'InfoOut:L1TK')
process.ak4L1IPF     = ak4PFJets.clone(src = 'InfoOut:L1PF')
process.ak4L1IPuppi  = ak4PFJets.clone(src = 'InfoOut:L1Puppi')

process.jets = cms.Sequence( 
    process.ak4L1RawCalo + process.ak4L1Calo + process.ak4L1TK + process.ak4L1PF + process.ak4L1Puppi + 
    process.ak4L1ICalo + process.ak4L1ITK + process.ak4L1IPF + process.ak4L1IPuppi 
)

process.ntuple = cms.EDAnalyzer("JetNTuplizer",
    jets = cms.PSet(
        Gen = cms.InputTag("ak4GenJetsNoNu"),
        RawCalo = cms.InputTag("ak4L1RawCalo"),
        Calo = cms.InputTag("ak4L1Calo"),
        TK = cms.InputTag("ak4L1TK"),
        PF = cms.InputTag("ak4L1PF"),
        Puppi = cms.InputTag("ak4L1Puppi"),
        L1Calo = cms.InputTag("ak4L1ICalo"),
        L1TK = cms.InputTag("ak4L1ITK"),
        L1PF = cms.InputTag("ak4L1IPF"),
        L1Puppi = cms.InputTag("ak4L1IPuppi"),
    ),
    sels = cms.PSet(
        E24Pt15 = cms.string("pt > 15 && abs(eta) < 2.4"),
        E24Pt20 = cms.string("pt > 20 && abs(eta) < 2.4"),
        E24Pt30 = cms.string("pt > 30 && abs(eta) < 2.4"),
        E24Pt40 = cms.string("pt > 40 && abs(eta) < 2.4"),
        E24Pt50 = cms.string("pt > 50 && abs(eta) < 2.4"),
        E24Pt60 = cms.string("pt > 60 && abs(eta) < 2.4"),
        E24Pt80 = cms.string("pt > 80 && abs(eta) < 2.4"),
        E30Pt15 = cms.string("pt > 15 && abs(eta) < 3.0"),
        E30Pt20 = cms.string("pt > 20 && abs(eta) < 3.0"),
        E30Pt30 = cms.string("pt > 30 && abs(eta) < 3.0"),
        E30Pt40 = cms.string("pt > 40 && abs(eta) < 3.0"),
        E30Pt50 = cms.string("pt > 50 && abs(eta) < 3.0"),
        E30Pt60 = cms.string("pt > 60 && abs(eta) < 3.0"),
        E30Pt80 = cms.string("pt > 80 && abs(eta) < 3.0"),
        E47Pt15 = cms.string("pt > 15 && abs(eta) < 4.7"),
        E47Pt20 = cms.string("pt > 20 && abs(eta) < 4.7"),
        E47Pt30 = cms.string("pt > 30 && abs(eta) < 4.7"),
        E47Pt40 = cms.string("pt > 40 && abs(eta) < 4.7"),
        E47Pt50 = cms.string("pt > 50 && abs(eta) < 4.7"),
        E47Pt60 = cms.string("pt > 50 && abs(eta) < 4.7"),
        E47Pt80 = cms.string("pt > 80 && abs(eta) < 4.7"),
    )
)

process.p = cms.Path(process.l1tPFMuProducerFromL1Mu + process.InfoOut + process.jets + process.ntuple)
process.TFileService = cms.Service("TFileService", fileName = cms.string("jetTupleNew.root"))

 


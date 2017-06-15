import FWCore.ParameterSet.Config as cms

process = cms.Process("RESP")

process.load('Configuration.StandardSequences.Services_cff')
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True), allowUnscheduled = cms.untracked.bool(False) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1))
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('/store/cmst3/user/gpetrucc/l1phase2/Spring17D/200517/inputs_17D_TTbar_PU140_job4.root')
)
process.source.duplicateCheckMode = cms.untracked.string("noDuplicateCheck")

process.load('FastPUPPI.NtupleProducer.caloNtupleProducer_cfi')
process.load('FastPUPPI.NtupleProducer.ntupleProducer_cfi')
process.CaloInfoOut.outputName = ""; # turn off Ntuples
process.InfoOut.outputName = ""; # turn off Ntuples

from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets
process.ak4L1RawCalo = ak4PFJets.clone(src = 'InfoOut:RawCalo')
process.ak4L1Calo    = ak4PFJets.clone(src = 'InfoOut:Calo')
process.ak4L1TK      = ak4PFJets.clone(src = 'InfoOut:TK')
process.ak4L1TKV     = ak4PFJets.clone(src = 'InfoOut:TKVtx')
process.ak4L1PF      = ak4PFJets.clone(src = 'InfoOut:PF')
process.ak4L1Puppi   = ak4PFJets.clone(src = 'InfoOut:Puppi')
process.ak4L1ICalo   = ak4PFJets.clone(src = 'InfoOut:L1Calo')
process.ak4L1ITK     = ak4PFJets.clone(src = 'InfoOut:L1TK')
process.ak4L1IPF     = ak4PFJets.clone(src = 'InfoOut:L1PF')
process.ak4L1IPuppi  = ak4PFJets.clone(src = 'InfoOut:L1Puppi')

process.jets = cms.Sequence( 
    process.ak4L1RawCalo + process.ak4L1Calo + process.ak4L1TK + process.ak4L1TKV + process.ak4L1PF + process.ak4L1Puppi + 
    process.ak4L1ICalo + process.ak4L1ITK + process.ak4L1IPF + process.ak4L1IPuppi 
)

JEC_PU140 = {
    'L1Calo' : cms.PSet(
                        etaBins = cms.vdouble( 0.500,  1.000,  1.500,  2.000,  2.500,  3.000,  3.500,  4.000,  4.500,  5.000),
                        offset  = cms.vdouble( 12.164,  15.991,  18.123,  29.138,  40.576,  37.452,  48.410,  60.027,  45.062,  53.013),
                        scale   = cms.vdouble( 0.938,  0.922,  0.881,  0.818,  0.947,  1.023,  0.748,  0.765,  1.213,  1.276),
            ),
    'L1TK' : cms.PSet(
                        etaBins = cms.vdouble( 0.500,  1.000,  1.500,  2.000,  2.500),
                        offset  = cms.vdouble( 0.440, -0.510,  0.974, -2.687,  0.977),
                        scale   = cms.vdouble( 0.593,  0.580,  0.553,  0.714,  0.537),
            ),
    'L1TKV' : cms.PSet(
                        etaBins = cms.vdouble( 0.500,  1.000,  1.500,  2.000,  2.500),
                        offset  = cms.vdouble(-1.546, -1.858,  0.633,  1.144,  0.618),
                        scale   = cms.vdouble( 0.507,  0.497,  0.395,  0.314,  0.253),
            ),
    'L1PF' : cms.PSet(
                        etaBins = cms.vdouble( 0.500,  1.000,  1.500,  2.000,  2.500,  3.000,  3.500,  4.000,  4.500,  5.000),
                        offset  = cms.vdouble( 12.622,  15.549,  18.607,  26.369,  41.212,  38.112,  48.410,  60.027,  45.062,  53.013),
                        scale   = cms.vdouble( 0.988,  0.967,  0.928,  0.926,  0.986,  1.025,  0.748,  0.765,  1.213,  1.276),
            ),
    'L1Puppi' : cms.PSet(
                        etaBins = cms.vdouble( 0.500,  1.000,  1.500,  2.000,  2.500,  3.000,  3.500,  4.000,  4.500,  5.000),
                        offset  = cms.vdouble(-10.071, -10.927, -7.854, -2.480, -7.974, -2.991,  1.054,  7.412,  0.000,  0.000),
                        scale   = cms.vdouble( 0.922,  0.959,  0.821,  0.571,  0.769,  0.363,  0.289,  0.226,  1.000,  1.000),
            ),
}
JEC = JEC_PU140;

process.ntuple = cms.EDAnalyzer("JetNTuplizer",
    jets = cms.PSet(
        Gen = cms.InputTag("ak4GenJetsNoNu"),
        #RawCalo = cms.InputTag("ak4L1RawCalo"),
        Calo = cms.InputTag("ak4L1Calo"),
        TK = cms.InputTag("ak4L1TK"),
        TKV = cms.InputTag("ak4L1TKV"),
        PF = cms.InputTag("ak4L1PF"),
        Puppi = cms.InputTag("ak4L1Puppi"),
        L1Calo = cms.InputTag("ak4L1ICalo"),
        L1TK = cms.InputTag("ak4L1ITK"),
        L1PF = cms.InputTag("ak4L1IPF"),
        L1Puppi = cms.InputTag("ak4L1IPuppi"),
    ),
    jecs = cms.PSet(
        Calo = JEC['L1Calo'],
        TK = JEC['L1TK'],
        TKV = JEC['L1TKV'],
        PF = JEC['L1PF'],
        Puppi = JEC['L1Puppi'],
        L1Calo = JEC['L1Calo'],
        L1TK = JEC['L1TK'],
        L1PF = JEC['L1PF'],
        L1Puppi = JEC['L1Puppi'],
    ),
    sels = cms.PSet(
        #E24Pt15 = cms.string("pt > 15 && abs(eta) < 2.4"),
        E24Pt20 = cms.string("pt > 20 && abs(eta) < 2.4"),
        E24Pt30 = cms.string("pt > 30 && abs(eta) < 2.4"),
        E24Pt40 = cms.string("pt > 40 && abs(eta) < 2.4"),
        E24Pt50 = cms.string("pt > 50 && abs(eta) < 2.4"),
        #E24Pt60 = cms.string("pt > 60 && abs(eta) < 2.4"),
        #E24Pt80 = cms.string("pt > 80 && abs(eta) < 2.4"),
        #E30Pt15 = cms.string("pt > 15 && abs(eta) < 3.0"),
        E30Pt20 = cms.string("pt > 20 && abs(eta) < 3.0"),
        E30Pt30 = cms.string("pt > 30 && abs(eta) < 3.0"),
        E30Pt40 = cms.string("pt > 40 && abs(eta) < 3.0"),
        E30Pt50 = cms.string("pt > 50 && abs(eta) < 3.0"),
        #E30Pt60 = cms.string("pt > 60 && abs(eta) < 3.0"),
        #E30Pt80 = cms.string("pt > 80 && abs(eta) < 3.0"),
        #E47Pt15 = cms.string("pt > 15 && abs(eta) < 4.7"),
        #E47Pt20 = cms.string("pt > 20 && abs(eta) < 4.7"),
        E47Pt30 = cms.string("pt > 30 && abs(eta) < 4.7"),
        E47Pt40 = cms.string("pt > 40 && abs(eta) < 4.7"),
        E47Pt50 = cms.string("pt > 50 && abs(eta) < 4.7"),
        #E47Pt60 = cms.string("pt > 50 && abs(eta) < 4.7"),
        #E47Pt80 = cms.string("pt > 80 && abs(eta) < 4.7"),
    )
)

process.p = cms.Path(process.CaloInfoOut + process.InfoOut + process.jets + process.ntuple)
process.TFileService = cms.Service("TFileService", fileName = cms.string("jetTupleNew.root"))

 


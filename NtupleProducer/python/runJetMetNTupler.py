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
#tighter track cuts as in track trigger TDR
#process.InfoOut.trkMinStubs = cms.uint32(5)
#process.InfoOut.trkMaxChi2 = cms.double(20)
	 
from RecoMET.METProducers.PFMET_cfi import pfMet
pfMet.calculateSignificance = False
process.l1MetRawCalo = pfMet.clone(src = "InfoOut:RawCalo")
process.l1MetCalo    = pfMet.clone(src = "InfoOut:Calo")
process.l1MetTK      = pfMet.clone(src = "InfoOut:TK")
process.l1MetTKV     = pfMet.clone(src = "InfoOut:TKVtx")
process.l1MetPF      = pfMet.clone(src = "InfoOut:PF")
process.l1MetPuppi   = pfMet.clone(src = "InfoOut:Puppi")

#process.METntuple = cms.EDAnalyzer("MetNTuplizer",
#    mets = cms.PSet(
#        Gen = cms.InputTag("genMetTrue"),
#        RawCalo = cms.InputTag("l1MetRawCalo"),
#        Calo = cms.InputTag("l1MetCalo"),
#        TK = cms.InputTag("l1MetTK"),
#        TKV = cms.InputTag("l1MetTKV"),
#        PF = cms.InputTag("l1MetPF"),
#        Puppi = cms.InputTag("l1MetPuppi"),
#    ),
#)

process.mets = cms.Sequence( process.l1MetRawCalo + process.l1MetCalo + process.l1MetTK + process.l1MetTKV + process.l1MetPF + process.l1MetPuppi)

from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets
process.ak4L1RawCalo = ak4PFJets.clone(src = 'InfoOut:RawCalo')
process.ak4L1Calo    = ak4PFJets.clone(src = 'InfoOut:Calo')
process.ak4L1TK      = ak4PFJets.clone(src = 'InfoOut:TK')
process.ak4L1TKV     = ak4PFJets.clone(src = 'InfoOut:TKVtx')
process.ak4L1PF      = ak4PFJets.clone(src = 'InfoOut:PF')
process.ak4L1Puppi   = ak4PFJets.clone(src = 'InfoOut:Puppi')

process.jets = cms.Sequence(     
    process.ak4L1RawCalo + process.ak4L1Calo + process.ak4L1TK + process.ak4L1TKV + process.ak4L1PF + process.ak4L1Puppi
)

JEC_PU140_pt2_trkCutsLoose = {
    'L1Calo' : cms.PSet(
			etaBins = cms.vdouble( 0.500,  1.000,  1.500,  2.000,  2.500,  3.000,  3.500,  4.000,  4.500,  5.000),
			offset  = cms.vdouble( 0.214,  2.048,  4.319,  14.301,  20.314,  16.401,  38.697,  49.053,  48.612,  47.783),
			scale   = cms.vdouble( 0.993,  0.990,  0.909,  0.931,  1.005,  0.993,  0.808,  0.921,  0.987,  1.312),
            ),
    'L1TK' : cms.PSet(
			etaBins = cms.vdouble( 0.500,  1.000,  1.500,  2.000,  2.500),
			offset  = cms.vdouble( 1.054,  1.465,  2.814,  3.392,  5.913),
			scale   = cms.vdouble( 0.557,  0.536,  0.481,  0.477,  0.397),
            ),
    'L1TKV' : cms.PSet(
			etaBins = cms.vdouble( 0.500,  1.000,  1.500,  2.000,  2.500),
			offset  = cms.vdouble(-0.675, -0.821,  0.434,  1.229,  3.230),
			scale   = cms.vdouble( 0.492,  0.487,  0.441,  0.437,  0.375),
            ),
    'L1PF' : cms.PSet(
			etaBins = cms.vdouble( 0.500,  1.000,  1.500,  2.000,  2.500,  3.000,  3.500,  4.000,  4.500,  5.000),
			offset  = cms.vdouble( 1.445,  2.989,  5.236,  15.621,  27.823,  18.562,  38.697,  49.053,  48.612,  47.783),
			scale   = cms.vdouble( 1.078,  1.072,  0.993,  1.046,  1.078,  1.003,  0.808,  0.921,  0.987,  1.312),
            ),
    'L1Puppi' : cms.PSet(
			etaBins = cms.vdouble( 0.500,  1.000,  1.500,  2.000,  2.500,  3.000,  3.500,  4.000,  4.500,  5.000),
			offset  = cms.vdouble(-12.054, -12.002, -10.493, -10.274, -8.143,  4.279,  4.964,  11.538,  10.682, -7.683),
			scale   = cms.vdouble( 1.042,  1.048,  0.966,  0.982,  1.070,  0.367,  0.885,  0.986,  1.097,  1.520),
            ),
}

JEC_PU140_pt2_trkCutsTight = {
    'L1Calo' : cms.PSet(
			etaBins = cms.vdouble( 0.500,  1.000,  1.500,  2.000,  2.500,  3.000,  3.500,  4.000,  4.500,  5.000),
			offset  = cms.vdouble( 0.214,  2.048,  4.319,  14.301,  20.314,  16.401,  38.697,  49.053,  48.612,  47.783),
			scale   = cms.vdouble( 0.993,  0.990,  0.909,  0.931,  1.005,  0.993,  0.808,  0.921,  0.987,  1.312),
            ),
    'L1TK' : cms.PSet(
			etaBins = cms.vdouble( 0.500,  1.000,  1.500,  2.000,  2.500),
			offset  = cms.vdouble( 1.329,  1.607,  3.287,  4.462,  8.821),
			scale   = cms.vdouble( 0.429,  0.403,  0.364,  0.401,  0.238),
            ),
    'L1TKV' : cms.PSet(
			etaBins = cms.vdouble( 0.500,  1.000,  1.500,  2.000,  2.500),
			offset  = cms.vdouble(-0.704, -0.178,  1.115,  1.921,  6.593),
			scale   = cms.vdouble( 0.414,  0.388,  0.354,  0.389,  0.232),
            ),
    'L1PF' : cms.PSet(
			etaBins = cms.vdouble( 0.500,  1.000,  1.500,  2.000,  2.500,  3.000,  3.500,  4.000,  4.500,  5.000),
			offset  = cms.vdouble(-0.917,  0.708,  3.730,  15.081,  27.327,  18.220,  38.697,  49.053,  48.612,  47.783),
			scale   = cms.vdouble( 1.066,  1.060,  0.977,  1.036,  1.070,  1.002,  0.808,  0.921,  0.987,  1.312),
            ),
    'L1Puppi' : cms.PSet(
			etaBins = cms.vdouble( 0.500,  1.000,  1.500,  2.000,  2.500,  3.000,  3.500,  4.000,  4.500,  5.000),
			offset  = cms.vdouble(-16.913, -16.867, -14.144, -12.415, -10.023,  4.423,  4.915,  11.427,  10.613, -7.733),
			scale   = cms.vdouble( 1.058,  1.053,  0.952,  0.974,  1.016,  0.346,  0.884,  0.987,  1.094,  1.520),
            ),
}

JEC_PU140_pt2 = {
    'L1Calo' : cms.PSet(
			etaBins = cms.vdouble( 0.500,  1.000,  1.500,  2.000,  2.500,  3.000,  3.500,  4.000,  4.500,  5.000),
			offset  = cms.vdouble( 0.225,  2.078,  4.325,  14.328,  20.329,  16.589,  38.821,  49.950,  50.680,  49.587),
			scale   = cms.vdouble( 0.993,  0.990,  0.909,  0.931,  1.005,  0.992,  0.807,  0.905,  0.948,  1.286),
            ),
    'L1TK' : cms.PSet(
			etaBins = cms.vdouble( 0.500,  1.000,  1.500,  2.000,  2.500),
			offset  = cms.vdouble( 0.428,  1.030,  1.410, -3.163,  4.369),
			scale   = cms.vdouble( 0.604,  0.579,  0.565,  0.747,  0.501),
            ),
    'L1TKV' : cms.PSet(
			etaBins = cms.vdouble( 0.500,  1.000,  1.500,  2.000,  2.500),
			offset  = cms.vdouble(-0.778, -0.889,  0.081, -0.034,  2.775),
			scale   = cms.vdouble( 0.509,  0.504,  0.471,  0.511,  0.403),
            ),
    'L1PF' : cms.PSet(
			etaBins = cms.vdouble( 0.500,  1.000,  1.500,  2.000,  2.500,  3.000,  3.500,  4.000,  4.500,  5.000),
			offset  = cms.vdouble( 1.749,  3.472,  5.346,  14.246, 28.300,  18.998,  38.821,  49.950,  50.680,  49.587),
			scale   = cms.vdouble( 1.090,  1.080,  1.021,  1.142,  1.091,  0.999,  0.807,  0.905,  0.948,  1.286),
            ),
    'L1Puppi' : cms.PSet(
			etaBins = cms.vdouble( 0.500,  1.000,  1.500,  2.000,  2.500,  3.000,  3.500,  4.000,  4.500,  5.000),
			offset  = cms.vdouble(-11.052, -10.943, -8.822, -5.544, -6.385,  3.602,  5.024,  11.132,  11.561, -0.803),
			scale   = cms.vdouble( 1.028,  1.033,  0.938,  0.889,  1.016,  0.382,  0.887,  0.994,  1.076,  1.379),
            ),
}

JEC_PU140_pt3 = {
    'L1Calo' : cms.PSet(
			etaBins = cms.vdouble( 0.500,  1.000,  1.500,  2.000,  2.500,  3.000,  3.500,  4.000,  4.500,  5.000),
			offset  = cms.vdouble( 0.214,  2.048,  4.319,  14.301,  20.314,  16.401,  38.697,  49.053,  48.612,  47.783),
			scale   = cms.vdouble( 0.993,  0.990,  0.909,  0.931,  1.005,  0.993,  0.808,  0.921,  0.987,  1.312),
            ),
    'L1TK' : cms.PSet(
			etaBins = cms.vdouble( 0.500,  1.000,  1.500,  2.000,  2.500),
			offset  = cms.vdouble(-3.279, -2.823, -1.039, -0.378,  2.416),
			scale   = cms.vdouble( 0.548,  0.531,  0.473,  0.466,  0.389),
            ),
    'L1TKV' : cms.PSet(
			etaBins = cms.vdouble( 0.500,  1.000,  1.500,  2.000,  2.500),
			offset  = cms.vdouble(-2.318, -2.450, -1.001, -0.459,  1.542),
			scale   = cms.vdouble( 0.484,  0.480,  0.432,  0.429,  0.369),
            ),
    'L1PF' : cms.PSet(
			etaBins = cms.vdouble( 0.500,  1.000,  1.500,  2.000,  2.500,  3.000,  3.500,  4.000,  4.500,  5.000),
			offset  = cms.vdouble(-2.264, -0.342,  2.326,  12.990,  26.232,  18.224,  38.697,  49.053,  48.612,  47.783),
			scale   = cms.vdouble( 1.078,  1.072,  0.990,  1.047,  1.074,  1.001,  0.808,  0.921,  0.987,  1.312),
            ),
    'L1Puppi' : cms.PSet(
			etaBins = cms.vdouble( 0.500,  1.000,  1.500,  2.000,  2.500,  3.000,  3.500,  4.000,  4.500,  5.000),
			offset  = cms.vdouble(-14.675, -14.877, -12.686, -12.862, -14.130,  3.085,  4.940,  11.376,  10.500, -7.733),
			scale   = cms.vdouble( 0.889,  0.890,  0.805,  0.812,  0.890,  0.383,  0.884,  0.986,  1.097,  1.520),
            ),
}

JEC = JEC_PU140_pt2_trkCutsLoose;

process.ntuple = cms.EDAnalyzer("JetMetNTuplizer",
    jets = cms.PSet(
        AK4GenJets = cms.InputTag("ak4GenJetsNoNu"),
        AK4RawCaloJets = cms.InputTag("ak4L1RawCalo"),
        AK4CaloJets = cms.InputTag("ak4L1Calo"),
        AK4TKJets = cms.InputTag("ak4L1TK"),
        AK4TKVJets = cms.InputTag("ak4L1TKV"),
        AK4PFJets = cms.InputTag("ak4L1PF"),
        AK4PuppiJets = cms.InputTag("ak4L1Puppi"),
    ),
    jecs = cms.PSet(
        AK4CaloJets = JEC['L1Calo'],
        AK4TKJets = JEC['L1TK'],
        AK4TKVJets = JEC['L1TKV'],
        AK4PFJets = JEC['L1PF'],
        AK4PuppiJets = JEC['L1Puppi'],
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
    ),
    mets = cms.PSet(
        METGen = cms.InputTag("genMetTrue"),
        METRawCalo = cms.InputTag("l1MetRawCalo"),
        METCalo = cms.InputTag("l1MetCalo"),
        METTK = cms.InputTag("l1MetTK"),
        METTKV = cms.InputTag("l1MetTKV"),
        METPF = cms.InputTag("l1MetPF"),
        METPuppi = cms.InputTag("l1MetPuppi"),
    ),
)


process.p = cms.Path(process.l1tPFHGCalProducerFrom3DTPsEM + process.CaloInfoOut + process.InfoOut + process.mets + process.jets + process.ntuple)
process.TFileService = cms.Service("TFileService", fileName = cms.string("jetmetTuple.root"))

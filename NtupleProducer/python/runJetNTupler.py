import FWCore.ParameterSet.Config as cms

process = cms.Process("RESP")

process.load('Configuration.StandardSequences.Services_cff')
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True), allowUnscheduled = cms.untracked.bool(False) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1))
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:/eos/cms/store/cmst3/user/jngadiub/L1PFInputs/SingleNeutrino_PU140/inputs_17D_SingleNeutrino_PU140_job1.root'),
    eventsToProcess = cms.untracked.VEventRange('1:973:48615'),
)
process.source.duplicateCheckMode = cms.untracked.string("noDuplicateCheck")

process.load('FastPUPPI.NtupleProducer.l1tPFHGCalProducerFrom3DTPsEM_cfi')
process.load('FastPUPPI.NtupleProducer.caloNtupleProducer_cfi')
process.load('FastPUPPI.NtupleProducer.ntupleProducer_cfi')
process.CaloInfoOut.outputName = ""; # turn off Ntuples
process.InfoOut.outputName = ""; # turn off Ntuples
process.InfoOut.debug = 2;

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
    process.ak4L1RawCalo + process.ak4L1Calo + process.ak4L1TK + process.ak4L1TKV + process.ak4L1PF + process.ak4L1Puppi
    #process.ak4L1ICalo + process.ak4L1ITK + process.ak4L1IPF + process.ak4L1IPuppi 
)

JEC_PU140 = {
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

JEC_PU0 = {
    'L1Calo' : cms.PSet(
			etaBins = cms.vdouble( 0.500,  1.000,  1.500,  2.000,  2.500,  3.000,  3.500,  4.000,  4.500,  5.000),
			offset  = cms.vdouble(-12.704, -12.376, -10.322, -9.258, -10.556, -11.380, -4.990, -6.409, -2.709, -9.033),
			scale   = cms.vdouble( 0.918,  0.911,  0.818,  0.852,  0.931,  0.899,  0.734,  0.878,  0.868,  1.088),
            ),
    'L1TK' : cms.PSet(
			etaBins = cms.vdouble( 0.500,  1.000,  1.500,  2.000,  2.500),
			offset  = cms.vdouble(-4.484, -4.062, -3.586, -8.047, -1.065),
			scale   = cms.vdouble( 0.546,  0.534,  0.522,  0.659,  0.482),
            ),
    'L1TKV' : cms.PSet(
			etaBins = cms.vdouble( 0.500,  1.000,  1.500,  2.000,  2.500),
			offset  = cms.vdouble(-2.068, -2.033, -1.172, -1.511,  1.154),
			scale   = cms.vdouble( 0.490,  0.484,  0.451,  0.483,  0.398),
            ),
    'L1PF' : cms.PSet(
                        etaBins = cms.vdouble( 0.500,  1.000,  1.500,  2.000,  2.500,  3.000,  3.500,  4.000,  4.500,  5.000),
                        offset  = cms.vdouble( 12.622,  15.549,  18.607,  26.369,  41.212,  38.112,  48.410,  60.027,  45.062,  53.013),
                        scale   = cms.vdouble( 0.988,  0.967,  0.928,  0.926,  0.986,  1.025,  0.748,  0.765,  1.213,  1.276),
            ),
    'L1Puppi' : cms.PSet(
			etaBins = cms.vdouble( 0.500,  1.000,  1.500,  2.000,  2.500,  3.000,  3.500,  4.000,  4.500,  5.000),
			offset  = cms.vdouble(-13.687, -13.801, -10.823, -11.606, -12.974,  4.456,  11.672,  5.632,  15.133, -6.580),
			scale   = cms.vdouble( 0.797,  0.792,  0.706,  0.739,  0.754,  0.131,  0.050,  0.182, -0.004,  0.619),
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
        #L1Calo = cms.InputTag("ak4L1ICalo"),
        #L1TK = cms.InputTag("ak4L1ITK"),
        #L1PF = cms.InputTag("ak4L1IPF"),
        #L1Puppi = cms.InputTag("ak4L1IPuppi"),
        L1Calo = cms.InputTag("ak4L1Calo"),
        L1TK = cms.InputTag("ak4L1TK"),
        L1PF = cms.InputTag("ak4L1PF"),
        L1Puppi = cms.InputTag("ak4L1Puppi"),
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

process.p = cms.Path(process.l1tPFHGCalProducerFrom3DTPsEM + process.CaloInfoOut + process.InfoOut + process.jets + process.ntuple)
process.TFileService = cms.Service("TFileService", fileName = cms.string("jetTupleNew.root"))

 


import FWCore.ParameterSet.Config as cms

process = cms.Process("RESP")

process.load('Configuration.StandardSequences.Services_cff')
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True), allowUnscheduled = cms.untracked.bool(False) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(200))
process.MessageLogger.cerr.FwkReport.reportEvery = 10

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:/eos/cms/store/cmst3/user/gpetrucc/l1phase2/101X/NewInputs/080818/TTbar_PU200/inputs_TTbar_PU200_job1.root'),
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck")
)

process.load('Configuration.Geometry.GeometryExtended2023D17Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff') # needed to read HCal TPs
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '100X_upgrade2023_realistic_v1', '')

process.load("L1Trigger.Phase2L1ParticleFlow.l1ParticleFlow_cff")
process.l1ParticleFlow.remove(process.l1EGammaCrystalsProducer)
process.l1pfProducerTightTK = process.l1pfProducer.clone(trkMinStubs = 6)

process.caloStage2 = cms.EDProducer("CandProducerFromStage2",
    srcCluster = cms.InputTag("simCaloStage2Digis","MP"),
    srcTower = cms.InputTag("simCaloStage2Digis","MP"),
    srcJet = cms.InputTag("simCaloStage2Digis","MP"),
    MP = cms.bool(True),
)


from RecoMET.METProducers.PFMET_cfi import pfMet
pfMet.calculateSignificance = False
process.l1MetCalo    = pfMet.clone(src = "l1pfProducer:Calo")
process.l1MetTK      = pfMet.clone(src = "l1pfProducer:TK")
process.l1MetTKV     = pfMet.clone(src = "l1pfProducer:TKVtx")
process.l1MetTightTK      = pfMet.clone(src = "l1pfProducerTightTK:TK")
process.l1MetTightTKV     = pfMet.clone(src = "l1pfProducerTightTK:TKVtx")
process.l1MetPF      = pfMet.clone(src = "l1pfProducer:PF")
process.l1MetPuppi   = pfMet.clone(src = "l1pfProducer:Puppi")

process.mets = cms.Sequence( process.l1MetCalo + process.l1MetTK + process.l1MetTKV + process.l1MetPF + process.l1MetPuppi + process.l1MetTightTK + process.l1MetTightTKV)

from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets
process.ak4L1Calo    = ak4PFJets.clone(src = 'l1pfProducer:Calo')
process.ak4L1TK      = ak4PFJets.clone(src = 'l1pfProducer:TK')
process.ak4L1TKV     = ak4PFJets.clone(src = 'l1pfProducer:TKVtx')
process.ak4L1TightTK      = ak4PFJets.clone(src = 'l1pfProducerTightTK:TK')
process.ak4L1TightTKV     = ak4PFJets.clone(src = 'l1pfProducerTightTK:TKVtx')
process.ak4L1PF      = ak4PFJets.clone(src = 'l1pfProducer:PF')
process.ak4L1Puppi   = ak4PFJets.clone(src = 'l1pfProducer:Puppi')

process.jets = cms.Sequence(     
    process.ak4L1Calo + process.ak4L1TK + process.ak4L1TKV + process.ak4L1PF + process.ak4L1Puppi  + process.ak4L1TightTK + process.ak4L1TightTKV
)

JEC_PU200 = dict(
### You can create this below with
# for X in L1{Calo,TK,TKV,TightTK,TightTKV,PF,Puppi} Stage2CaloJets RefL1TkCaloJets RefL1TrackerJets; do python scripts/respCorrSimple.py respTupleNew_TTbar_PU200.root plots_dir --fitrange 30 250 -p jet -w ${X}_pt -e ${X}_pt | grep -v [Pp]lot | sed -e "s/simpleCorrEm =/    $X = /" -e "s/)$/    ),"/; done
    L1Calo = cms.PSet(
                        etaBins = cms.vdouble( 0.500,  1.000,  1.500,  2.000,  2.500,  3.000,  3.500,  4.000,  4.500,  5.000),
                        offset  = cms.vdouble( 9.233,  12.095,  12.911,  36.007,  63.148,  84.710,  63.820,  69.563,  60.367,  67.407),
                        scale   = cms.vdouble( 1.025,  1.026,  1.073,  0.956,  1.022,  1.061,  0.897,  0.941,  1.092,  1.324),
    ),
    L1TK =  cms.PSet(
                        etaBins = cms.vdouble( 0.500,  1.000,  1.500,  2.000,  2.500,  3.000,  3.500,  4.000,  4.500,  5.000),
                        offset  = cms.vdouble( 1.605,  0.868,  0.838, -0.293, -1.913,  3.362,  0.000,  0.000,  0.000,  0.000),
                        scale   = cms.vdouble( 0.624,  0.617,  0.582,  0.603,  0.582,  0.010,  1.000,  1.000,  1.000,  1.000),
    ),
    L1TKV =  cms.PSet(
                        etaBins = cms.vdouble( 0.500,  1.000,  1.500,  2.000,  2.500,  3.000,  3.500,  4.000,  4.500,  5.000),
                        offset  = cms.vdouble(-3.009, -3.711, -3.537, -3.644, -5.207,  3.173,  0.000,  0.000,  0.000,  0.000),
                        scale   = cms.vdouble( 0.553,  0.557,  0.538,  0.556,  0.548,  0.012,  1.000,  1.000,  1.000,  1.000),
    ),
    L1TightTK =  cms.PSet(
                        etaBins = cms.vdouble( 0.500,  1.000,  1.500,  2.000,  2.500,  3.000,  3.500,  4.000,  4.500,  5.000),
                        offset  = cms.vdouble(-0.793, -0.516, -0.415, -1.505, -1.413,  2.669,  0.000,  0.000,  0.000,  0.000),
                        scale   = cms.vdouble( 0.345,  0.292,  0.305,  0.426,  0.314,  0.005,  1.000,  1.000,  1.000,  1.000),
    ),
    L1TightTKV =  cms.PSet(
                        etaBins = cms.vdouble( 0.500,  1.000,  1.500,  2.000,  2.500,  3.000,  3.500,  4.000,  4.500,  5.000),
                        offset  = cms.vdouble(-1.676, -1.209, -1.579, -3.231, -2.825,  3.076,  0.000,  0.000,  0.000,  0.000),
                        scale   = cms.vdouble( 0.335,  0.290,  0.303,  0.417,  0.320,  0.005,  1.000,  1.000,  1.000,  1.000),
    ),
    L1PF =  cms.PSet(
                        etaBins = cms.vdouble( 0.500,  1.000,  1.500,  2.000,  2.500,  3.000,  3.500,  4.000,  4.500,  5.000),
                        offset  = cms.vdouble( 12.454,  14.691,  15.236,  40.116,  73.101,  86.003,  63.820,  69.563,  60.367,  67.407),
                        scale   = cms.vdouble( 1.126,  1.126,  1.174,  1.101,  1.131,  1.090,  0.897,  0.941,  1.092,  1.324),
    ),
    L1Puppi =  cms.PSet(
                        etaBins = cms.vdouble( 0.500,  1.000,  1.500,  2.000,  2.500,  3.000,  3.500,  4.000,  4.500,  5.000),
                        offset  = cms.vdouble(-9.922, -10.116, -12.620, -9.839, -3.331,  33.446,  10.677,  10.830,  5.794, -33.316),
                        scale   = cms.vdouble( 1.123,  1.143,  1.175,  1.097,  1.199,  0.672,  0.953,  1.153,  1.343,  1.812),
    ),
    Stage2CaloJets =  cms.PSet(
                        etaBins = cms.vdouble( 0.500,  1.000,  1.500,  2.000,  2.500,  3.000,  3.500,  4.000,  4.500,  5.000),
                        offset  = cms.vdouble( 51.627,  62.332,  66.170,  95.151,  0.000,  76.440,  54.312,  56.799,  58.494,  66.112),
                        scale   = cms.vdouble( 1.164,  1.362,  1.750,  0.165,  1.000,  0.099,  0.834,  0.999,  0.983,  0.711),
    ),
    RefL1TkCaloJets =  cms.PSet(
                        etaBins = cms.vdouble( 0.500,  1.000,  1.500,  2.000,  2.500,  3.000,  3.500,  4.000,  4.500,  5.000),
                        offset  = cms.vdouble( 51.624,  62.463,  66.137,  95.174,  0.000,  72.103,  54.872,  57.316,  58.744,  66.112),
                        scale   = cms.vdouble( 1.163,  1.359,  1.750,  0.164,  1.000,  0.164,  0.826,  0.993,  0.980,  0.711),
    ),
    RefL1TrackerJets =  cms.PSet(
                        etaBins = cms.vdouble( 0.500,  1.000,  1.500,  2.000,  2.500,  3.000,  3.500,  4.000,  4.500,  5.000),
                        offset  = cms.vdouble(-2.067, -2.448, -2.396, -2.586, -2.195,  3.286,  0.000,  0.000,  0.000,  0.000),
                        scale   = cms.vdouble( 0.533,  0.535,  0.506,  0.466,  0.360,  0.012,  1.000,  1.000,  1.000,  1.000),
    ),
)


JEC = JEC_PU200;

process.ntuple = cms.EDAnalyzer("JetMetNTuplizer",
    jets = cms.PSet(
        AK4GenJets = cms.InputTag("ak4GenJetsNoNu"),
        AK4Stage2Calo = cms.InputTag("caloStage2:Jet"),
        AK4CaloJets = cms.InputTag("ak4L1Calo"),
        AK4TKJets = cms.InputTag("ak4L1TK"),
        AK4TKVJets = cms.InputTag("ak4L1TKV"),
        AK4TightTKJets = cms.InputTag("ak4L1TightTK"),
        AK4TightTKVJets = cms.InputTag("ak4L1TightTKV"),
        AK4PFJets = cms.InputTag("ak4L1PF"),
        AK4PuppiJets = cms.InputTag("ak4L1Puppi"),
    ),
    jecs = cms.PSet(
        AK4Stage2CaloJets = JEC['Stage2CaloJets'],
        AK4CaloJets = JEC['L1Calo'],
        AK4TKJets = JEC['L1TK'],
        AK4TKVJets = JEC['L1TKV'],
        AK4TightTKJets = JEC['L1TightTK'],
        AK4TightTKVJets = JEC['L1TightTKV'],
        AK4PFJets = JEC['L1PF'],
        AK4PuppiJets = JEC['L1Puppi'],
    ),
    sels = cms.PSet(
        E13Pt30 = cms.string("pt > 30 && abs(eta) < 1.3"),
        E13Pt40 = cms.string("pt > 40 && abs(eta) < 1.3"),
        #E24Pt15 = cms.string("pt > 15 && abs(eta) < 2.4"),
        #E24Pt20 = cms.string("pt > 20 && abs(eta) < 2.4"),
        E24Pt30 = cms.string("pt > 30 && abs(eta) < 2.4"),
        E24Pt40 = cms.string("pt > 40 && abs(eta) < 2.4"),
        #E24Pt50 = cms.string("pt > 50 && abs(eta) < 2.4"),
        #E24Pt60 = cms.string("pt > 60 && abs(eta) < 2.4"),
        #E24Pt80 = cms.string("pt > 80 && abs(eta) < 2.4"),
        #E30Pt15 = cms.string("pt > 15 && abs(eta) < 3.0"),
        #E30Pt20 = cms.string("pt > 20 && abs(eta) < 3.0"),
        E30Pt30 = cms.string("pt > 30 && abs(eta) < 3.0"),
        E30Pt40 = cms.string("pt > 40 && abs(eta) < 3.0"),
        #E30Pt50 = cms.string("pt > 50 && abs(eta) < 3.0"),
        #E30Pt60 = cms.string("pt > 60 && abs(eta) < 3.0"),
        #E30Pt80 = cms.string("pt > 80 && abs(eta) < 3.0"),
        #E47Pt15 = cms.string("pt > 15 && abs(eta) < 4.7"),
        #E47Pt20 = cms.string("pt > 20 && abs(eta) < 4.7"),
        E47Pt30 = cms.string("pt > 30 && abs(eta) < 4.7"),
        E47Pt40 = cms.string("pt > 40 && abs(eta) < 4.7"),
        #E47Pt50 = cms.string("pt > 50 && abs(eta) < 4.7"),
        #E47Pt60 = cms.string("pt > 50 && abs(eta) < 4.7"),
        #E47Pt80 = cms.string("pt > 80 && abs(eta) < 4.7"),
    ),
    mets = cms.PSet(
        METGen = cms.InputTag("genMetTrue"),
        METCalo = cms.InputTag("l1MetCalo"),
        METTK = cms.InputTag("l1MetTK"),
        METTKV = cms.InputTag("l1MetTKV"),
        METTightTK = cms.InputTag("l1MetTightTK"),
        METTightTKV = cms.InputTag("l1MetTightTKV"),
        METPF = cms.InputTag("l1MetPF"),
        METPuppi = cms.InputTag("l1MetPuppi"),
    ),
    specials = cms.PSet(
        TP_TkEtMiss = cms.PSet( 
            src = cms.InputTag("L1TrackerEtMiss","MET"),
            cut = cms.string(""),
            expr = cms.string("pt")),
        TP_TkHTMissVtx = cms.PSet( 
            src = cms.InputTag("L1TkCaloHTMissVtx","L1TkCaloHTMiss"),
            cut = cms.string(""),
            expr = cms.string("et")),
        TP_TkCaloHTVtx = cms.PSet( 
            src = cms.InputTag("L1TkCaloHTMissVtx","L1TkCaloHTMiss"),
            cut = cms.string(""),
            expr = cms.string("EtTotal")),
        TP_TrackerHTVtx = cms.PSet( 
            src = cms.InputTag("L1TrackerHTMiss","L1TrackerHTMiss"),
            cut = cms.string(""),
            expr = cms.string("EtTotal")),
    )
)


process.p = cms.Path(
        process.caloStage2 +
        process.l1ParticleFlow + process.l1pfProducerTightTK +
        process.mets + process.jets +
        process.ntuple)
process.TFileService = cms.Service("TFileService", fileName = cms.string("jetmetTuple.root"))

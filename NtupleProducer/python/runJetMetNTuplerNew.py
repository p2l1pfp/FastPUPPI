import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras

process = cms.Process("RESP", eras.Phase2C4_trigger)

process.load('Configuration.StandardSequences.Services_cff')
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True), allowUnscheduled = cms.untracked.bool(False) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(200))
process.MessageLogger.cerr.FwkReport.reportEvery = 10

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:/eos/cms/store/cmst3/group/hzz/gpetrucc/tmp/prod104X/TTbar_14TeV_TuneCP5_Pythia8_PU0/TTbar_14TeV_TuneCP5_Pythia8_PU0.batch1.job40.root'),
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck")
)

process.load('Configuration.Geometry.GeometryExtended2023D35Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2023D35_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff') # needed to read HCal TPs
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('RecoMET.Configuration.GenMETParticles_cff')
process.load('RecoMET.METProducers.genMetTrue_cfi')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '103X_upgrade2023_realistic_v2', '') 

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

process.centralPF = cms.EDFilter("CandPtrSelector", src = cms.InputTag("l1pfProducer:PF"), cut = cms.string("abs(eta) < 2.4"))
process.centralCalo = cms.EDFilter("CandPtrSelector", src = cms.InputTag("l1pfProducer:Calo"), cut = cms.string("abs(eta) < 2.4"))
process.centralPuppi = cms.EDFilter("CandPtrSelector", src = cms.InputTag("l1pfProducer:Puppi"), cut = cms.string("abs(eta) < 2.4"))
process.centralGen = cms.EDFilter("CandPtrSelector", src = cms.InputTag("genParticlesForMETAllVisible"), cut = cms.string("abs(eta) < 2.4"))
process.genMetCentralTrue = process.genMetTrue.clone(src = cms.string("centralGen"))
process.l1MetCentralCalo  = pfMet.clone(src = "centralCalo")
process.l1MetCentralPF    = pfMet.clone(src = "centralPF")
process.l1MetCentralPuppi = pfMet.clone(src = "centralPuppi")


process.mets = cms.Sequence( process.l1MetCalo + process.l1MetTK + process.l1MetTKV + process.l1MetPF + process.l1MetPuppi + process.l1MetTightTK + process.l1MetTightTKV +
        process.genParticlesForMETAllVisible +
        process.centralGen + 
        process.centralCalo + process.centralPF + process.centralPuppi +
        process.genMetCentralTrue +
        process.l1MetCentralCalo + process.l1MetCentralPF + process.l1MetCentralPuppi 
)

        

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
    L1Calo =  cms.PSet(
                        etaBins = cms.vdouble( 0.500,  1.000,  1.500,  2.000,  2.500,  3.000,  3.500,  4.000,  4.500,  5.000),
                        offset  = cms.vdouble( 14.041,  18.421,  27.169,  70.760,  100.682,  99.327,  67.345,  72.802,  67.144,  56.801),
                        scale   = cms.vdouble( 0.913,  0.909,  0.933,  0.809,  0.792,  0.808,  0.868,  0.862,  0.960,  1.542),
    ),
    L1PF =  cms.PSet(
                        etaBins = cms.vdouble( 0.500,  1.000,  1.500,  2.000,  2.500,  3.000,  3.500,  4.000,  4.500,  5.000),
                        offset  = cms.vdouble( 17.467,  22.076,  31.229,  78.241,  115.035,  102.950,  67.345,  72.802,  67.144,  56.801),
                        scale   = cms.vdouble( 1.041,  1.030,  1.064,  0.997,  0.956,  0.842,  0.868,  0.862,  0.960,  1.542),
    ),
    L1Puppi =  cms.PSet(
                        etaBins = cms.vdouble( 0.500,  1.000,  1.500,  2.000,  2.500,  3.000,  3.500,  4.000,  4.500,  5.000),
                        offset  = cms.vdouble(-6.664, -5.996, -6.798,  3.459,  13.866,  55.320,  10.236,  11.906,  7.550, -43.093),
                        scale   = cms.vdouble( 1.052,  1.065,  1.121,  1.100,  1.157,  0.264,  0.978,  1.100,  1.276,  2.046),
    ),
    L1TKV =  cms.PSet(
                        etaBins = cms.vdouble( 0.500,  1.000,  1.500,  2.000,  2.500,  3.000,  3.500,  4.000,  4.500,  5.000),
                        offset  = cms.vdouble(-3.428, -3.409, -3.242, -3.440, -4.713,  3.134,  0.000,  0.000,  0.000,  0.000),
                        scale   = cms.vdouble( 0.560,  0.552,  0.533,  0.550,  0.541,  0.013,  1.000,  1.000,  1.000,  1.000),
    ),
    L1TightTKV =  cms.PSet(
                        etaBins = cms.vdouble( 0.500,  1.000,  1.500,  2.000,  2.500,  3.000,  3.500,  4.000,  4.500,  5.000),
                        offset  = cms.vdouble(-2.218, -1.328, -1.397, -3.178, -1.838,  2.940,  0.000,  0.000,  0.000,  0.000),
                        scale   = cms.vdouble( 0.345,  0.292,  0.299,  0.415,  0.305,  0.004,  1.000,  1.000,  1.000,  1.000),
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
        AK4TKVJets = cms.InputTag("ak4L1TKV"),
        AK4TightTKVJets = cms.InputTag("ak4L1TightTKV"),
        AK4PFJets = cms.InputTag("ak4L1PF"),
        AK4PuppiJets = cms.InputTag("ak4L1Puppi"),
    ),
    jecs = cms.PSet(
        AK4Stage2CaloJets = JEC['Stage2CaloJets'],
        AK4CaloJets = JEC['L1Calo'],
        AK4TKVJets = JEC['L1TKV'],
        AK4TightTKVJets = JEC['L1TightTKV'],
        AK4PFJets = JEC['L1PF'],
        AK4PuppiJets = JEC['L1Puppi'],
    ),
    sels = cms.PSet(
        E13Pt30 = cms.string("pt > 30 && abs(eta) < 1.3"),
        E13Pt40 = cms.string("pt > 40 && abs(eta) < 1.3"),
        #E24Pt15 = cms.string("pt > 15 && abs(eta) < 2.4"),
        #E24Pt20 = cms.string("pt > 20 && abs(eta) < 2.4"),
        E21Pt30 = cms.string("pt > 30 && abs(eta) < 2.1"),
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
        METGenCentral = cms.InputTag("genMetCentralTrue"),
        METCaloCentral = cms.InputTag("l1MetCentralCalo"),
        METPFCentral = cms.InputTag("l1MetCentralPF"),
        METPuppiCentral = cms.InputTag("l1MetCentralPuppi"),
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

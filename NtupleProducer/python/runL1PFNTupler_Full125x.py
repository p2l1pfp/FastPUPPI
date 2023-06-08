import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras

process = cms.Process("IN", eras.Phase2C17I13M9)
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.Geometry.GeometryExtended2026D88Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2026D88_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '125X_mcRun4_realistic_v2', '')

process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff')
process.load('CalibCalorimetry.CaloTPG.CaloTPGTranscoder_cfi')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('L1Trigger.TrackTrigger.TrackTrigger_cff')
process.load("L1Trigger.TrackFindingTracklet.L1HybridEmulationTracks_cff") 
process.load("L1Trigger.TrackerDTC.ProducerES_cff") 
process.load("L1Trigger.TrackerDTC.ProducerED_cff") 
process.load("RecoVertex.BeamSpotProducer.BeamSpot_cfi")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
         #'root://cms-xrd-global.cern.ch//store/mc/Phase2HLTTDRSummer20ReRECOMiniAOD/TT_TuneCP5_14TeV-powheg-pythia8/FEVT/PU200_111X_mcRun4_realistic_T15_v1-v2/280000/003ACFBC-23B2-EA45-9A12-BECFF07760FC.root',
         'file:/afs/cern.ch/work/d/ddiaz/L1SampleGen/L1_sampleGen/12_5_x/GEN-SIMDIGIRAW.root',
    )
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1))
process.options = cms.untracked.PSet( 
        wantSummary = cms.untracked.bool(True),
        #numberOfThreads = cms.untracked.uint32(4),
        #numberOfStreams = cms.untracked.uint32(4),
)

process.load('L1Trigger.Phase2L1GMT.gmt_cfi')

process.PFInputsTask = cms.Task(
    process.TTClustersFromPhase2TrackerDigis,
    process.TTStubsFromPhase2TrackerDigis,
    process.TrackerDTCProducer,
    process.offlineBeamSpot,
    process.l1tTTTracksFromTrackletEmulation,
    process.l1tTTTracksFromExtendedTrackletEmulation,
    process.TTTrackAssociatorFromPixelDigis,
    process.TTTrackAssociatorFromPixelDigisExtended,
    process.SimL1EmulatorTask
)

process.p = cms.Path(
    process.TrackTriggerClustersStubs +
    process.offlineBeamSpot +
    process.l1tTTTracksFromExtendedTrackletEmulation +
    process.SimL1Emulator
)
process.p.associate(process.PFInputsTask)
# next two lines not needed anynmore
# updated code does this, just take the proper collections
# for example l1tLayer1 --> l1tLayer1Extended
##process.l1tTTTracksFromTrackletEmulation.Hnpar = 5
##process.l1tPFTracksFromL1Tracks.nParam = 5

###process.L1TkPrimaryVertex.nVtx = cms.int32(5)
###process.Vertex.nVtx = cms.int32(5)

process.ntuple0 = cms.EDAnalyzer("L1PFCompare",
    emcalo = cms.InputTag(""),
    egcalo = cms.InputTag(""),
    calo = cms.InputTag(""),
    pf = cms.InputTag("l1tLayer1:PF"),
    pup = cms.InputTag("l1tLayer1:Puppi"),
    vtx = cms.InputTag("l1tVertexProducer", "l1vertices"),
    generator = cms.InputTag('genParticles'),
    l1jet = cms.InputTag(""),
    l1bName = cms.InputTag(""),
    recojet = cms.InputTag(""),
    recobName = cms.string(""),
    minPt = cms.double(2.),
    maxEta = cms.double(3.),
    maxN = cms.uint32(9999),
    genIDs = cms.vint32(1000006,-1000006,5,-5,-4,4),
    addGenIDs = cms.vint32(),
    genStatuses = cms.vint32(22,22,23,23,23,23),
)
process.p += process.ntuple0

process.schedule = cms.Schedule([process.p])

process.TFileService = cms.Service("TFileService", fileName = cms.string("pfTuple.root"))

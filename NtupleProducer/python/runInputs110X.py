import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras

process = cms.Process("IN", eras.Phase2C9)
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.Geometry.GeometryExtended2026D49Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2026D49_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '123X_mcRun4_realistic_v3', '')

process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff')
process.load('CalibCalorimetry.CaloTPG.CaloTPGTranscoder_cfi')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('L1Trigger.TrackTrigger.TrackTrigger_cff')
process.load("L1Trigger.TrackFindingTracklet.L1HybridEmulationTracks_cff") 
process.load("L1Trigger.TrackerDTC.ProducerES_cff") 
process.load("L1Trigger.TrackerDTC.ProducerED_cff") 
process.load("RecoVertex.BeamSpotProducer.BeamSpot_cfi")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('/store/mc/Phase2HLTTDRWinter20DIGI/TT_TuneCP5_14TeV-powheg-pythia8/GEN-SIM-DIGI-RAW/PU200_110X_mcRun4_realistic_v3-v2/110000/005E74D6-B50E-674E-89E6-EAA9A617B476.root',)
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(5))
process.options = cms.untracked.PSet( 
        wantSummary = cms.untracked.bool(True),
        numberOfThreads = cms.untracked.uint32(4),
        numberOfStreams = cms.untracked.uint32(4),
)

process.p = cms.Path(
    process.TrackTriggerClustersStubs +
    process.TrackerDTCProducer +
    process.offlineBeamSpot +
    process.TTTracksFromTrackletEmulation
    + process.TrackTriggerAssociatorTracks
    + process.TTTracksFromExtendedTrackletEmulation
    + process.TTTrackAssociatorFromPixelDigisExtended
)
process.p.associate(process.SimL1EmulatorTask)

process.out = cms.OutputModule("PoolOutputModule",
        fileName = cms.untracked.string("inputs110X.root"),
        outputCommands = cms.untracked.vstring("drop *",
            # --- GEN
            "keep *_genParticles_*_*",
            "keep *_ak4GenJetsNoNu_*_*",
            "keep *_genMetTrue_*_*",
            # --- Track TPs
            "keep *_TTTracksFromTrackletEmulation_*_*",
            "keep *_TrackTriggerAssociatorTracks_*_*",
            "keep *_TTTracksFromExtendedTrackletEmulation_*_*",
            "keep *_TTTrackAssociatorFromPixelDigisExtended_*_*",
            # --- Calo TPs
            "keep *_simHcalTriggerPrimitiveDigis_*_*",
            "keep *_simCaloStage2Layer1Digis_*_*",
            "keep *_simCaloStage2Digis_*_*",
            # --- Muon TPs
            "keep *_simMuonRPCDigis_*_*",
            "keep *_simMuonGEMPadDigis_*_*",
            "keep *_simMuonGEMPadDigiClusters_*_*",
            "keep *_simDtTriggerPrimitiveDigis_*_*",
            "keep *_simCscTriggerPrimitiveDigis_*_*",
            "keep *_simTwinMuxDigis_*_*",
            "keep *_simBmtfDigis_*_*",
            "keep *_simKBmtfStubs_*_*",
            "keep *_simKBmtfDigis_*_*",
            "keep *_simEmtfDigis_*_*",
            "keep *_simOmtfDigis_*_*",
            "keep *_simGmtCaloSumDigis_*_*",
            "keep *_simGmtStage2Digis_*_*",
            "keep *_simEmtfShowers_*_*",
            "keep *_simGmtShowerDigis_*_*",
            "keep *_simCscTriggerPrimitiveDigisRun3_*_*",
            "keep *_simMuonME0PadDigis_*_*",
            "keep *_me0TriggerDigis_*_*",
            "keep *_simMuonME0PseudoReDigisCoarse_*_*",
            "keep *_me0RecHitsCoarse_*_*",
            "keep *_me0TriggerPseudoDigis_*_*",
            "keep *_me0RecHits_*_*",
            "keep *_me0Segments_*_*",
            "keep *_me0TriggerConvertedPseudoDigis_*_*",
            "keep *_simCscTriggerPrimitiveDigisPhase2_*_*",
            "keep *_simGtExtFakeStage2Digis_*_*",
            "keep *_simGtStage2Digis_*_*",
            "keep *_CalibratedDigis_*_*",
            "keep *_dtTriggerPhase2PrimitiveDigis_*_*",
            # --- HGCal TPs
            #"keep *_hgcalVFEProducer_*_*",
            "keep l1tHGCalTriggerCellBXVector_hgcalConcentratorProducer_*_*",
            "keep *_hgcalBackEndLayer2Producer_HGCalBackendLayer2Processor3DClustering_*",
            "keep *_hgcalTowerProducer_*_*",
            # --- GCT reconstruction
            "keep *_L1EGammaClusterEmuProducer_*_*",
            "keep *_l1EGammaEEProducer_*_*",
            "keep *_L1TowerCalibration_*_*",
            "keep *_L1CaloJet_*_*",
            "keep *_L1CaloJetHTT_*_*",
            # --- GTT reconstruction
            "keep *_L1VertexFinder_*_*",
            #"keep *_L1GTTInputProducer_*_*",
            #"keep *_L1GTTInputProducerExtended_*_*",
            "keep *_L1VertexFinderEmulator_*_*",
            "keep *_L1TkPrimaryVertex_*_*",
            #"keep *_L1TrackJets_*_*", # can't be saved: The class "l1t::TkJetWord" is compiled and for its data member "tkJetWord_" we do not have a dictionary for the collection "bitset<70>". Because of this, we will not be able to read or write this data member.
            #"keep *_L1TrackJetsExtended_*_*",
            #"keep *_L1TrackFastJets_*_*",
            "keep *_L1TrackerEtMiss_*_*",
            "keep *_L1TrackerHTMiss_*_*",
            #"keep *_L1TrackJetsEmulation_*_*",
            #"keep *_L1TrackJetsExtendedEmulation_*_*",
            "keep *_L1TrackerEmuEtMiss_*_*",
            "keep *_L1TrackerEmuHTMiss_*_*",
            "keep *_L1TrackerEmuHTMissExtended_*_*",
            # --- GMT reconstruction
            "keep *_L1TkMuons_*_*",
            "keep *_L1TkMuonsTP_*_*",
            "keep *_L1TkGlbMuons_*_*",
            "keep *_L1TkStubsGmt_*_*",
            "keep *_L1TkMuonsGmt_*_*",
            "keep *_L1SAMuonsGmt_*_*",
            # --- TDR (obsolete) Tk-Egamma reconstruction
            "keep *_L1TkElectronsCrystal_*_*",
            "keep *_L1TkElectronsLooseCrystal_*_*",
            "keep *_L1TkElectronsEllipticMatchCrystal_*_*",
            "keep *_L1TkIsoElectronsCrystal_*_*",
            "keep *_L1TkPhotonsCrystal_*_*",
            "keep *_L1TkElectronsHGC_*_*",
            "keep *_L1TkElectronsEllipticMatchHGC_*_*",
            "keep *_L1TkElectronsLooseHGC_*_*",
            "keep *_L1TkIsoElectronsHGC_*_*",
            "keep *_L1TkPhotonsHGC_*_*",
        ),
        compressionAlgorithm = cms.untracked.string('LZMA'),
        compressionLevel = cms.untracked.int32(4),
        dropMetaData = cms.untracked.string('ALL'),
        fastCloning = cms.untracked.bool(False),
        overrideInputFileSplitLevels = cms.untracked.bool(True),
        eventAutoFlushCompressedSize = cms.untracked.int32(15728640),
        SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring("p")),
)
process.e = cms.EndPath(process.out)

process.schedule = cms.Schedule([process.p,process.e])

def goSlim():
    process.out.outputCommands += [ "drop *_hgcalConcentratorProducer*_*_*", "drop *_hgcalTowerMapProducer_*_*", ]

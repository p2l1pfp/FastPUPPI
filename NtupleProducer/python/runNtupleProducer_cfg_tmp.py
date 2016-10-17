import os
import FWCore.ParameterSet.Config as cms

process = cms.Process("OUT")

# process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                            # replace 'myfile.root' with the source file you want to use
                            fileNames = cms.untracked.vstring(
        'XXXX'
        )
                            )


process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedGauss_cfi')
process.load('Configuration.Geometry.GeometryExtended2023TTIReco_cff')
process.load('Geometry.TrackerGeometryBuilder.StackedTrackerGeometry_cfi')

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('L1Trigger.Configuration.SimL1Emulator_cff')
process.load('L1Trigger.Configuration.CaloTriggerPrimitives_cff')
process.load("SLHCUpgradeSimulations.L1TrackTrigger.L1TTrack_cfi")

process.load('Configuration.StandardSequences.L1TrackTrigger_cff')
process.BeamSpotFromSim =cms.EDProducer("BeamSpotFromSimProducer")
process.L1Tracks.geometry = cms.untracked.string('BE5D')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS3', '')
cmssw_base = os.environ['CMSSW_BASE']
process.load('FastPUPPI.NtupleProducer.Puppi_cff')
process.InfoOut = cms.EDProducer('NtupleProducer',
                                 zeroSuppress = cms.bool(True),
                                 L1TrackTag  = cms.InputTag('TTTracksFromPixelDigis','Level1TTTracks'),
                                 EcalTPTag   = cms.InputTag('simEcalTriggerPrimitiveDigis'),
                                 HcalTPTag   = cms.InputTag('simHcalTriggerPrimitiveDigis'), 
                                 MuonTPTag   = cms.InputTag('simGmtDigis'), 
                                 genParTag   = cms.InputTag('genParticles'), 
                                 corrector   = cms.InputTag(cmssw_base+"/src/FastPUPPI/NtupleProducer/data/pion_eta_phi.root"),
                                 corrector2  = cms.InputTag(cmssw_base+"/src/FastPUPPI/NtupleProducer/data/pion_eta_phi_hpu.root"),
                                 ecorrector  = cms.InputTag(cmssw_base+"/src/FastPUPPI/NtupleProducer/data/ecorr.root"),
                                 trackres    = cms.InputTag(cmssw_base+"/src/FastPUPPI/NtupleProducer/data/tkres.root"),
                                 eleres      = cms.InputTag(cmssw_base+"/src/FastPUPPI/NtupleProducer/data/eres.root"),
                                 pionres     = cms.InputTag(cmssw_base+"/src/FastPUPPI/NtupleProducer/data/pionres.root"),
                                 puppi       = process.puppi
                                 )

process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string('TESTFILE.root'),
                               outputCommands = cms.untracked.vstring('drop *',
                                                                      "keep *_simEcalTriggerPrimitiveDigis_*_*",
                                                                      "keep *_*_ecalET_*",
                                                                      "keep *_*_trkPtPerp_*",
                                                                      "keep *_*_trkPtEta_*",
                                                                      "keep *_*_trkPtPhi_*",
                                                                      "keep *_*_trkPOCAz_*",)
                               )


process.TT_step           = cms.Path(process.TrackTriggerTTTracks)
process.TTAssociator_step = cms.Path(process.TrackTriggerAssociatorTracks)
process.ana               = cms.Path(process.InfoOut)
process.p = cms.Schedule(process.TT_step,process.TTAssociator_step,process.ana)
#process.e = cms.EndPath(process.out)

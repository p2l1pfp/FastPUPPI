import FWCore.ParameterSet.Config as cms

process = cms.Process("OUT")

# process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.source = cms.Source("PoolSource",
                            # replace 'myfile.root' with the source file you want to use
                            fileNames = cms.untracked.vstring(
        '/store/mc/TTI2023Upg14D/Zmumu_TuneZ2star_14TeV_Eta4-pythia6/GEN-SIM-DIGI-RAW/PU140bx25_PH2_1K_FB_V3-v2/00000/08B6D2DB-D5E6-E311-9356-002618943975.root'
        )
                            )

# process.load('Configuration/StandardSequences/GeometryIdeal_cff')
# process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.Geometry.GeometryExtended2023TTIReco_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.L1TrackTrigger_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.load('L1Trigger.Configuration.SimL1Emulator_cff')
process.load('L1Trigger.Configuration.CaloTriggerPrimitives_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS3', '')
# process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgrade2019', '')
process.InfoOut = cms.EDProducer('NtupleProducer',
                                 L1TrackTag  = cms.InputTag('TTTracksFromPixelDigis','Level1TTTracks'),
                                 EcalTPTag   = cms.InputTag('simEcalTriggerPrimitiveDigis'),
                                 HcalTPTag   = cms.InputTag('simHcalTriggerPrimitiveDigis'), 
                                 genParTag   = cms.InputTag('genParticles'), 
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


process.p = cms.Path(process.InfoOut)

process.e = cms.EndPath(process.out)

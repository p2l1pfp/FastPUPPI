import FWCore.ParameterSet.Config as cms

process = cms.Process("OUT")

# process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(500) )

process.source = cms.Source("PoolSource",
                            # replace 'myfile.root' with the source file you want to use
                            fileNames = cms.untracked.vstring(
        '/store/mc/TTI2023Upg14D/Zmumu_TuneZ2star_14TeV_Eta4-pythia6/GEN-SIM-DIGI-RAW/PU140bx25_PH2_1K_FB_V3-v2/00000/F6C61A52-E6E6-E311-B866-002618FDA259.root'
        #'/store/mc/TTI2023Upg14D/SinglePionPlusFlatPt0p2To50/GEN-SIM-DIGI-RAW/PU140bx25_PH2_1K_FB_V3-v2/00000/0065654D-E9E6-E311-A6E1-0026189438B0.root',
        #'file:/tmp/pharris/0065654D-E9E6-E311-A6E1-0026189438B0.root',
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
                                 zeroSuppress = cms.bool(True),
                                 L1TrackTag  = cms.InputTag('TTTracksFromPixelDigis','Level1TTTracks'),
                                 EcalTPTag   = cms.InputTag('simEcalTriggerPrimitiveDigis'),
                                 HcalTPTag   = cms.InputTag('simHcalTriggerPrimitiveDigis'), 
                                 genParTag   = cms.InputTag('genParticles'),
                                 corrector   = cms.InputTag("/afs/cern.ch/user/p/pharris/pharris/public/bacon/prod/CMSSW_6_2_0_SLHC12_patch1/src/FastPUPPI/NtupleProducer/data/corr.root"),
                                 ecorrector  = cms.InputTag("/afs/cern.ch/user/p/pharris/pharris/public/bacon/prod/CMSSW_6_2_0_SLHC12_patch1/src/FastPUPPI/NtupleProducer/data/ecorr.root"),
                                 trackres    = cms.InputTag("/afs/cern.ch/user/p/pharris/pharris/public/bacon/prod/CMSSW_6_2_0_SLHC12_patch1/src/FastPUPPI/NtupleProducer/data/tkres.root"),
                                 eleres      = cms.InputTag("/afs/cern.ch/user/p/pharris/pharris/public/bacon/prod/CMSSW_6_2_0_SLHC12_patch1/src/FastPUPPI/NtupleProducer/data/eres.root"),
                                 pionres     = cms.InputTag("/afs/cern.ch/user/p/pharris/pharris/public/bacon/prod/CMSSW_6_2_0_SLHC12_patch1/src/FastPUPPI/NtupleProducer/data/pionres.root")
                                 )


process.load('RecoMET.METProducers.PFMET_cfi')
process.pfMet.src = cms.InputTag('InfoOut')
process.pfMet.calculateSignificance = False

process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string('TESTFILE.root'),
                               outputCommands = cms.untracked.vstring('drop *',
                                                                      "keep *_*_*Met_*",
                                                                      "keep *_simEcalTriggerPrimitiveDigis_*_*",
                                                                      "keep *_*_ecalET_*",
                                                                      "keep *_*_trkPtPerp_*",
                                                                      "keep *_*_trkPtEta_*",
                                                                      "keep *_*_trkPtPhi_*",
                                                                      "keep *_*_trkPOCAz_*",
                                                                      "keep *_*_*_OUT")
                               )


process.p = cms.Path(process.InfoOut*process.pfMet)
process.e = cms.EndPath(process.out)

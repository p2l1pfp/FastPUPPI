import FWCore.ParameterSet.Config as cms

process = cms.Process("OUT")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                            # replace 'myfile.root' with the source file you want to use
                            fileNames = cms.untracked.vstring(
        'file:/afs/cern.ch/work/s/ssevova/public/eos/cms/store/user/ssevova/cms/store/mc/TTI2023Upg14D/Zmumu_TuneZ2star_14TeV_Eta4-pythia6/GEN-SIM-DIGI-RAW/PU140bx25_PH2_1K_FB_V3-v2/00000/0CBC405D-E1E6-E311-9F5A-001A92971B72.root'#08B6D2DB-D5E6-E311-9356-002618943975.root'  
        )
                            )

process.InfoOut = cms.EDProducer('NtupleProducer',
                                 L1TrackTag  = cms.InputTag('TTTracksFromPixelDigis','Level1TTTracks'),
                                 EcalTPTag   = cms.InputTag('simEcalTriggerPrimitiveDigis'),
                                 HcalTPTag   = cms.InputTag('simHcalTriggerPrimitiveDigis'), 
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

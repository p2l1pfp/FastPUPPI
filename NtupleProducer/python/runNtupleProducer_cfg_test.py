import FWCore.ParameterSet.Config as cms

process = cms.Process("OUT")

# process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(20) )

process.source = cms.Source("PoolSource",
                            # replace 'myfile.root' with the source file you want to use
                            fileNames = cms.untracked.vstring(
'/store/mc/TTI2023Upg14D/Zmumu_TuneZ2star_14TeV_Eta4-pythia6/GEN-SIM-DIGI-RAW/PU140bx25_PH2_1K_FB_V3-v2/00000/08B6D2DB-D5E6-E311-9356-002618943975.root',
'/store/mc/TTI2023Upg14D/Zmumu_TuneZ2star_14TeV_Eta4-pythia6/GEN-SIM-DIGI-RAW/PU140bx25_PH2_1K_FB_V3-v2/00000/0CBC405D-E1E6-E311-9F5A-001A92971B72.root',
'/store/mc/TTI2023Upg14D/Zmumu_TuneZ2star_14TeV_Eta4-pythia6/GEN-SIM-DIGI-RAW/PU140bx25_PH2_1K_FB_V3-v2/00000/0E451F74-D9E6-E311-82F8-0025905A6134.root',
'/store/mc/TTI2023Upg14D/Zmumu_TuneZ2star_14TeV_Eta4-pythia6/GEN-SIM-DIGI-RAW/PU140bx25_PH2_1K_FB_V3-v2/00000/0E8A6B69-DAE6-E311-8F69-003048FFCB8C.root',
'/store/mc/TTI2023Upg14D/Zmumu_TuneZ2star_14TeV_Eta4-pythia6/GEN-SIM-DIGI-RAW/PU140bx25_PH2_1K_FB_V3-v2/00000/0EB4C094-D5E6-E311-A842-0025905A608A.root',
'/store/mc/TTI2023Upg14D/Zmumu_TuneZ2star_14TeV_Eta4-pythia6/GEN-SIM-DIGI-RAW/PU140bx25_PH2_1K_FB_V3-v2/00000/1E48918B-D8E6-E311-89CA-003048FFD760.root',
'/store/mc/TTI2023Upg14D/Zmumu_TuneZ2star_14TeV_Eta4-pythia6/GEN-SIM-DIGI-RAW/PU140bx25_PH2_1K_FB_V3-v2/00000/1E597BA3-D4E6-E311-B9CB-00261894386E.root',
'/store/mc/TTI2023Upg14D/Zmumu_TuneZ2star_14TeV_Eta4-pythia6/GEN-SIM-DIGI-RAW/PU140bx25_PH2_1K_FB_V3-v2/00000/1EBE196A-D3E6-E311-BE65-002618FDA262.root',
'/store/mc/TTI2023Upg14D/Zmumu_TuneZ2star_14TeV_Eta4-pythia6/GEN-SIM-DIGI-RAW/PU140bx25_PH2_1K_FB_V3-v2/00000/2220B0EB-CDE6-E311-A118-003048678FEA.root',
'/store/mc/TTI2023Upg14D/Zmumu_TuneZ2star_14TeV_Eta4-pythia6/GEN-SIM-DIGI-RAW/PU140bx25_PH2_1K_FB_V3-v2/00000/24B321CD-DBE6-E311-8297-003048678DA2.root',
'/store/mc/TTI2023Upg14D/Zmumu_TuneZ2star_14TeV_Eta4-pythia6/GEN-SIM-DIGI-RAW/PU140bx25_PH2_1K_FB_V3-v2/00000/28BC6679-CEE6-E311-85AF-0025905938A8.root',
'/store/mc/TTI2023Upg14D/Zmumu_TuneZ2star_14TeV_Eta4-pythia6/GEN-SIM-DIGI-RAW/PU140bx25_PH2_1K_FB_V3-v2/00000/2A73BD15-DCE6-E311-A93B-00304867D836.root',
'/store/mc/TTI2023Upg14D/Zmumu_TuneZ2star_14TeV_Eta4-pythia6/GEN-SIM-DIGI-RAW/PU140bx25_PH2_1K_FB_V3-v2/00000/32197DBA-CEE6-E311-A909-00304867BED8.root',
'/store/mc/TTI2023Upg14D/Zmumu_TuneZ2star_14TeV_Eta4-pythia6/GEN-SIM-DIGI-RAW/PU140bx25_PH2_1K_FB_V3-v2/00000/3224F1DA-CAE6-E311-B22E-0030486790A0.root',
'/store/mc/TTI2023Upg14D/Zmumu_TuneZ2star_14TeV_Eta4-pythia6/GEN-SIM-DIGI-RAW/PU140bx25_PH2_1K_FB_V3-v2/00000/422A08A6-D7E6-E311-BB3B-0025905A48D6.root',
'/store/mc/TTI2023Upg14D/Zmumu_TuneZ2star_14TeV_Eta4-pythia6/GEN-SIM-DIGI-RAW/PU140bx25_PH2_1K_FB_V3-v2/00000/44424C0A-D4E6-E311-B6C6-003048678B3C.root',
'/store/mc/TTI2023Upg14D/Zmumu_TuneZ2star_14TeV_Eta4-pythia6/GEN-SIM-DIGI-RAW/PU140bx25_PH2_1K_FB_V3-v2/00000/44542505-DAE6-E311-9129-0025905A612E.root',
'/store/mc/TTI2023Upg14D/Zmumu_TuneZ2star_14TeV_Eta4-pythia6/GEN-SIM-DIGI-RAW/PU140bx25_PH2_1K_FB_V3-v2/00000/46FDB35B-E0E6-E311-9A91-003048678B38.root',
'/store/mc/TTI2023Upg14D/Zmumu_TuneZ2star_14TeV_Eta4-pythia6/GEN-SIM-DIGI-RAW/PU140bx25_PH2_1K_FB_V3-v2/00000/4A430606-E0E6-E311-A8D4-003048678C06.root',
'/store/mc/TTI2023Upg14D/Zmumu_TuneZ2star_14TeV_Eta4-pythia6/GEN-SIM-DIGI-RAW/PU140bx25_PH2_1K_FB_V3-v2/00000/5416165C-CFE6-E311-A209-0025905A6126.root',
'/store/mc/TTI2023Upg14D/Zmumu_TuneZ2star_14TeV_Eta4-pythia6/GEN-SIM-DIGI-RAW/PU140bx25_PH2_1K_FB_V3-v2/00000/5648EF5B-E5E6-E311-B41E-003048678FD6.root',
'/store/mc/TTI2023Upg14D/Zmumu_TuneZ2star_14TeV_Eta4-pythia6/GEN-SIM-DIGI-RAW/PU140bx25_PH2_1K_FB_V3-v2/00000/5A88D387-D2E6-E311-91E4-00261894390E.root',
'/store/mc/TTI2023Upg14D/Zmumu_TuneZ2star_14TeV_Eta4-pythia6/GEN-SIM-DIGI-RAW/PU140bx25_PH2_1K_FB_V3-v2/00000/5AD04CDF-D2E6-E311-B175-003048678E8A.root',
'/store/mc/TTI2023Upg14D/Zmumu_TuneZ2star_14TeV_Eta4-pythia6/GEN-SIM-DIGI-RAW/PU140bx25_PH2_1K_FB_V3-v2/00000/5EFC27E6-D4E6-E311-9EBB-0026189437FA.root',
'/store/mc/TTI2023Upg14D/Zmumu_TuneZ2star_14TeV_Eta4-pythia6/GEN-SIM-DIGI-RAW/PU140bx25_PH2_1K_FB_V3-v2/00000/64F284C7-D6E6-E311-A1DF-001A92810AB2.root',
#        '/store/mc/TTI2023Upg14D/Zmumu_TuneZ2star_14TeV_Eta4-pythia6/GEN-SIM-DIGI-RAW/PU140bx25_PH2_1K_FB_V3-v2/00000/F6C61A52-E6E6-E311-B866-002618FDA259.root'
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
                                 MuonTPTag   = cms.InputTag('simGmtDigis'), 
                                 genParTag   = cms.InputTag('genParticles'),
                                 corrector   = cms.InputTag("/afs/cern.ch/user/p/pharris/pharris/public/bacon/prod/CMSSW_6_2_0_SLHC12_patch1/src/FastPUPPI/NtupleProducer/data/corr.root"),
                                 ecorrector  = cms.InputTag("/afs/cern.ch/user/p/pharris/pharris/public/bacon/prod/CMSSW_6_2_0_SLHC12_patch1/src/FastPUPPI/NtupleProducer/data/ecorr.root"),
                                 trackres    = cms.InputTag("/afs/cern.ch/user/p/pharris/pharris/public/bacon/prod/CMSSW_6_2_0_SLHC12_patch1/src/FastPUPPI/NtupleProducer/data/tkres.root"),
                                 eleres      = cms.InputTag("/afs/cern.ch/user/p/pharris/pharris/public/bacon/prod/CMSSW_6_2_0_SLHC12_patch1/src/FastPUPPI/NtupleProducer/data/eres.root"),
                                 pionres     = cms.InputTag("/afs/cern.ch/user/p/pharris/pharris/public/bacon/prod/CMSSW_6_2_0_SLHC12_patch1/src/FastPUPPI/NtupleProducer/data/pionres.root")
                                 )


process.load('RecoMET.METProducers.PFMET_cfi')
process.pfMet.src = cms.InputTag('InfoOut','PF')
process.pfMet.calculateSignificance = False
process.pfMetCalo = process.pfMet.clone()
process.pfMetCalo.src = cms.InputTag('InfoOut','Calo')
process.pfMetCalo.calculateSignificance = False
process.pfMetRawCalo = process.pfMet.clone()
process.pfMetRawCalo.src = cms.InputTag('InfoOut','RawCalo')
process.pfMetRawCalo.calculateSignificance = False
process.pfMetTK = process.pfMet.clone()
process.pfMetTK.src = cms.InputTag('InfoOut','TK')
process.pfMetTK.calculateSignificance = False

process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string('TESTFILE.root'),
                               outputCommands = cms.untracked.vstring('drop *') # killing my disk space
                               )


process.p = cms.Path(process.InfoOut*process.pfMet*process.pfMetCalo*process.pfMetRawCalo*process.pfMetTK)
process.e = cms.EndPath(process.out)

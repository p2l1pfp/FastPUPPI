import FWCore.ParameterSet.Config as cms

process = cms.Process("OUT")

# process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(50) )

process.source = cms.Source("PoolSource",
                            # replace 'myfile.root' with the source file you want to use
                            fileNames = cms.untracked.vstring(
        #'/store/mc/TTI2023Upg14D/PYTHIA6_Tauola_TTbar_TuneZ2star_14TeV/GEN-SIM-DIGI-RAW/PU140bx25_PH2_1K_FB_V3-v2/00000/00F0F213-B1E5-E311-9C9A-002618943977.root'
        #'/store/group/dpg_trigger/comm_trigger/L1TrackTrigger/620_SLHC12/Extended2023TTI/SingleTau1p/NoPU/SingleTau1p_E2023TTI_NoPU.root'
        '/store/mc/TTI2023Upg14D/Zmumu_TuneZ2star_14TeV_Eta4-pythia6/GEN-SIM-DIGI-RAW/PU140bx25_PH2_1K_FB_V3-v2/00000/C050E16A-D1E6-E311-A244-0025905A60D2.root',
        '/store/mc/TTI2023Upg14D/Zmumu_TuneZ2star_14TeV_Eta4-pythia6/GEN-SIM-DIGI-RAW/PU140bx25_PH2_1K_FB_V3-v2/00000/CCBD6681-D0E6-E311-9747-0025905964CC.root',
        '/store/mc/TTI2023Upg14D/Zmumu_TuneZ2star_14TeV_Eta4-pythia6/GEN-SIM-DIGI-RAW/PU140bx25_PH2_1K_FB_V3-v2/00000/F4A1119F-DDE6-E311-A9D4-0025905AA9F0.root',
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

process.InfoOut = cms.EDProducer('NtupleProducer',
                                 zeroSuppress = cms.bool(False),
                                 L1TrackTag  = cms.InputTag('TTTracksFromPixelDigis','Level1TTTracks'),
                                 EcalTPTag   = cms.InputTag('simEcalTriggerPrimitiveDigis'),
                                 HcalTPTag   = cms.InputTag('simHcalTriggerPrimitiveDigis'), 
                                 MuonTPTag   = cms.InputTag('simGmtDigis'), 
                                 genParTag   = cms.InputTag('genParticles'),
                                 corrector   = cms.InputTag("/afs/cern.ch/work/n/ntran/private/Correlator/analysis/go7/CMSSW_6_2_0_SLHC12/src/FastPUPPI/NtupleProducer/data/pion_eta_phi.root"),
                                 corrector2  = cms.InputTag("/afs/cern.ch/work/n/ntran/private/Correlator/analysis/go7/CMSSW_6_2_0_SLHC12/src/FastPUPPI/NtupleProducer/data/pion_eta_phi_hpu.root"),
                                 ecorrector  = cms.InputTag("/afs/cern.ch/work/n/ntran/private/Correlator/analysis/go7/CMSSW_6_2_0_SLHC12/src/FastPUPPI/NtupleProducer/data/ecorr.root"),
                                 trackres    = cms.InputTag("/afs/cern.ch/work/n/ntran/private/Correlator/analysis/go7/CMSSW_6_2_0_SLHC12/src/FastPUPPI/NtupleProducer/data/tkres.root"),
                                 eleres      = cms.InputTag("/afs/cern.ch/work/n/ntran/private/Correlator/analysis/go7/CMSSW_6_2_0_SLHC12/src/FastPUPPI/NtupleProducer/data/eres.root"),
                                 pionres     = cms.InputTag("/afs/cern.ch/work/n/ntran/private/Correlator/analysis/go7/CMSSW_6_2_0_SLHC12/src/FastPUPPI/NtupleProducer/data/pionres.root")
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

process.TT_step           = cms.Path(process.TrackTriggerTTTracks)
process.TTAssociator_step = cms.Path(process.TrackTriggerAssociatorTracks)
process.ana               = cms.Path(process.InfoOut*process.pfMet*process.pfMetCalo*process.pfMetRawCalo*process.pfMetTK)
process.p = cms.Schedule(process.TT_step,process.TTAssociator_step,process.ana)
process.e = cms.EndPath(process.out)

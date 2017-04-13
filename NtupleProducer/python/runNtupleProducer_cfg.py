import FWCore.ParameterSet.Config as cms

process = cms.Process("OUT")
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.Geometry.GeometryExtended2023D3Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True), allowUnscheduled = cms.untracked.bool(False) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1))
process.MessageLogger.cerr.FwkReport.reportEvery = 10

process.source = cms.Source("PoolSource",
                            # replace 'myfile.root' with the source file you want to use
                            fileNames = cms.untracked.vstring(
        #'/store/relval/CMSSW_8_1_0_pre15/RelValSingleMuPt100Extended/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D1PU140-v1/10000/204D0758-109A-E611-856C-003048FF9ABC.root'
        #'/store/relval/CMSSW_8_1_0_pre15/RelValSingleMuPt100Extended/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3PU140-v1/10000/0824BBCA-119A-E611-B39E-0025905A48F0.root'
        #'/store/relval/CMSSW_8_1_0_pre15/RelValSingleMuPt100Extended/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D4PU140-v1/10000/025FF191-119A-E611-B262-0025905B8594.root'
        '/store/relval/CMSSW_8_1_0_pre15/RelValZMM_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU140r1-v1/10000/DEF317C3-28A1-E611-9428-0CC47A4C8E2E.root'
        #/store/relval/CMSSW_8_1_0_pre15/RelValZMM_13/GEN-SIM-RECO/81X_upgrade2023_realistic_v3_2023D3Timing13TeV-v1/10000/44798C42-6F99-E611-A504-0CC47A7C354C.root'
        )
)

process.load('FastPUPPI.NtupleProducer.l1tPFCaloProducersFromOfflineRechits_cff')
process.load('FastPUPPI.NtupleProducer.l1tPFTkProducersFromOfflineTracks_cfi')

process.InfoOut = cms.EDProducer('NtupleProducer',
                                 L1TrackTag  = cms.InputTag('l1tPFTkProducersFromOfflineTracksStrips'),
                                 EcalTPTag   = cms.InputTag('l1tPFEcalProducerFromOfflineRechits','towers'),
                                 HGEcalTPTag = cms.InputTag('l1tPFHGCalEEProducerFromOfflineRechits','towers'),
                                 HcalTPTag   = cms.InputTag('l1tPFHcalProducerFromOfflineRechits','towers'),
                                 HGHcalTPTag = cms.InputTag('l1tPFHGCalFHProducerFromOfflineRechits','towers'),
                                 BHHcalTPTag = cms.InputTag('l1tPFHGCalBHProducerFromOfflineRechits','towers'),
                                 HFTPTag     = cms.InputTag('l1tPFHFProducerFromOfflineRechits','towers'),
                                 MuonTPTag   = cms.InputTag('simGmtDigis'), 
                                 genParTag   = cms.InputTag('genParticles'),
                                 zeroSuppress = cms.bool(False),
                                 corrector   = cms.InputTag("/afs/cern.ch/user/p/pharris/pharris/public/bacon/prod/CMSSW_8_1_0_pre16/src/FastPUPPI/NtupleProducer/data/pion_eta_phi.root"),
                                 corrector2  = cms.InputTag("/afs/cern.ch/user/p/pharris/pharris/public/bacon/prod/CMSSW_8_1_0_pre16/src/FastPUPPI/NtupleProducer/data/pion_eta_phi_res_old.root"),
                                 ecorrector  = cms.InputTag("/afs/cern.ch/user/p/pharris/pharris/public/bacon/prod/CMSSW_8_1_0_pre16/src/FastPUPPI/NtupleProducer/data/ecorr.root"),
                                 trackres    = cms.InputTag("/afs/cern.ch/user/p/pharris/pharris/public/bacon/prod/CMSSW_8_1_0_pre16/src/FastPUPPI/NtupleProducer/data/tkres.root"),
                                 eleres      = cms.InputTag("/afs/cern.ch/user/p/pharris/pharris/public/bacon/prod/CMSSW_8_1_0_pre16/src/FastPUPPI/NtupleProducer/data/eres.root"),
                                 pionres     = cms.InputTag("/afs/cern.ch/user/p/pharris/pharris/public/bacon/prod/CMSSW_8_1_0_pre16/src/FastPUPPI/NtupleProducer/data/pionres.root"),
                                 trkPtCut    = cms.double(4.0),
                                 metRate     = cms.bool(True),
                                 etaCharged  = cms.double(2.5),
                                 puppiPtCut  = cms.double(4.0),
                                 vtxRes      = cms.double(0.333),
                                 debug       = cms.untracked.int32(1),
                                 )


process.l1Puppi = cms.Sequence(process.l1tPFCaloProducersFromOfflineRechits+process.l1tPFTkProducersFromOfflineTracksStrips)
process.p = cms.Path(process.l1Puppi*process.InfoOut)
process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string("l1pf_out.root"),
)
#process.e = cms.EndPath(process.out)

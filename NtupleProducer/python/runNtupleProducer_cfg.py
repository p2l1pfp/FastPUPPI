import FWCore.ParameterSet.Config as cms

scenario = "D12"
#scenario = "D4T"
from Configuration.StandardSequences.Eras import eras

if "T" in scenario:
    process = cms.Process("OUT", eras.Phase2C2_timing)
else:
    process = cms.Process("OUT", eras.Phase2C2)

process.load('Configuration.StandardSequences.Services_cff')
if "D12" in scenario: process.load('Configuration.Geometry.GeometryExtended2023D12Reco_cff')
if "D11" in scenario: process.load('Configuration.Geometry.GeometryExtended2023D11Reco_cff')
if "D4"  in scenario: process.load('Configuration.Geometry.GeometryExtended2023D4Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True), allowUnscheduled = cms.untracked.bool(False) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10))
process.MessageLogger.cerr.FwkReport.reportEvery = 10

process.source = cms.Source("PoolSource",
                            # replace 'myfile.root' with the source file you want to use
                            fileNames = cms.untracked.vstring('file:inputs_17D.root')
                            #fileNames = cms.untracked.vstring(
        #'/store/relval/CMSSW_9_1_0_pre1/RelValZMM_14/GEN-SIM-RECO/90X_upgrade2023_realistic_v9_D12-v1/00000/0A8B6505-A215-E711-B25B-0CC47A4C8F26.root'
                            #)
)
process.source.duplicateCheckMode = cms.untracked.string("noDuplicateCheck")

process.load('FastPUPPI.NtupleProducer.l1tPFCaloProducersFromOfflineRechits_cff')
process.load('FastPUPPI.NtupleProducer.l1tPFTkProducersFromOfflineTracks_cfi')
process.load('FastPUPPI.NtupleProducer.l1tPFMuProducerFromL1Mu_cfi')
process.load('FastPUPPI.NtupleProducer.ntupleProducer_cfi')

process.l1Puppi = cms.Sequence(process.l1tPFCaloProducersFromOfflineRechits+process.l1tPFTkProducersFromOfflineTracksStrips+process.l1tPFMuProducerFromL1Mu)

process.p = cms.Path(process.l1Puppi * process.InfoOut )

def goMuonGun(process):
    process.load('FastPUPPI.NtupleProducer.muonGunSelector_cfi')
    process.l1tPFTkProducersFromOfflineTracksStrips.MuonGunVeto = True 
    process.l1tPFMuProducerFromL1Mu.MuonGunVeto = True
    process.l1tPFTkProducersFromOfflineTracksStrips.GenTagForVeto = cms.InputTag("genMuonsFromGun")
    process.l1tPFMuProducerFromL1Mu.GenTagForVeto = cms.InputTag("genMuonsFromGun")
    process.p = cms.Path(process.genMuonsFromGun + process.l1Puppi * process.InfoOut )
    
    
if False: # turn on CMSSW downstream processing and output
    process.InfoOut.outputName = ""; # turn off Ntuples
    process.InfoOut.debug = 0;

    from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets
    process.ak4L1RawCalo = ak4PFJets.clone(src = 'InfoOut:RawCalo')
    process.ak4L1Calo    = ak4PFJets.clone(src = 'InfoOut:Calo')
    process.ak4L1TK      = ak4PFJets.clone(src = 'InfoOut:TK')
    process.ak4L1PF      = ak4PFJets.clone(src = 'InfoOut:PF')
    process.ak4L1Puppi   = ak4PFJets.clone(src = 'InfoOut:Puppi')

    from RecoMET.METProducers.PFMET_cfi import pfMet
    pfMet.calculateSignificance = False
    process.l1MetRawCalo = pfMet.clone(src = "InfoOut:RawCalo")
    process.l1MetCalo    = pfMet.clone(src = "InfoOut:Calo")
    process.l1MetTK      = pfMet.clone(src = "InfoOut:TK")
    process.l1MetPF      = pfMet.clone(src = "InfoOut:PF")
    process.l1MetPuppi   = pfMet.clone(src = "InfoOut:Puppi")

    process.l1JetMET = cms.Sequence(
        process.ak4L1RawCalo + process.ak4L1Calo + process.ak4L1TK + process.ak4L1PF + process.ak4L1Puppi +
        process.l1MetRawCalo + process.l1MetCalo + process.l1MetTK + process.l1MetPF + process.l1MetPuppi
    )

    process.p = cms.Path(process.l1Puppi * process.InfoOut * process.l1JetMET )

    process.out = cms.OutputModule("PoolOutputModule",
            fileName = cms.untracked.string("l1pf_out.root"),
            outputCommands = cms.untracked.vstring("drop *",
                "keep *_gmtStage2Digis_*_*",
                "keep *_genParticles_*_*",
                "keep *_ak4GenJetsNoNu_*_*",
                "keep *_genMetTrue_*_*",
                "keep *_*_*_OUT",
            )
    )
    process.e = cms.EndPath(process.out)

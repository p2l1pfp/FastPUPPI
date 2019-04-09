import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras

process = cms.Process("ID", eras.Phase2C4_trigger)

process.load('Configuration.StandardSequences.Services_cff')
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1))
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.source = cms.Source("PoolSource",
                           fileNames = cms.untracked.vstring('file:/afs/cern.ch/user/p/pharris/pharris/public/bacon/prod/tmp/CMSSW_10_5_0_pre1/src/L1Trigger/Phase2L1Taus/test/tmp2/tau104X.root',
                                                             'file:/afs/cern.ch/user/p/pharris/pharris/public/bacon/prod/tmp/CMSSW_10_5_0_pre1/src/L1Trigger/Phase2L1Taus/test/tmp2/tau104X_v2.root'
                                                              ),

#                           fileNames = cms.untracked.vstring('file:/afs/cern.ch/user/p/pharris/pharris/public/bacon/prod/tmp/CMSSW_10_5_0_pre1/src/L1Trigger/Phase2L1Taus/test/tau104X.root'),
                            duplicateCheckMode = cms.untracked.string("noDuplicateCheck")
)
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '103X_upgrade2023_realistic_v2', '') 


ntuple = cms.EDAnalyzer("TauNTuplizer",
    src = cms.InputTag("L1PFTauProducer","L1PFTaus"),
    genParticles = cms.InputTag("genParticles"),
    drMax = cms.double(0.2),
    minRecoPtOverGenPt = cms.double(0.2),
    onlyMatched = cms.bool(False),
    variables = cms.PSet(
	pt = cms.string("pt"),
 	eta = cms.string("eta"),
 	phi = cms.string("phi"),
 	chargedIso = cms.string("chargedIso"),
 	fullIso    = cms.string("fullIso"),
        decayId    = cms.string("id"),
    ),
)

process.ntuplePF   = ntuple.clone()
process.ntupleNN   = ntuple.clone(src=cms.InputTag("L1NNTauProducer","L1PFTausNN") )
process.ntuplePF2  = ntuple.clone(src=cms.InputTag("L1PFTauProducer2","L1PFTaus") )
process.ntupleNN2  = ntuple.clone(src=cms.InputTag("L1NNTauProducer2","L1PFTausNN") )
process.ntupleNN3  = ntuple.clone(src=cms.InputTag("L1NNTauProducer3","L1PFTausNN") )
  
modules = [
    process.ntuplePF,
    process.ntupleNN,
    process.ntuplePF2,
    process.ntupleNN2,
    process.ntupleNN3
]

process.p = cms.Path(sum(modules[1:], modules[0]))
process.TFileService = cms.Service("TFileService", fileName = cms.string("idTupleNew.root"))


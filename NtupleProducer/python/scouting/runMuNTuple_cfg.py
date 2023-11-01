import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras
from PhysicsTools.NanoAOD.common_cff import Var, ExtVar 

process = cms.Process("RESP", eras.Phase2C17I13M9)

process.load('Configuration.StandardSequences.Services_cff')
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options = cms.untracked.PSet(
            wantSummary = cms.untracked.bool(True),
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000))
process.MessageLogger.cerr.FwkReport.reportEvery = 10

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('/store/cmst3/group/l1tr/gpetrucc/12_5_X/NewInputs125X/150223/DYToLL_M-50_PU200/inputs125X_1.root')
)

process.load('Configuration.Geometry.GeometryExtended2026D88Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2026D88_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff') # needed to read HCal TPs
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '125X_mcRun4_realistic_v2', '')

process.extraPFStuff = cms.Task(
        process.l1tTkMuonsGmt,
        process.l1tSAMuonsGmt,
        process.l1tGTTInputProducer,
        process.l1tVertexFinderEmulator,
        process.L1TLayer1TaskInputsTask,
        process.L1TLayer1Task,
        process.L1TPFJetsEmulationTask
)


process.tkMuTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
    src = cms.InputTag("l1tTkMuonsGmt"),
    cut = cms.string(""),
    name = cms.string("TkMu"),
    doc = cms.string("TkMuons from GMT"),
    singleton = cms.bool(False), # the number of entries is variable
    extension = cms.bool(False), # this is the main table
    variables = cms.PSet(
        pt   = Var("phPt",  float),
        eta  = Var("phEta", float),
        phi  = Var("phPhi", float),
        mass = Var("0.10566", float),
        vz   = Var("phZ0",  float, doc="Z coordinate of the reconstructed production vertex"),
        d0   = Var("phD0",  float, doc="transverse impact parameter (always zero currently)"),
        charge = Var("phCharge", int, doc="charge"),
        quality = Var("hwQual", int, doc="quality (TBD)"),
        hwPt   = Var("hwPt", int),
        hwEta  = Var("hwEta", int),
        hwPhi  = Var("hwPhi", int),
        hwZ0  = Var("hwZ0", int),
        hwD0  = Var("hwD0", int),
        hwCharge  = Var("hwCharge", int),
        hwIsoSum  = Var("hwIsoSum", int),
        hwIsoSumAp  = Var("hwIsoSumAp", int),
        hwQual  = Var("hwQual", int),
    )
)

process.genMu = cms.EDFilter("GenParticleSelector",
    src = cms.InputTag("genParticles"),
    cut = cms.string("abs(pdgId) = 13 && status == 1 && pt > 0.5 && abs(eta) < 2.7"),
)

process.genMuTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
    src = cms.InputTag("genMu"),
    cut = cms.string(""),
    name = cms.string("GenMu"),
    doc = cms.string("gen muons"),
    singleton = cms.bool(False), # the number of entries is variable
    extension = cms.bool(False), # this is the main table
    variables = cms.PSet(
        pt  = Var("pt",  float),
        phi = Var("phi", float),
        eta  = Var("eta", float),
        mass  = Var("mass", float),
        vz   = Var("vz",  float, doc="Production point along the beam axis"),
        d0   = Var("vertex.Rho",  float, doc="transverse distance of production point from the beam axis"),
        charge  = Var("charge", int, doc="charge"),
        isPrompt  = Var("statusFlags().isPrompt()", int, doc="Prompt muon"),
        isFromTau  = Var("statusFlags().isDirectPromptTauDecayProduct()", int, doc="Muon from prompt tau decay"),
    )
)
process.tkMuMCMatch = cms.EDProducer("MCMatcher",       # cut on deltaR, deltaPt/Pt; pick best by deltaR
    src         = process.tkMuTable.src,                         # final reco collection
    matched     = cms.InputTag("genMu"),     # final mc-truth particle collection
    mcPdgId     = cms.vint32(13),               # one or more PDG ID (13 = mu); absolute values (see below)
    checkCharge = cms.bool(False),              # True = require RECO and MC objects to have the same charge
    mcStatus    = cms.vint32(1),                # PYTHIA status code (1 = stable, 2 = shower, 3 = hard scattering)
    maxDeltaR   = cms.double(0.1),              # Minimum deltaR for the match
    maxDPtRel   = cms.double(0.5),              # Minimum deltaPt/Pt for the match
    resolveAmbiguities    = cms.bool(True),     # Forbid two RECO objects to match to the same GEN object
    resolveByMatchQuality = cms.bool(True),    # False = just match input in order; True = pick lowest deltaR pair first
)

process.tkMuMCTable = cms.EDProducer("CandMCMatchTableProducer",
    src     = process.tkMuTable.src,
    mcMap   = cms.InputTag("tkMuMCMatch"),
    objName = process.tkMuTable.name,
    objType = cms.string("Other"),
    branchName = cms.string("GenMu"),
    docString = cms.string("MC matching"),
)
process.extraPFStuff.add(
        process.genMu,
        process.tkMuMCMatch,
)

process.p = cms.Path(
    process.tkMuTable  + process.genMuTable + process.tkMuMCTable
)
process.p.associate(process.extraPFStuff)

process.outnano = cms.OutputModule("NanoAODOutputModule",
    fileName = cms.untracked.string("l1MuNano.root"),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('p')),
    outputCommands = cms.untracked.vstring("drop *", "keep nanoaodFlatTable_*Table_*_*"),
    compressionLevel = cms.untracked.int32(4),
    compressionAlgorithm = cms.untracked.string("ZLIB"),
)
process.end = cms.EndPath(process.outnano)
process.schedule = cms.Schedule(process.p, process.end)
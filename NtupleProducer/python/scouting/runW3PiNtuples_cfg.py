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
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10))
process.MessageLogger.cerr.FwkReport.reportEvery = 10

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:/tmp/gpetrucc/005E74D6-B50E-674E-89E6-EAA9A617B476.root')
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

process.l1tTrackSelectionProducer.processSimulatedTracks = False # these would need stubs, and are not used anyway

process.extraPFStuff = cms.Task(
        process.l1tSAMuonsGmt,
        process.l1tGTTInputProducer,
        process.l1tTrackSelectionProducer,
        process.l1tVertexFinderEmulator,
        process.L1TLayer1TaskInputsTask,
        process.L1TLayer1Task,
        process.L1TPFJetsEmulationTask
)

process.l1tLayer2Deregionizer.nPuppiFinalBuffer = 200

process.puppiTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
                src = cms.InputTag("l1tLayer2Deregionizer:Puppi"),
                cut = cms.string(""),
                name = cms.string("L1Puppi"),
                doc = cms.string("L1Puppi candidates"),
                singleton = cms.bool(False), # the number of entries is variable
                extension = cms.bool(False), # this is the main table
                variables = cms.PSet(
                    pt   = Var("pt",  float,precision=8),
                    phi  = Var("phi", float,precision=8),
                    eta  = Var("eta", float,precision=8),
                    mass = Var("mass", float,precision=8),
                    vz   = Var("vz",  float,precision=8),
                    charge = Var("charge", int, doc="charge"),
                    pdgId  = Var("pdgId", int, doc="PDG id"),
                )
            )

process.p = cms.Path(
    process.puppiTable
)
process.p.associate(process.extraPFStuff)

def goSignal():
    ## Add to the configuration:
    ##   * the generator-level W boson and pions
    ##   * for each Puppi candidate, information on whether it's matched to a gen pion
    ##   * the table of PF candidates, so that one can 
    process.genW = cms.EDFilter("GenParticleSelector",
            src = cms.InputTag("genParticles"),
            cut = cms.string("abs(pdgId) = 24 && numberOfDaughters = 3 && abs(daughter(0).pdgId) = 211 && abs(daughter(1).pdgId) = 211 && abs(daughter(2).pdgId) = 211"),
            filter = cms.bool(True),
    )
    process.genPiFromW = cms.EDFilter("GenParticleSelector",
            src = cms.InputTag("genParticles"),
            cut = cms.string("abs(pdgId) = 211 && numberOfMothers > 0 && abs(motherRef.pdgId) = 24"),
    )
    process.genWTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
                    src = cms.InputTag("genW"),
                    cut = cms.string(""),
                    name = cms.string("GenW"),
                    doc = cms.string("gen W boson"),
                    singleton = cms.bool(False), # the number of entries is variable
                    extension = cms.bool(False), # this is the main table
                    variables = cms.PSet(
                        pt  = Var("pt",  float,precision=8),
                        phi = Var("phi", float,precision=8),
                        eta  = Var("eta", float,precision=8),
                        mass  = Var("mass", float,precision=8),
                        vz   = Var("vz",  float,precision=8, doc="Production point along the beam axis"),
                        charge  = Var("charge", int, doc="charge"),
                        pdgId  = Var("pdgId", int, doc="PDG id"),
                    )
                )
    process.genPiTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
                    src = cms.InputTag("genPiFromW"),
                    cut = cms.string(""),
                    name = cms.string("GenPi"),
                    doc = cms.string("gen pions from W boson decay"),
                    singleton = cms.bool(False), # the number of entries is variable
                    extension = cms.bool(False), # this is the main table
                    variables = cms.PSet(
                        pt  = Var("pt",  float,precision=8),
                        phi = Var("phi", float,precision=8),
                        eta  = Var("eta", float,precision=8),
                        mass  = Var("mass", float,precision=8),
                        vz   = Var("vz",  float,precision=8, doc="Production point along the beam axis"),
                        charge  = Var("charge", int, doc="charge"),
                        pdgId  = Var("pdgId", int, doc="PDG id"),
                        prompt  = Var("2*statusFlags().isPrompt() + statusFlags().isDirectPromptTauDecayProduct()", int, doc="Particle status."),
                    )
                )
    process.pfTable = process.puppiTable.clone(
            src = cms.InputTag("l1ctLayer1","PF"),
            name = cms.string("L1PF"),
            doc = cms.string("L1PF candidates"),
    )
    process.puppiMCMatch = cms.EDProducer("MCMatcher",       # cut on deltaR, deltaPt/Pt; pick best by deltaR
        src         = process.puppiTable.src,                         # final reco collection
        matched     = cms.InputTag("genPiFromW"),     # final mc-truth particle collection
        mcPdgId     = cms.vint32(211),               # one or more PDG ID (13 = mu); absolute values (see below)
        checkCharge = cms.bool(False),              # True = require RECO and MC objects to have the same charge
        mcStatus    = cms.vint32(1),                # PYTHIA status code (1 = stable, 2 = shower, 3 = hard scattering)
        maxDeltaR   = cms.double(0.1),              # Minimum deltaR for the match
        maxDPtRel   = cms.double(0.5),              # Minimum deltaPt/Pt for the match
        resolveAmbiguities    = cms.bool(True),     # Forbid two RECO objects to match to the same GEN object
        resolveByMatchQuality = cms.bool(True),    # False = just match input in order; True = pick lowest deltaR pair first
    )
    process.pfMCMatch = process.puppiMCMatch.clone(
        src = process.pfTable.src 
    )
    process.puppiMCTable = cms.EDProducer("CandMCMatchTableProducer",
        src     = process.puppiTable.src,
        mcMap   = cms.InputTag("puppiMCMatch"),
        objName = process.puppiTable.name,
        objType = cms.string("Other"),
        branchName = cms.string("GenPi"),
        docString = cms.string("MC matching"),
    )
    process.pfMCTable = process.puppiMCTable.clone(
        src     = process.pfTable.src,
        mcMap   = cms.InputTag("pfMCMatch"),
        objName = process.pfTable.name
    )
    process.extraPFStuff.add(
            process.genW,
            process.genPiFromW,
            process.puppiMCMatch,
    )
    process.p.replace(process.puppiTable, process.puppiTable + process.genWTable + process.genPiTable + process.puppiMCTable)

def noPU():
    for X in "", "EM", "Raw", "EMRaw":
        pfc = getattr(process, 'l1tPFClustersFromHGC3DClusters'+X, None)
        if not pfc: continue
        pfc.emVsPUID.wp = "-1.0"

process.outnano = cms.OutputModule("NanoAODOutputModule",
    fileName = cms.untracked.string("l1Nano.root"),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('p')),
    outputCommands = cms.untracked.vstring("drop *", "keep nanoaodFlatTable_*Table_*_*"),
    compressionLevel = cms.untracked.int32(4),
    compressionAlgorithm = cms.untracked.string("ZLIB"),
)
process.end = cms.EndPath(process.outnano)
process.schedule = cms.Schedule(process.p, process.end)

process.source.fileNames = [
        '/store/cmst3/group/l1tr/gpetrucc/12_5_X/NewInputs125X/150223/WTo3Pion_PU200/inputs125X_WTo3Pion_PU200_job1.root',
]

#goSignal()

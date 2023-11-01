import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras
from PhysicsTools.NanoAOD.common_cff import Var, ExtVar 

process = cms.Process("L1Dump", eras.Phase2C17I13M9)

process.load('Configuration.StandardSequences.Services_cff')
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100))
process.MessageLogger.cerr.FwkReport.reportEvery = 10

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:/tmp/gpetrucc/005E74D6-B50E-674E-89E6-EAA9A617B476.root')
)

process.load('Configuration.Geometry.GeometryExtended2026D88Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2026D88_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff') # needed to read HCal TPs
process.load('SimCalorimetry.HGCalSimProducers.hgcalDigitizer_cfi') # needed for HGCAL_noise_fC
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '125X_mcRun4_realistic_v2', '')

process.l1tTrackSelectionProducer.processSimulatedTracks = False # these would need stubs, and are not used anyway

process.deps = cms.Task(
    process.l1tTkMuonsGmt,
    process.l1tSAMuonsGmt,
    process.l1tGTTInputProducer,
    process.l1tTrackSelectionProducer,
    process.l1tVertexFinderEmulator,
    process.L1TLayer1TaskInputsTask,
    process.L1TLayer1Task,
    process.l1tLayer2EG,
    process.L1TPFJetsEmulationTask,
    process.L1TPFJetsExtendedTask,
    process.L1TBJetsTask
)

process.puppiDump = cms.EDAnalyzer("L1PuppiBinaryDumper",
                src = cms.InputTag("l1tLayer2Deregionizer:Puppi"),
                outName = cms.string("puppi.dump"))

process.jetDump = cms.EDAnalyzer("L1JetBinaryDumper",
                src = cms.InputTag("l1tSCPFL1PuppiExtendedCorrectedEmulator"),
                ptMin = cms.double(15),
                outName = cms.string("puppiJets.dump"))

process.tkMuDump = cms.EDAnalyzer("L1TrackerMuonBinaryDumper",
                src = cms.InputTag("l1tTkMuonsGmt"),
                outName = cms.string("tkMuons.dump"))

process.egDump = cms.EDAnalyzer("L1CTL2EgammaBinaryDumper",
                srcEle = cms.InputTag("l1tLayer2EG", "L1CtTkElectron"),
                srcEm = cms.InputTag("l1tLayer2EG", "L1CtTkEm"),
                interleaveOutputs = cms.bool(False), # False = first 12 photons, then electrons; True = pho1, ele1, pho2, ele2, ...
                outName = cms.string("egamma.dump"))

process.p_dumps = cms.EndPath(
    process.puppiDump +
    process.jetDump +
    process.tkMuDump +
    process.egDump
)
process.p_dumps.associate(process.deps)

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
        z0   = Var("phZ0",  float, doc="Z coordinate of the reconstructed production vertex"),
        dxy   = Var("phD0",  float, doc="transverse impact parameter (always zero currently)"),
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
        z0   = Var("vz",  float, doc="Production point along the beam axis"),
        dxy   = Var("vertex.Rho",  float, doc="transverse distance of production point from the beam axis"),
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

process.p_mu = cms.Path(process.tkMuTable)
process.p_muMC = cms.Path(process.genMu + process.genMuTable + process.tkMuMCMatch + process.tkMuMCTable)

process.puppiTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
        src = cms.InputTag("l1tLayer2Deregionizer:Puppi"),
        cut = cms.string(""),
        name = cms.string("Puppi"),
        doc = cms.string("L1Puppi candidates"),
        singleton = cms.bool(False), # the number of entries is variable
        extension = cms.bool(False), # this is the main table
        variables = cms.PSet(
            pt   = Var("pt",  float),
            phi  = Var("phi", float),
            eta  = Var("eta", float),
            mass = Var("mass", float),
            z0   = Var("vz",  float),
            charge = Var("charge", int, doc="charge"),
            pdgId  = Var("pdgId", int, doc="PDG id"),
        )
    )

process.genW = cms.EDFilter("GenParticleSelector",
        src = cms.InputTag("genParticles"),
        cut = cms.string("abs(pdgId) = 24 && numberOfDaughters = 3 && abs(daughter(0).pdgId) = 211 && abs(daughter(1).pdgId) = 211 && abs(daughter(2).pdgId) = 211"),
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
                    z0   = Var("vz",  float,precision=8, doc="Production point along the beam axis"),
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
                    z0   = Var("vz",  float,precision=8, doc="Production point along the beam axis"),
                    charge  = Var("charge", int, doc="charge"),
                    pdgId  = Var("pdgId", int, doc="PDG id"),
                    prompt  = Var("2*statusFlags().isPrompt() + statusFlags().isDirectPromptTauDecayProduct()", int, doc="Particle status."),
                )
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
process.puppiMCTable = cms.EDProducer("CandMCMatchTableProducer",
    src     = process.puppiTable.src,
    mcMap   = cms.InputTag("puppiMCMatch"),
    objName = process.puppiTable.name,
    objType = cms.string("Other"),
    branchName = cms.string("GenPi"),
    docString = cms.string("MC matching"),
)

process.p_puppi = cms.Path(process.puppiTable)
process.p_puppiMC = cms.Path(process.genW + process.genPiFromW + process.genWTable + process.genPiTable + process.puppiMCMatch + process.puppiMCTable)

process.puppiJetsTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
        src = cms.InputTag("l1tSCPFL1PuppiExtendedCorrectedEmulator"),
        cut  = cms.string("pt > 10"),
        name = cms.string("PuppiJet"),
        doc = cms.string("Puppi Jets reconstructed by L1T (extended tracking)"),
        singleton = cms.bool(False), # the number of entries is variable
        extension = cms.bool(False), # this is the main table
        variables = cms.PSet(
            pt  = Var("pt",  float),
            phi = Var("phi", float),
            eta  = Var("eta", float),
            mass  = Var("mass", float),
        ),
        externalVariables = cms.PSet(
            btagScore = ExtVar(cms.InputTag("l1tBJetProducerPuppiCorrectedEmulator", "L1PFBJets"), float, doc="b-tag score (L1 model)"),
        )
)

process.genJetsTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
        src = cms.InputTag("ak4GenJetsNoNu"),
        cut  = cms.string("pt > 10 && abs(eta) < 5.0"),
        name = cms.string("GenJet"),
        doc = cms.string("GenJets (ak4)"),
        singleton = cms.bool(False), # the number of entries is variable
        extension = cms.bool(False), # this is the main table
        variables = cms.PSet(
            pt  = Var("pt",  float),
            phi = Var("phi", float),
            eta  = Var("eta", float),
            mass  = Var("mass", float),
        )
)

process.p_jets = cms.Path(process.puppiJetsTable)
process.p_jetsMC = cms.Path(process.genJetsTable)

process.genVertexTable = cms.EDProducer("SimpleXYZPointFlatTableProducer",
    src = cms.InputTag("genParticles:xyz0"),
    cut = cms.string(""),
    name= cms.string("GenVtx"),
    doc = cms.string("Gen vertex"),
    variables = cms.PSet(
         x = Var("x", float, doc="gen vertex x", precision=10),
         y = Var("y", float, doc="gen vertex y", precision=10),
         z = Var("z", float, doc="gen vertex z", precision=16),
    )
)

process.l1VertexTable = cms.EDProducer("VertexWordFlatTableProducer",
        name = cms.string("L1Vtx"),
        cut  = cms.string(""),
        src = cms.InputTag("l1tVertexFinderEmulator","l1verticesEmulation"),
        doc = cms.string("Primary vertices reconstructed by L1T"),
        singleton = cms.bool(False), # the number of entries is variable
        extension = cms.bool(False), # this is the main table
        variables = cms.PSet(
            sumpt = Var("pt",  float,precision=10),
            z     = Var("z0",  float,precision=16),
        )
)
process.p_pv = cms.Path(process.l1VertexTable)
process.p_pvMC = cms.Path(process.genVertexTable)

process.phoTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
    src = cms.InputTag("l1tLayer2EG","L1CtTkEm"),
    cut = cms.string(""),
    name = cms.string("Pho"),
    doc = cms.string("Photons (TkEm) from CTL2"),
    singleton = cms.bool(False), # the number of entries is variable
    extension = cms.bool(False), # this is the main table
    variables = cms.PSet(
        pt   = Var("pt",  float),
        eta  = Var("eta", float),
        phi  = Var("phi", float),
        mass = Var("0", float),
        quality = Var("hwQual", int, doc="quality (TBD)"),
        trkIsol = Var("trkIsol", float),
        trkIsolPV = Var("trkIsolPV", float),
        puppiIsol = Var("puppiIsol", float),
        puppiIsolPV = Var("puppiIsolPV", float),
        hwPt   = Var("hwPt", int),
        hwEta  = Var("hwEta", int),
        hwPhi  = Var("hwPhi", int),
        hwQual  = Var("hwQual", int),
    )
)

process.eleTable = process.phoTable.clone(
    src = "l1tLayer2EG:L1CtTkElectron",
    cut = "",
    name = "Ele",
    doc = "TkElectrons from CTL2",
    variables = dict(
        z0   = Var("trkzVtx",  float),
        idScore = Var("idScore",  float),
    )
)
process.genPho = cms.EDFilter("GenParticleSelector",
    src = cms.InputTag("genParticles"),
    cut = cms.string("pdgId == 22 && status == 1 && (pt > 10 || isPromptFinalState())"),
)
process.genEle = process.genPho.clone(
    cut = "pdgId == 11 && status == 1",
)
process.genPhoTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
    src = cms.InputTag("genPho"),
    cut = cms.string(""),
    name = cms.string("GenPho"),
    doc = cms.string("gen photons (pt > 10 || prompt)"),
    singleton = cms.bool(False), # the number of entries is variable
    extension = cms.bool(False), # this is the main table
    variables = cms.PSet(
        pt  = Var("pt",  float),
        phi = Var("phi", float),
        eta  = Var("eta", float),
        mass  = Var("mass", float),
        z0   = Var("vz",  float, doc="Production point along the beam axis"),
        dxy   = Var("vertex.Rho",  float, doc="transverse distance of production point from the beam axis"),
        isPrompt  = Var("statusFlags().isPrompt()", int, doc="Prompt"),
        isFromTau  = Var("statusFlags().isDirectPromptTauDecayProduct()", int, doc="Electron from prompt tau decay"),
        motherId  = Var("? motherRef.isNonnull() ? motherRef.pdgId : 0", int, doc="Mother particle ID"),
    )
)
process.genEleTable = process.genPhoTable.clone(
    src = "genEle",
    cut = "",
    name = "GenEle",
    doc = "gen photons (pt > 10 || prompt)",
)
process.phoMCMatch = cms.EDProducer("MCMatcher",       # cut on deltaR, deltaPt/Pt; pick best by deltaR
    src         = process.phoTable.src,                         # final reco collection
    matched     = cms.InputTag("genPho"),     # final mc-truth particle collection
    mcPdgId     = cms.vint32(22),               # one or more PDG ID (13 = mu); absolute values (see below)
    checkCharge = cms.bool(False),              # True = require RECO and MC objects to have the same charge
    mcStatus    = cms.vint32(1),                # PYTHIA status code (1 = stable, 2 = shower, 3 = hard scattering)
    maxDeltaR   = cms.double(0.2),              # Minimum deltaR for the match
    maxDPtRel   = cms.double(0.5),              # Minimum deltaPt/Pt for the match
    resolveAmbiguities    = cms.bool(True),     # Forbid two RECO objects to match to the same GEN object
    resolveByMatchQuality = cms.bool(True),    # False = just match input in order; True = pick lowest deltaR pair first
)
process.eleMCMatch = process.phoMCMatch.clone(
    src         = process.eleTable.src,  # final reco collection
    matched     = "genEle",              # final mc-truth particle collection
    mcPdgId     = [11],                  # one or more PDG ID (13 = mu); absolute values (see below)
)

process.phoMCTable = cms.EDProducer("CandMCMatchTableProducer",
    src     = process.phoTable.src,
    mcMap   = cms.InputTag("phoMCMatch"),
    objName = process.phoTable.name,
    objType = cms.string("Other"),
    branchName = cms.string("GenPho"),
    docString = cms.string("MC matching"),
)
process.eleMCTable = process.phoMCTable.clone(
    src     = process.eleTable.src,
    mcMap   = "eleMCMatch",
    objName = process.eleTable.name,
    branchName = "GenEle",
)
process.p_pho = cms.Path(process.phoTable)
process.p_phoMC = cms.Path(process.genPho + process.genPhoTable + process.phoMCMatch + process.phoMCTable)
process.p_ele = cms.Path(process.eleTable)
process.p_eleMC = cms.Path(process.genEle + process.genEleTable + process.eleMCMatch + process.eleMCTable)

process.outnano = cms.OutputModule("NanoAODOutputModule",
    fileName = cms.untracked.string("l1Nano.root"),
    outputCommands = cms.untracked.vstring("drop *", "keep nanoaodFlatTable_*Table_*_*"),
    compressionLevel = cms.untracked.int32(4),
    compressionAlgorithm = cms.untracked.string("ZLIB"),
)
process.end = cms.EndPath(process.outnano)

def noPU():
    for X in "", "EM", "Raw", "EMRaw":
        pfc = getattr(process, 'pfClustersFromHGC3DClusters'+X, None)
        if not pfc: continue
        pfc.emVsPUID.wp = "-1.0"

def noNano():
    process.schedule = cms.Schedule(process.p_dumps)

process.source.fileNames = [
    '/store/cmst3/group/l1tr/gpetrucc/12_5_X/NewInputs125X/150223/TTbar_PU200/inputs125X_1.root',
]

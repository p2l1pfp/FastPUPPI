import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras
from PhysicsTools.NanoAOD.common_cff import Var, ExtVar 

process = cms.Process("RESP", eras.Phase2C17I13M9)

process.load('Configuration.StandardSequences.Services_cff')
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True), allowUnscheduled = cms.untracked.bool(False) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1))
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:inputs110X.root'),
    inputCommands = cms.untracked.vstring("keep *", 
            "drop l1tPFClusters_*_*_*",
            "drop l1tPFTracks_*_*_*",
            "drop l1tPFCandidates_*_*_*",
            "drop l1tTkPrimaryVertexs_*_*_*")
)

process.load('Configuration.Geometry.GeometryExtended2026D88Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2026D88_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff') # needed to read HCal TPs
process.load('SimCalorimetry.HGCalSimProducers.hgcalDigitizer_cfi') # needed for HGCAL_noise_fC
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('RecoMET.Configuration.GenMETParticles_cff')
process.load('RecoMET.METProducers.genMetTrue_cfi')


from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets
from RecoMET.METProducers.pfMet_cfi import pfMet

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '125X_mcRun4_realistic_v2', '')

process.l1tTrackSelectionProducer.processSimulatedTracks = False # these would need stubs, and are not used anyway

process.extraPFStuff = cms.Task(
        process.l1tTkMuonsGmt,
        process.genParticlesForMETAllVisible,
        process.l1tSAMuonsGmt,
        process.l1tGTTInputProducer,
        process.l1tTrackSelectionProducer,
        process.l1tVertexFinderEmulator,
        process.L1TLayer1TaskInputsTask,
        process.L1TLayer1Task,
        process.L1TPFJetsEmulationTask,
        process.L1TPFJetsExtendedTask,
        process.L1TBJetsTask
)


process.jetTable = cms.EDProducer("L1PFJetTableProducer",
    gen = cms.InputTag("ak4GenJetsNoNu"),
    commonSel = cms.string("pt > 10 && abs(eta) < 5.0"),
    drMax = cms.double(0.2),
    minRecoPtOverGenPt = cms.double(0.1),
    jets = cms.PSet(
        Gen = cms.InputTag("ak4GenJetsNoNu"),
    ),
    moreVariables = cms.PSet(
    ),
)
process.puppiJetsTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
        name = cms.string("L1PuppiJets"),
        cut  = cms.string(""),
        src = cms.InputTag("l1tSCPFL1PuppiCorrectedEmulator"),
        doc = cms.string("Puppi Jets reconstructed by L1T"),
        singleton = cms.bool(False), # the number of entries is variable
        extension = cms.bool(False), # this is the main table
        variables = cms.PSet(
            pt  = Var("pt",  float),
            phi = Var("phi", float),
            eta  = Var("eta", float),
            mass  = Var("mass", float),
        )
)
process.puppiExtJetsTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
        name = cms.string("L1PuppiExtJets"),
        cut  = cms.string(""),
        src = cms.InputTag("l1tSCPFL1PuppiExtendedCorrectedEmulator"),
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

process.genParticleTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
        name = cms.string("GenParticles"),
        cut  = cms.string("abs(eta) < 5"),
        src = cms.InputTag("genParticlesForMETAllVisible"),
        doc = cms.string("visible gen particles"),
        singleton = cms.bool(False), # the number of entries is variable
        extension = cms.bool(False), # this is the main table
        variables = cms.PSet(
            pt  = Var("pt",  float,precision=8),
            phi = Var("phi", float,precision=8),
            eta  = Var("eta", float,precision=8),
            vz   = Var("vz",  float,precision=8),
            pdgId  = Var("pdgId", int, doc="charge"),
            charge  = Var("charge", int, doc="charge"),
            prompt  = Var("2*statusFlags().isPrompt() + statusFlags().isDirectPromptTauDecayProduct()", int, doc="Particle status."),
        )
)

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


process.puppiParticleTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
        name = cms.string("L1PuppiParticles"),
        cut  = cms.string(""),
        src = cms.InputTag("l1tLayer2Deregionizer:Puppi"),
        doc = cms.string("Particles reconstructed by L1T"),
        singleton = cms.bool(False), # the number of entries is variable
        extension = cms.bool(False), # this is the main table
        variables = cms.PSet(
            pt  = Var("pt",  float,precision=8),
            phi = Var("phi", float,precision=8),
            eta  = Var("eta", float,precision=8),
            vz   = Var("vz",  float,precision=8),
            pdgId  = Var("pdgId", int, doc="charge"),
            charge  = Var("charge", int, doc="charge"),
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
    process.l1tBJetProducerPuppi +
    process.jetTable +
    process.puppiJetsTable +
    process.puppiExtJetsTable +
    process.genParticleTable +
    process.genVertexTable +
    process.l1VertexTable +
    process.puppiParticleTable +
    process.tkMuTable  + process.genMuTable + process.tkMuMCTable
)
process.p.associate(process.extraPFStuff)

process.outnano = cms.OutputModule("NanoAODOutputModule",
    fileName = cms.untracked.string("jetNTuple.root"),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('p')),
    outputCommands = cms.untracked.vstring("drop *", "keep nanoaodFlatTable_*Table_*_*"),
    compressionLevel = cms.untracked.int32(4),
    compressionAlgorithm = cms.untracked.string("ZLIB"),
)
process.end = cms.EndPath(process.outnano)

def noPU():
    for X in "", "EM", "Raw", "EMRaw":
        pfc = getattr(process, 'l1tPFClustersFromHGC3DClusters'+X, None)
        if not pfc: continue
        pfc.emVsPUID.wp = "-1.0"

def oldInputs_11_1_6():
    process.l1tTowerCalibration.l1CaloTowers = cms.InputTag("L1EGammaClusterEmuProducer","")
    process.l1tTowerCalibration.L1HgcalTowersInputTag = cms.InputTag("hgcalTowerProducer", "HGCalTowerProcessor")
    process.l1tPFClustersFromL1EGClusters.src = cms.InputTag("L1EGammaClusterEmuProducer",)
    process.l1tPFClustersFromCombinedCaloHCal.phase2barrelCaloTowers = [cms.InputTag("L1EGammaClusterEmuProducer",)]
    process.l1tPFClustersFromHGC3DClusters.src  = cms.InputTag("hgcalBackEndLayer2Producer","HGCalBackendLayer2Processor3DClustering")
    process.l1tPFClustersFromCombinedCaloHF.hcalCandidates = [ cms.InputTag("hgcalBackEndLayer2Producer","HGCalBackendLayer2Processor3DClustering")]
    process.l1tPFTracksFromL1Tracks.L1TrackTag = cms.InputTag("TTTracksFromTrackletEmulation","Level1TTTracks")
    process.l1tGTTInputProducer.l1TracksInputTag = cms.InputTag("TTTracksFromTrackletEmulation","Level1TTTracks")

def oldInputs_12_3_X():
    process.l1tTowerCalibration.l1CaloTowers = cms.InputTag("L1EGammaClusterEmuProducer","L1CaloTowerCollection")
    process.l1tTowerCalibration.L1HgcalTowersInputTag = cms.InputTag("hgcalTowerProducer", "HGCalTowerProcessor")
    process.l1tPFClustersFromL1EGClusters.src = cms.InputTag("L1EGammaClusterEmuProducer",)
    process.l1tPFClustersFromCombinedCaloHCal.phase2barrelCaloTowers = [cms.InputTag("L1EGammaClusterEmuProducer","L1CaloTowerCollection")]
    process.l1tPFClustersFromHGC3DClusters.src  = cms.InputTag("hgcalBackEndLayer2Producer","HGCalBackendLayer2Processor3DClustering")
    process.l1tPFClustersFromCombinedCaloHF.hcalCandidates = [ cms.InputTag("hgcalBackEndLayer2Producer","HGCalBackendLayer2Processor3DClustering")]
    process.l1tPFTracksFromL1Tracks.L1TrackTag = cms.InputTag("TTTracksFromTrackletEmulation","Level1TTTracks")
    process.l1tGTTInputProducer.l1TracksInputTag = cms.InputTag("TTTracksFromTrackletEmulation","Level1TTTracks")

def isSignal():
    process.genHard = cms.EDProducer("GenParticlePruner",
        src = cms.InputTag("genParticles"),
        select = cms.vstring(
            "drop  *", # this is the default
            "keep+ abs(pdgId) == 1000022",
            "keep+ 9900001 <= abs(pdgId) <= 9900008",
        )
    )
    process.genHardParticleTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
            name = cms.string("GenHard"),
            cut  = cms.string(""),
            src = cms.InputTag("genHard"),
            doc = cms.string("Hard scattering and direct decay particles"),
            singleton = cms.bool(False), # the number of entries is variable
            extension = cms.bool(False), # this is the main table
            variables = cms.PSet(
                pt  = Var("pt",  float,precision=8),
                phi = Var("phi", float,precision=8),
                eta  = Var("eta", float,precision=8),
                mass  = Var("mass", float,precision=8),
                pdgId  = Var("pdgId", int, doc="PDG id"),
                genPartIdxMother = Var("?numberOfMothers>0?motherRef(0).key():-1", int, doc="index of the mother particle"),
                )
    )    
    process.p.replace(process.genParticleTable, process.genParticleTable + process.genHardParticleTable)
    process.extraPFStuff.add(process.genHard)
process.source.fileNames  = [ '/store/cmst3/group/l1tr/gpetrucc/12_5_X/NewInputs125X/150223/Disorder_Cascade_Phi3_M20R_PU200/inputs125X_Disorder_Cascade_Phi3_M20R_PU200_job1.root'] 
isSignal()
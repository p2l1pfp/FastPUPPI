import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras
from PhysicsTools.NanoAOD.common_cff import Var, ExtVar
import FWCore.ParameterSet.VarParsing as VarParsing
import os
import sys
import glob
from pathlib import Path


options = VarParsing.VarParsing ('analysis')

options.register ("inPath",
                  "",
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.string,
                  "Path to the input files")

options.register ("outPath",
                  ".",
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.string,
                  "Path of the output file")

options.register ("nThreads",
                  1,
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.int,
                  "Number of threads")

options.register ("nStreams",
                  0,
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.int,
                  "Number of streams")

options.register ("signal",
                  "",
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.string,
                  "Signal to consider")

options.register ("binaryDump",
                  "",
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.string,
                  "Collections to dump in binary format [puppi+tkmu+eg+pfcandall+pfcandbarrel+pfcandhgcal+trk]")

options.parseArguments()

out_dir = Path(options.outPath) / options.signal
out_dir.mkdir(parents=True, exist_ok=True)

process = cms.Process("NTUPLIZE", eras.Phase2C17I13M9)
process.load('Configuration.StandardSequences.Services_cff')
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load('Configuration.Geometry.GeometryExtendedRun4D110Reco_cff')
process.load('Configuration.Geometry.GeometryExtendedRun4D110_cff')
#process.load('Configuration.Geometry.GeometryExtended2026D110Reco_cff')
#process.load('Configuration.Geometry.GeometryExtended2026D110_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff') # needed to read HCal TPs
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('L1Trigger.Phase2L1ParticleFlow.l1tPFTracksFromL1Tracks_cfi')

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)
process.MessageLogger.cerr.FwkReport.reportEvery = 10

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        f'file:{f}' for f in glob.glob(f"{options.inPath}/*.root")
    ),
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck")
)

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '141X_mcRun4_realistic_v3', '')

process.l1tTrackSelectionProducer.processSimulatedTracks = False # these would need stubs, and are not used anyway
process.l1tPFTracksFromL1Tracks.L1TrackTag = cms.InputTag("l1tTTTracksFromTrackletEmulation","Level1TTTracks")

from L1Trigger.L1CaloTrigger.l1tPhase2L1CaloEGammaEmulator_cfi import l1tPhase2L1CaloEGammaEmulator
process.l1tPhase2L1CaloEGammaEmulator = l1tPhase2L1CaloEGammaEmulator.clone()

from L1Trigger.L1CaloTrigger.l1tPhase2CaloPFClusterEmulator_cfi import l1tPhase2CaloPFClusterEmulator
process.l1tPhase2CaloPFClusterEmulator = l1tPhase2CaloPFClusterEmulator.clone()

from L1Trigger.L1CaloTrigger.l1tPhase2GCTBarrelToCorrelatorLayer1Emulator_cfi import l1tPhase2GCTBarrelToCorrelatorLayer1Emulator
process.l1tPhase2GCTBarrelToCorrelatorLayer1Emulator = l1tPhase2GCTBarrelToCorrelatorLayer1Emulator.clone()

process.deps = cms.Task(
    process.l1tPhase2L1CaloEGammaEmulator,
    process.l1tPhase2CaloPFClusterEmulator,
    process.l1tPhase2GCTBarrelToCorrelatorLayer1Emulator,
    process.l1tTkMuonsGmt,
    process.l1tSAMuonsGmt,
    process.l1tGTTInputProducer,
    process.l1tTrackSelectionProducer,
    process.l1tVertexFinderEmulator,
    process.L1TLayer1TaskInputsTask,
    process.L1TLayer1Task,
    process.l1tLayer2EG,
    process.L1TPFJetsEmulationTask,
    # process.L1TPFJetsExtendedTask,
    # process.L1TBJetsTask,
    process.l1tPFTracksFromL1Tracks
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
        pt = Var("pt", float, precision=8),
        phi = Var("phi", float, precision=8),
        eta = Var("eta", float, precision=8),
        mass = Var("mass", float, precision=8),
        vz = Var("vz", float, precision=8),
        charge = Var("charge", int, doc="charge"),
        pdgId  = Var("pdgId", int, doc="PDG id"),
    )
)

process.pfTable = process.puppiTable.clone(
    src = cms.InputTag("l1tLayer1","PF"),
    name = cms.string("L1PF"),
    doc = cms.string("L1PF candidates"),
)

process.pfBarrelTable = process.puppiTable.clone(
    src = cms.InputTag("l1tLayer1Barrel","PF"),
    name = cms.string("L1PFBarrel"),
    doc = cms.string("L1PF candidates (barrel)"),
)

process.pfHFTable = process.puppiTable.clone(
    src = cms.InputTag("l1tLayer1HF","PF"),
    name = cms.string("L1PFHF"),
    doc = cms.string("L1PF candidates (HF)"),
)

process.pfHGCalTable = process.puppiTable.clone(
    src = cms.InputTag("l1tLayer1HGCal","PF"),
    name = cms.string("L1PFHGCal"),
    doc = cms.string("L1PF candidates (HGCal)"),
)

process.pfHGCalNoTKTable = process.puppiTable.clone(
    src = cms.InputTag("l1tLayer1HGCalNoTK","PF"),
    name = cms.string("L1PFHGCalNoTK"),
    doc = cms.string("L1PF candidates (HGCal) outside tracker"),
)

process.pfTracksTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
    src = cms.InputTag("l1tPFTracksFromL1Tracks"),
    cut = cms.string(""),
    name = cms.string("Trk"),
    doc = cms.string("L1PFTrack candidates"),
    singleton = cms.bool(False), # the number of entries is variable
    extension = cms.bool(False), # this is the main table
    variables = cms.PSet(
        pt = Var("pt", float, precision=8),
        eta = Var("eta", float, precision=8),
        phi = Var("phi", float, precision=8),
        charge = Var("charge", int, doc="charge"),
        vx = Var("vx", float, precision=8),
        vy = Var("vy", float, precision=8),
        vz = Var("vz", float, precision=8),
        # qual   = Var("quality", int, doc="quality")
    )
)

process.phoTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
    src = cms.InputTag("l1tLayer2EG", "L1CtTkEm"),
    cut = cms.string(""),
    name = cms.string("Pho"),
    doc = cms.string("Photons (TkEm) from CTL2"),
    singleton = cms.bool(False), # the number of entries is variable
    extension = cms.bool(False), # this is the main table
    variables = cms.PSet(
        pt = Var("pt", float),
        eta = Var("eta", float),
        phi = Var("phi", float),
        mass = Var("0", float),
        quality = Var("hwQual", int, doc="quality (TBD)", lazyEval=True),
        trkIsol = Var("trkIsol", float, lazyEval=True),
        trkIsolPV = Var("trkIsolPV", float, lazyEval=True),
        puppiIsol = Var("puppiIsol", float, lazyEval=True),
        puppiIsolPV = Var("puppiIsolPV", float, lazyEval=True),
        hwPt = Var("hwPt", int, lazyEval=True),
        hwEta = Var("hwEta", int, lazyEval=True),
        hwPhi = Var("hwPhi", int, lazyEval=True),
        hwQual = Var("hwQual", int, lazyEval=True),
    )
)

process.eleTable = process.phoTable.clone(
    src = "l1tLayer2EG:L1CtTkElectron",
    cut = "",
    name = "Ele",
    doc = "TkElectrons from CTL2",
    variables = dict(
        z0 = Var("trkzVtx", float, lazyEval=True),
        idScore = Var("idScore", float, lazyEval=True),
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
        pt = Var("phPt",  float, lazyEval=True),
        eta = Var("phEta", float, lazyEval=True),
        phi = Var("phPhi", float, lazyEval=True),
        mass = Var("0.10566", float),
        z0 = Var("phZ0",  float, doc="Z coordinate of the reconstructed production vertex", lazyEval=True),
        dxy = Var("phD0",  float, doc="transverse impact parameter (always zero currently)", lazyEval=True),
        charge = Var("phCharge", int, doc="charge", lazyEval=True),
        quality = Var("hwQual", int, doc="quality (TBD)", lazyEval=True),
        hwPt = Var("hwPt", int, lazyEval=True),
        hwEta = Var("hwEta", int, lazyEval=True),
        hwPhi = Var("hwPhi", int, lazyEval=True),
        hwZ0 = Var("hwZ0", int, lazyEval=True),
        hwD0 = Var("hwD0", int, lazyEval=True),
        hwCharge = Var("hwCharge", int, lazyEval=True),
        hwIsoSum = Var("hwIsoSum", int, lazyEval=True),
        hwIsoSumAp = Var("hwIsoSumAp", int, lazyEval=True),
        hwQual = Var("hwQual", int, lazyEval=True),
    )
)

process.l1VertexTable = cms.EDProducer("VertexWordFlatTableProducer",
        name = cms.string("L1Vtx"),
        cut  = cms.string(""),
        src = cms.InputTag("l1tVertexFinderEmulator","L1VerticesEmulation"),
        doc = cms.string("Primary vertices reconstructed by L1T"),
        singleton = cms.bool(False), # the number of entries is variable
        extension = cms.bool(False), # this is the main table
        variables = cms.PSet(
            sumpt = Var("pt",  float,precision=10),
            z = Var("z0",  float,precision=16),
        )
)

process.p = cms.Path(
    process.puppiTable +
    process.pfTable +
    process.pfBarrelTable +
    process.pfHFTable +
    process.pfHGCalTable +
    process.pfHGCalNoTKTable +
    process.pfTracksTable +
    process.phoTable +
    process.eleTable +
    process.tkMuTable +
    process.l1VertexTable
)

if "puppi" in options.binaryDump:
    process.puppiDump = cms.EDAnalyzer("L1PuppiBinaryDumper",
                    src = cms.InputTag("l1tLayer2Deregionizer:Puppi"),
                    outName = cms.string(f"{str(out_dir)}/{options.signal}_puppi.dump"))
    process.p += process.puppiDump

if "pfcandall" in options.binaryDump:
    process.pfDump = cms.EDAnalyzer("L1PuppiBinaryDumper",
                    src = cms.InputTag("l1tLayer1","PF"),
                    outName = cms.string(f"{str(out_dir)}/{options.signal}_pfcandall.dump"))
    process.p += process.pfDump

if "pfcandbarrel" in options.binaryDump:
    process.pfBarrelDump = cms.EDAnalyzer("L1PuppiBinaryDumper",
                    src = cms.InputTag("l1tLayer1Barrel","PF"),
                    outName = cms.string(f"{str(out_dir)}/{options.signal}_pfcandbarrel.dump"))
    process.p += process.pfBarrelDump

if "pfcandhgcal" in options.binaryDump:
    process.pfHGCalDump = cms.EDAnalyzer("L1PuppiBinaryDumper",
                    src = cms.InputTag("l1tLayer1HGCal","PF"),
                    outName = cms.string(f"{str(out_dir)}/{options.signal}_pfcandhgcal.dump"))
    process.p += process.pfHGCalDump

if "tkmu" in options.binaryDump:
    process.tkMuDump = cms.EDAnalyzer("L1TrackerMuonBinaryDumper",
                    src = cms.InputTag("l1tTkMuonsGmt"),
                    outName = cms.string(f"{str(out_dir)}/{options.signal}_tkMuons.dump"))
    process.p += process.tkMuDump

if "eg" in options.binaryDump:
    process.egDump = cms.EDAnalyzer("L1CTL2EgammaBinaryDumper",
                    srcEle = cms.InputTag("l1tLayer2EG", "L1CtTkElectron"),
                    srcEm = cms.InputTag("l1tLayer2EG", "L1CtTkEm"),
                    interleaveOutputs = cms.bool(False), # False = first 12 photons, then electrons; True = pho1, ele1, pho2, ele2, ...
                    outName = cms.string(f"{str(out_dir)}/{options.signal}_egamma.dump"))
    process.p += process.egDump

if "trk" in options.binaryDump:
    process.trkDump = cms.EDAnalyzer("L1TrackerTrackBinaryDumper",
                    src = cms.InputTag("l1tPFTracksFromL1Tracks"),
                    outName = cms.string(f"{str(out_dir)}/{options.signal}_trk.dump"))
    process.p += process.trkDump

process.p.associate(process.deps)

def genXToQGamma(PQ="", PD=""):
    if PQ=="Phi":
        PQid = 333
        PDid = 321
    elif PQ=="Rho":
        PQid = 113
        PDid = 211
    elif PQ=="JPsi":
        PQid = 443
        if PD=="Mu":
            PDid = 13
        elif PD=="El":
            PDid = 11

    genH_cut = f"abs(pdgId)==25 && numberOfDaughters==2 && ((abs(daughter(0).pdgId)==22 && abs(daughter(1).pdgId)=={PQid}) || (abs(daughter(0).pdgId)=={PQid} && abs(daughter(1).pdgId)==22))"
    genQFromH_cut = f"abs(pdgId)=={PQid} && numberOfMothers>0 && abs(motherRef.pdgId)==25"
    genQ_cut = f"abs(pdgId)=={PQid} && numberOfDaughters==2 && abs(daughter(0).pdgId)=={PDid} && abs(daughter(1).pdgId)=={PDid}"
    genPDFromPQ_cut = f"abs(pdgId)=={PDid} && numberOfMothers>0 && abs(motherRef.pdgId)=={PQid} && abs(motherRef.motherRef.pdgId)==25"
    genGammaFromH_cut = f"abs(pdgId)==22 && numberOfMothers>0 && abs(motherRef.pdgId)==25"

    process.genH = cms.EDFilter("GenParticleSelector",
        src = cms.InputTag("genParticles"),
        cut = cms.string(genH_cut),
        filter = cms.bool(True),
    )
    process.genPQFromH = cms.EDFilter("GenParticleSelector",
        src = cms.InputTag("genParticles"),
        cut = cms.string(genQFromH_cut),
    )
    process.genPQ = cms.EDFilter("GenParticleSelector",
        src = cms.InputTag("genParticles"),
        cut = cms.string(genQ_cut),
    )
    process.genPDFromPQ = cms.EDFilter("GenParticleSelector",
        src = cms.InputTag("genParticles"),
        cut = cms.string(genPDFromPQ_cut),
    )
    process.genGammaFromH = cms.EDFilter("GenParticleSelector",
        src = cms.InputTag("genParticles"),
        cut = cms.string(genGammaFromH_cut),
    )
    process.genHTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
        src = cms.InputTag("genH"),
        cut = cms.string(""),
        name = cms.string("GenH"),
        doc = cms.string("gen H boson"),
        singleton = cms.bool(False), # the number of entries is variable
        extension = cms.bool(False), # this is the main table
        variables = cms.PSet(
            pt     = Var("pt",  float,precision=8),
            phi    = Var("phi", float,precision=8),
            eta    = Var("eta", float,precision=8),
            mass   = Var("mass", float,precision=8),
            vz     = Var("vz", float,precision=8, doc="Production point along the beam axis"),
            charge = Var("charge", int, doc="charge"),
            pdgId  = Var("pdgId", int, doc="PDG id"),
        )
    )
    process.genPQTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
        src = cms.InputTag("genPQFromH"),
        cut = cms.string(""),
        name = cms.string("GenPQ"),
        doc = cms.string("gen Q from H boson decay"),
        singleton = cms.bool(False), # the number of entries is variable
        extension = cms.bool(False), # this is the main table
        variables = cms.PSet(
            pt     = Var("pt",  float,precision=8),
            phi    = Var("phi", float,precision=8),
            eta    = Var("eta", float,precision=8),
            mass   = Var("mass", float,precision=8),
            vz     = Var("vz", float,precision=8, doc="Production point along the beam axis"),
            charge = Var("charge", int, doc="charge"),
            pdgId  = Var("pdgId", int, doc="PDG id"),
            # prompt = Var("2*statusFlags().isPrompt() + statusFlags().isDirectPromptTauDecayProduct()", int, doc="Particle status."),
        )
    )
    process.genPDTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
        src = cms.InputTag("genPDFromPQ"),
        cut = cms.string(""),
        name = cms.string("GenPD"),
        doc = cms.string("gen daughters from Q meson decay"),
        singleton = cms.bool(False), # the number of entries is variable
        extension = cms.bool(False), # this is the main table
        variables = cms.PSet(
            pt     = Var("pt",  float,precision=8),
            phi    = Var("phi", float,precision=8),
            eta    = Var("eta", float,precision=8),
            mass   = Var("mass", float,precision=8),
            vz     = Var("vz", float,precision=8, doc="Production point along the beam axis"),
            charge = Var("charge", int, doc="charge"),
            pdgId  = Var("pdgId", int, doc="PDG id"),
            # prompt = Var("2*statusFlags().isPrompt() + statusFlags().isDirectPromptTauDecayProduct()", int, doc="Particle status."),
        )
    )
    process.genGammaTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
        src = cms.InputTag("genGammaFromH"),
        cut = cms.string(""),
        name = cms.string("GenGamma"),
        doc = cms.string("gen photon from Higgs decay"),
        singleton = cms.bool(False), # the number of entries is variable
        extension = cms.bool(False), # this is the main table
        variables = cms.PSet(
            pt     = Var("pt",  float,precision=8),
            phi    = Var("phi", float,precision=8),
            eta    = Var("eta", float,precision=8),
            mass   = Var("mass", float,precision=8),
            vz     = Var("vz",  float,precision=8, doc="Production point along the beam axis"),
            charge = Var("charge", int, doc="charge"),
            pdgId  = Var("pdgId", int, doc="PDG id"),
            # prompt = Var("2*statusFlags().isPrompt() + statusFlags().isDirectPromptTauDecayProduct()", int, doc="Particle status."),
        )
    )
    process.deps.add(
        process.genH,
        process.genPQFromH,
        process.genGammaFromH,
        process.genPQ,
        process.genPDFromPQ,
    )
    process.p.replace(process.puppiTable, process.puppiTable + process.genHTable + process.genPQTable + process.genPDTable + process.genGammaTable)


def genXToQQ(PQ1="", PQ2="", PD12="", PD34=""):
    if PQ1=="Phi":
        PQ1id = 333
        PD12id = 321
    elif PQ1=="Rho":
        PQ1id = 113
        PD12id = 211
    elif PQ1=="JPsi":
        PQ1id = 443
        if PD12=="Mu":
            PD12id = 13
        elif PD12=="El":
            PD12id = 11

    if PQ2=="Phi":
        PQ2id = 333
        PD34id = 321
    elif PQ2=="Rho":
        PQ2id = 113
        PD34id = 211
    elif PQ2=="JPsi":
        PQ2id = 443
        if PD34=="Mu":
            PD34id = 13
        elif PD34=="El":
            PD34id = 11

    if PQ1==PQ2:
        genH_cut = f"abs(pdgId)==25 && numberOfDaughters==2 && (abs(daughter(0).pdgId)=={PQ1id} && abs(daughter(1).pdgId)=={PQ2id})"
        genQFromH_cut = f"abs(pdgId)=={PQ1id} && numberOfMothers>0 && abs(motherRef.pdgId)==25"
        genQ_cut = f"abs(pdgId)=={PQ1id} && numberOfDaughters==2 && abs(daughter(0).pdgId)=={PD12id} && abs(daughter(1).pdgId)=={PD12id}"
        genPDFromPQ_cut = f"abs(pdgId)=={PD12id} && numberOfMothers>0 && abs(motherRef.pdgId)=={PQ1id} && abs(motherRef.motherRef.pdgId)==25"
    else:
        genH_cut = f"abs(pdgId)==25 && numberOfDaughters==2 && ((abs(daughter(0).pdgId)=={PQ1id} && abs(daughter(1).pdgId)=={PQ2id}) || (abs(daughter(0).pdgId)=={PQ2id} && abs(daughter(1).pdgId)=={PQ1id}))"
        genQ1FromH_cut = f"abs(pdgId)=={PQ1id} && numberOfMothers>0 && abs(motherRef.pdgId)==25"
        genQ2FromH_cut = f"abs(pdgId)=={PQ2id} && numberOfMothers>0 && abs(motherRef.pdgId)==25"
        genQ_cut = f"abs(pdgId)=={PQ1id} && numberOfDaughters==2 && abs(daughter(0).pdgId)=={PD12id} && abs(daughter(1).pdgId)=={PD12id}"
        genPD12FromPQ1_cut = f"abs(pdgId)=={PD12id} && numberOfMothers>0 && abs(motherRef.pdgId)=={PQ1id} && abs(motherRef.motherRef.pdgId)==25"
        genPD34FromPQ2_cut = f"abs(pdgId)=={PD34id} && numberOfMothers>0 && abs(motherRef.pdgId)=={PQ2id} && abs(motherRef.motherRef.pdgId)==25"

    process.genH = cms.EDFilter("GenParticleSelector",
        src = cms.InputTag("genParticles"),
        cut = cms.string(genH_cut),
        filter = cms.bool(True),
    )
    if PQ1==PQ2:
        process.genPQFromH = cms.EDFilter("GenParticleSelector",
            src = cms.InputTag("genParticles"),
            cut = cms.string(genQFromH_cut),
        )
        process.genPQ = cms.EDFilter("GenParticleSelector",
            src = cms.InputTag("genParticles"),
            cut = cms.string(genQ_cut),
        )
        process.genPDFromPQ = cms.EDFilter("GenParticleSelector",
            src = cms.InputTag("genParticles"),
            cut = cms.string(genPDFromPQ_cut),
        )
    else:
        process.genPQ1FromH = cms.EDFilter("GenParticleSelector",
            src = cms.InputTag("genParticles"),
            cut = cms.string(genQ1FromH_cut),
        )
        process.genPQ2FromH = cms.EDFilter("GenParticleSelector",
            src = cms.InputTag("genParticles"),
            cut = cms.string(genQ2FromH_cut),
        )
        process.genPD12FromPQ1 = cms.EDFilter("GenParticleSelector",
            src = cms.InputTag("genParticles"),
            cut = cms.string(genPD12FromPQ1_cut),
        )
        process.genPD34FromPQ2 = cms.EDFilter("GenParticleSelector",
            src = cms.InputTag("genParticles"),
            cut = cms.string(genPD34FromPQ2_cut),
        )


    process.genHTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
        src = cms.InputTag("genH"),
        cut = cms.string(""),
        name = cms.string("GenH"),
        doc = cms.string("gen H boson"),
        singleton = cms.bool(False), # the number of entries is variable
        extension = cms.bool(False), # this is the main table
        variables = cms.PSet(
            pt     = Var("pt",  float,precision=8),
            phi    = Var("phi", float,precision=8),
            eta    = Var("eta", float,precision=8),
            mass   = Var("mass", float,precision=8),
            vz     = Var("vz", float,precision=8, doc="Production point along the beam axis"),
            charge = Var("charge", int, doc="charge"),
            pdgId  = Var("pdgId", int, doc="PDG id"),
        )
    )

    if PQ1==PQ2:
        process.genPQTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
            src = cms.InputTag("genPQFromH"),
            cut = cms.string(""),
            name = cms.string(f"Gen{PQ1}"),
            doc = cms.string(f"gen {PQ1} from H boson decay"),
            singleton = cms.bool(False), # the number of entries is variable
            extension = cms.bool(False), # this is the main table
            variables = cms.PSet(
                pt     = Var("pt",  float,precision=8),
                phi    = Var("phi", float,precision=8),
                eta    = Var("eta", float,precision=8),
                mass   = Var("mass", float,precision=8),
                vz     = Var("vz", float,precision=8, doc="Production point along the beam axis"),
                charge = Var("charge", int, doc="charge"),
                pdgId  = Var("pdgId", int, doc="PDG id"),
                # prompt = Var("2*statusFlags().isPrompt() + statusFlags().isDirectPromptTauDecayProduct()", int, doc="Particle status."),
            )
        )
        process.genPDTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
            src = cms.InputTag("genPDFromPQ"),
            cut = cms.string(""),
            name = cms.string(f"Gen{PD12}"),
            doc = cms.string(f"gen daughters from {PQ1} meson decay"),
            singleton = cms.bool(False), # the number of entries is variable
            extension = cms.bool(False), # this is the main table
            variables = cms.PSet(
                pt     = Var("pt",  float,precision=8),
                phi    = Var("phi", float,precision=8),
                eta    = Var("eta", float,precision=8),
                mass   = Var("mass", float,precision=8),
                vz     = Var("vz", float,precision=8, doc="Production point along the beam axis"),
                charge = Var("charge", int, doc="charge"),
                pdgId  = Var("pdgId", int, doc="PDG id"),
                # prompt = Var("2*statusFlags().isPrompt() + statusFlags().isDirectPromptTauDecayProduct()", int, doc="Particle status."),
            )
        )
    else:
        process.genPQ1Table = cms.EDProducer("SimpleCandidateFlatTableProducer",
            src = cms.InputTag("genPQ1FromH"),
            cut = cms.string(""),
            name = cms.string(f"Gen{PQ1}"),
            doc = cms.string(f"gen {PQ1} from H boson decay"),
            singleton = cms.bool(False), # the number of entries is variable
            extension = cms.bool(False), # this is the main table
            variables = cms.PSet(
                pt     = Var("pt",  float,precision=8),
                phi    = Var("phi", float,precision=8),
                eta    = Var("eta", float,precision=8),
                mass   = Var("mass", float,precision=8),
                vz     = Var("vz", float,precision=8, doc="Production point along the beam axis"),
                charge = Var("charge", int, doc="charge"),
                pdgId  = Var("pdgId", int, doc="PDG id"),
                # prompt = Var("2*statusFlags().isPrompt() + statusFlags().isDirectPromptTauDecayProduct()", int, doc="Particle status."),
            )
        )
        process.genPQ2Table = cms.EDProducer("SimpleCandidateFlatTableProducer",
            src = cms.InputTag("genPQ2FromH"),
            cut = cms.string(""),
            name = cms.string(f"Gen{PQ2}"),
            doc = cms.string(f"gen {PQ2} from H boson decay"),
            singleton = cms.bool(False), # the number of entries is variable
            extension = cms.bool(False), # this is the main table
            variables = cms.PSet(
                pt     = Var("pt",  float,precision=8),
                phi    = Var("phi", float,precision=8),
                eta    = Var("eta", float,precision=8),
                mass   = Var("mass", float,precision=8),
                vz     = Var("vz", float,precision=8, doc="Production point along the beam axis"),
                charge = Var("charge", int, doc="charge"),
                pdgId  = Var("pdgId", int, doc="PDG id"),
                # prompt = Var("2*statusFlags().isPrompt() + statusFlags().isDirectPromptTauDecayProduct()", int, doc="Particle status."),
            )
        )
        process.genPD12Table = cms.EDProducer("SimpleCandidateFlatTableProducer",
            src = cms.InputTag("genPD12FromPQ1"),
            cut = cms.string(""),
            name = cms.string(f"Gen{PD12}"),
            doc = cms.string(f"gen daughters from {PQ1} meson decay"),
            singleton = cms.bool(False), # the number of entries is variable
            extension = cms.bool(False), # this is the main table
            variables = cms.PSet(
                pt     = Var("pt",  float,precision=8),
                phi    = Var("phi", float,precision=8),
                eta    = Var("eta", float,precision=8),
                mass   = Var("mass", float,precision=8),
                vz     = Var("vz", float,precision=8, doc="Production point along the beam axis"),
                charge = Var("charge", int, doc="charge"),
                pdgId  = Var("pdgId", int, doc="PDG id"),
                # prompt = Var("2*statusFlags().isPrompt() + statusFlags().isDirectPromptTauDecayProduct()", int, doc="Particle status."),
            )
        )
        process.genPD34able = cms.EDProducer("SimpleCandidateFlatTableProducer",
            src = cms.InputTag("genPD34FromPQ2"),
            cut = cms.string(""),
            name = cms.string(f"Gen{PD34}"),
            doc = cms.string(f"gen daughters from {PQ2} meson decay"),
            singleton = cms.bool(False), # the number of entries is variable
            extension = cms.bool(False), # this is the main table
            variables = cms.PSet(
                pt     = Var("pt",  float,precision=8),
                phi    = Var("phi", float,precision=8),
                eta    = Var("eta", float,precision=8),
                mass   = Var("mass", float,precision=8),
                vz     = Var("vz", float,precision=8, doc="Production point along the beam axis"),
                charge = Var("charge", int, doc="charge"),
                pdgId  = Var("pdgId", int, doc="PDG id"),
                # prompt = Var("2*statusFlags().isPrompt() + statusFlags().isDirectPromptTauDecayProduct()", int, doc="Particle status."),
            )
        )

    if PQ1==PQ2:
        process.deps.add(
            process.genH,
            process.genPQFromH,
            process.genPQ,
            process.genPDFromPQ,
        )
        process.p.replace(process.puppiTable, process.puppiTable + process.genHTable + process.genPQTable + process.genPDTable)
    else:
        process.deps.add(
            process.genH,
            process.genPQ1FromH,
            process.genPQ2FromH,
            process.genPQ1,
            process.genPQ2,
            process.genPD12FromPQ,
            process.genPD34FromPQ,
        )
        process.p.replace(
            process.puppiTable,
            process.puppiTable + process.genHTable + process.genPQ1Table + process.genPQ2Table + process.genPD12Table + process.genPD34Table
        )





process.outnano = cms.OutputModule("NanoAODOutputModule",
    fileName = cms.untracked.string(f"{str(out_dir)}/{options.signal}_L1NANO.root"),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('p')),
    outputCommands = cms.untracked.vstring(
        "drop *",
        "keep nanoaodFlatTable_*Table_*_*"
    ),
    compressionLevel = cms.untracked.int32(4),
    compressionAlgorithm = cms.untracked.string("ZLIB"),
)
process.end = cms.EndPath(process.outnano)
process.schedule = cms.Schedule(process.p, process.end)

process.options.numberOfThreads = options.nThreads
process.options.numberOfStreams = options.nStreams

if options.signal[1:]=="ToRhoGammaTo2PiGamma":
    genXToQGamma(PQ="Rho", PD="Pi")
elif options.signal[1:]=="ToPhiGammaTo2KGamma":
    genXToQGamma(PQ="Phi", PD="K")
elif options.signal[1:]=="ToJPsiGammaTo2MuGamma":
    genXToQGamma(PQ="JPsi", PD="Mu")
elif options.signal[1:]=="ToJPsiGammaTo2ElGamma":
    genXToQGamma(PQ="JPsi", PD="El")
elif options.signal[1:]=="To2PhiTo4K":
    genXToQQ(PQ1="Phi", PQ2="Phi", PD12="K", PD34="K")
elif "SingleNeutrino" in options.signal:
    pass
elif options.signal=="TTBar":
    pass
else:
    print(f"Wrong signal chosen: {options.signal}")
    exit()

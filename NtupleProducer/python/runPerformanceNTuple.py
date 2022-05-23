import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras
from PhysicsTools.NanoAOD.common_cff import Var, ExtVar 

process = cms.Process("RESP", eras.Phase2C9)

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
            "drop l1tPFCandidates_*_*_*")
)

process.load('Configuration.Geometry.GeometryExtended2026D49Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2026D49_cff')
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
process.GlobalTag = GlobalTag(process.GlobalTag, '123X_mcRun4_realistic_v3', '')


process.extraPFStuff = cms.Task(
        process.L1SAMuonsGmt,
        process.L1GTTInputProducer,
        process.L1VertexFinderEmulator,
        process.l1ctLayer1TaskInputsTask,
        process.l1ctLayer1Task)

process.centralGen = cms.EDFilter("CandPtrSelector", src = cms.InputTag("genParticlesForMETAllVisible"), cut = cms.string("abs(eta) < 2.4"))
process.barrelGen = cms.EDFilter("CandPtrSelector", src = cms.InputTag("genParticlesForMETAllVisible"), cut = cms.string("abs(eta) < 1.5"))
process.genMetCentralTrue = process.genMetTrue.clone(src = cms.InputTag("centralGen"))
process.extraPFStuff.add(
    process.genParticlesForMETAllVisible,
    process.centralGen,
    process.barrelGen,
    process.genMetCentralTrue
)

def monitorPerf(label, tag, makeResp=True, makeRespSplit=True, makeJets=True, makeMET=True, makeCentralMET=True,
                makeInputMultiplicities=False, makeOutputMultiplicities=False, saveCands=False):
    def _add(name, what):
        setattr(process, name, what)
        process.extraPFStuff.add(what)
    if type(tag) != str and len(tag) > 1:
        _add('merged'+label, cms.EDProducer("L1TPFCandMerger", src = cms.VInputTag(cms.InputTag(x) for x in tag)))
        tag = 'merged'+label
    if makeResp:
        setattr(process.ntuple.objects, label, cms.VInputTag(cms.InputTag(tag)))
        if makeRespSplit:
            setattr(process.ntuple.objects, label+"Charged", cms.VInputTag(cms.InputTag(tag)))
            setattr(process.ntuple.objects, label+"Charged_sel", cms.string("charge != 0"))
            setattr(process.ntuple.objects, label+"Photon",  cms.VInputTag(cms.InputTag(tag)))
            setattr(process.ntuple.objects, label+"Photon_sel", cms.string("pdgId == 22"))
            setattr(process.ntuple.objects, label+"Neutral",  cms.VInputTag(cms.InputTag(tag)))
            setattr(process.ntuple.objects, label+"Neutral_sel", cms.string("charge == 0"))
            setattr(process.ntuple.objects, label+"NeutralHad",  cms.VInputTag(cms.InputTag(tag)))
            setattr(process.ntuple.objects, label+"NeutralHad_sel", cms.string("charge == 0 && pdgId != 22"))
    if makeJets:
        _add('ak4'+label, ak4PFJets.clone(src = tag, doAreaFastjet = False))
        setattr(process.l1pfjetTable.jets, label, cms.InputTag('ak4'+label))
    if saveCands:
        setattr(process.l1pfcandTable.cands, label, cms.InputTag(tag))
    if makeMET:
        _add('met'+label, pfMet.clone(src = tag, calculateSignificance = False))
        setattr(process.l1pfmetTable.mets, label, cms.InputTag('met'+label))
        if makeCentralMET:
            _add('central'+label, cms.EDFilter("CandPtrSelector", src = cms.InputTag(tag), cut = cms.string("abs(eta) < 2.4")))
            _add('met'+label+'Central', pfMet.clone(src = 'central'+label, calculateSignificance = False))
            setattr(process.l1pfmetCentralTable.mets, label, cms.InputTag('met'+label+'Central'))
    if makeInputMultiplicities == "CTL1":
        D = tag.split(":")[0] # l1ctLayer1[Barrel,HGCal,HF] usually
        I = tag.split(":")[1] # Calo, EmCalo, TK, or Mu, usually
        for W in ("Reg", "Sec"):
            for X in ("tot","max"):
                process.ntuple.copyUInts.append( "%s:%sN%s%s" % (D,X,W,I))
            process.ntuple.copyVecUInts.append( "%s:vecN%s%s" % (D,W,I))
    if makeOutputMultiplicities == "CTL1":
        D = tag.split(":")[0] # l1pfProducer[Barrel,HGCal,HF] usually
        P = tag.split(":")[1] # PF or Puppi, usually
        for O in ["", "Charged", "Neutral", "Electron", "Muon", "ChargedHadron", "NeutralHadron", "Photon"]:
            for X in ("tot","max"):
                process.ntuple.copyUInts.append( "%s:%sN%s%s" % (D,X,P,O))
            process.ntuple.copyVecUInts.append( "%s:vecN%s%s" % (D,P,O))    

process.ntuple = cms.EDAnalyzer("ResponseNTuplizer",
    genJets = cms.InputTag("ak4GenJetsNoNu"),
    genParticles = cms.InputTag("genParticles"),
    isParticleGun = cms.bool(False),
    writeExtraInfo = cms.bool(False),
    doRandom = cms.bool(False),
    objects = cms.PSet(
        # -- inputs and PF --
        RawTK  = cms.VInputTag('pfTracksFromL1Tracks',),
        # outputs
    ),
    copyUInts = cms.VInputTag(),
    copyFloats = cms.VInputTag(),
    copyVecUInts = cms.VInputTag(),
)

process.extraPFStuff.add(process.pfTracksFromL1Tracks)

process.l1pfjetTable = cms.EDProducer("L1PFJetTableProducer",
    gen = cms.InputTag("ak4GenJetsNoNu"),
    commonSel = cms.string("pt > 5 && abs(eta) < 5.0"),
    drMax = cms.double(0.2),
    minRecoPtOverGenPt = cms.double(0.1),
    jets = cms.PSet(
        Gen = cms.InputTag("ak4GenJetsNoNu"),
        Gen_sel = cms.string("pt > 15"),
    ),
    moreVariables = cms.PSet(
    ),
)

process.l1pfmetTable = cms.EDProducer("L1PFMetTableProducer",
    genMet = cms.InputTag("genMetTrue"), 
    flavour = cms.string(""),
    mets = cms.PSet(
    ),
)
process.l1pfmetCentralTable = process.l1pfmetTable.clone(genMet = "genMetCentralTrue", flavour = "Central")

monitorPerf("L1Calo", "l1ctLayer1:Calo")
monitorPerf("L1TK",   "l1ctLayer1:TK")
monitorPerf("L1PF",    "l1ctLayer1:PF")
monitorPerf("L1Puppi", "l1ctLayer1:Puppi")

# to check available tags:
#process.content = cms.EDAnalyzer("EventContentAnalyzer")
process.p = cms.Path(
        process.ntuple + #process.content +
        process.l1pfjetTable + 
        process.l1pfmetTable + process.l1pfmetCentralTable
        )
process.p.associate(process.extraPFStuff)
process.TFileService = cms.Service("TFileService", fileName = cms.string("perfTuple.root"))

# for full debug:
#process.out = cms.OutputModule("PoolOutputModule",
#                               fileName = cms.untracked.string("debugPF.root"),
#                               SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring("p"))
#                           )
#process.end = cms.EndPath(process.out)

process.outnano = cms.OutputModule("NanoAODOutputModule",
    fileName = cms.untracked.string("perfNano.root"),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('p')),
    outputCommands = cms.untracked.vstring("drop *", "keep nanoaodFlatTable_*Table_*_*"),
    compressionLevel = cms.untracked.int32(4),
    compressionAlgorithm = cms.untracked.string("ZLIB"),
)
process.end = cms.EndPath(process.outnano)

# Below for more debugging
if True:
    process.genInAcceptance = cms.EDFilter("GenParticleSelector",
        src = cms.InputTag("genParticles"),
        cut = cms.string("status == 1 && (abs(pdgId) != 12 && abs(pdgId) != 14 && abs(pdgId) != 16) && "+
                         "(abs(eta) < 2.5 && pt > 2 && charge != 0 || "+
                         "abs(pdgId) == 22 && pt > 1 || "+
                         "charge == 0 && pt > 1 || "+
                         "charge != 0 && abs(eta) > 2.5 && pt > 2) ") # tracks below pT 2 bend by more than 0.4,
    )
    process.ntuple.objects.GenAcc = cms.VInputTag(cms.InputTag("genInAcceptance"))
    process.ntuple.objects.ChGenAcc = cms.VInputTag(cms.InputTag("genInAcceptance"))
    process.ntuple.objects.ChGenAcc_sel = cms.string("(abs(eta) < 2.5 && pt > 2 && charge != 0)")
    process.ntuple.objects.PhGenAcc = cms.VInputTag(cms.InputTag("genInAcceptance"))
    process.ntuple.objects.PhGenAcc_sel = cms.string("pdgId == 22")
    process.extraPFStuff.add(process.genInAcceptance)
def respOnly():
    process.p.remove(process.l1pfjetTable)
    process.p.remove(process.l1pfmetTable)
    process.p.remove(process.l1pfmetCentralTable)
    process.end.remove(process.outnano)
def noResp():
    process.p.remove(process.ntuple)

def addMult():
    for D in ['Barrel','HF','HGCal','HGCalNoTK']:
        monitorPerf("L1%sCalo"%D,  "l1ctLayer1%s:Calo"%D,   makeResp=False, makeRespSplit=False, makeJets=False, makeMET=False, makeCentralMET=False, makeInputMultiplicities="CTL1")
        monitorPerf("L1%sEmCalo"%D,"l1ctLayer1%s:EmCalo"%D, makeResp=False, makeRespSplit=False, makeJets=False, makeMET=False, makeCentralMET=False, makeInputMultiplicities="CTL1")
        monitorPerf("L1%sTK"%D,    "l1ctLayer1%s:TK"%D,     makeResp=False, makeRespSplit=False, makeJets=False, makeMET=False, makeCentralMET=False, makeInputMultiplicities="CTL1")
        monitorPerf("L1%sMu"%D,    "l1ctLayer1%s:Mu"%D,     makeResp=False, makeRespSplit=False, makeJets=False, makeMET=False, makeCentralMET=False, makeInputMultiplicities="CTL1")
        monitorPerf("L1%sPF"%D,    "l1ctLayer1%s:PF"%D,     makeResp=False, makeRespSplit=False, makeJets=False, makeMET=False, makeCentralMET=False, makeOutputMultiplicities="CTL1")
        monitorPerf("L1%sPuppi"%D, "l1ctLayer1%s:Puppi"%D,  makeResp=False, makeRespSplit=False, makeJets=False, makeMET=False, makeCentralMET=False, makeOutputMultiplicities="CTL1")


def addCHS():
    process.l1PuppiCharged = cms.EDFilter("L1TPFCandSelector",
        src = cms.InputTag("l1ctLayer1:Puppi"),
        cut = cms.string("charge != 0"))
    process.l1PFNeutral = cms.EDFilter("L1TPFCandSelector",
        src = cms.InputTag("l1ctLayer1:PF"),
        cut = cms.string("charge == 0"))
    process.extraPFStuff.add(process.l1PuppiCharged, process.l1PFNeutral)
    monitorPerf("L1CHS", [ "l1PuppiCharged", "l1PFNeutral" ], makeRespSplit = False)
def addCalib():
    process.load("L1Trigger.Phase2L1ParticleFlow.pfClustersFromHGC3DClustersEM_cfi")
    process.pfClustersFromL1EGClustersRaw    = process.pfClustersFromL1EGClusters.clone(corrector = "")
    process.pfClustersFromHGC3DClustersRaw   = process.pfClustersFromHGC3DClusters.clone(corrector = "")
    process.pfClustersFromHGC3DClustersEMRaw = process.pfClustersFromHGC3DClustersEM.clone(corrector = "")
    process.extraPFStuff.add(
            process.pfClustersFromL1EGClustersRaw, 
            process.pfClustersFromHGC3DClustersRaw, 
            process.pfClustersFromHGC3DClustersEM,
            process.pfClustersFromHGC3DClustersEMRaw)
    process.ntuple.objects.L1RawBarrelEcal   = cms.VInputTag('pfClustersFromL1EGClustersRaw' )
    process.ntuple.objects.L1RawBarrelCalo   = cms.VInputTag('pfClustersFromCombinedCaloHCal:uncalibrated')
    process.ntuple.objects.L1RawBarrelCaloEM = cms.VInputTag('pfClustersFromCombinedCaloHCal:emUncalibrated')
    process.ntuple.objects.L1RawHGCal   = cms.VInputTag('pfClustersFromHGC3DClustersRaw')
    process.ntuple.objects.L1RawHGCalEM = cms.VInputTag('pfClustersFromHGC3DClustersEMRaw')
    process.ntuple.objects.L1RawHFCalo  = cms.VInputTag('pfClustersFromCombinedCaloHF:uncalibrated')
    process.ntuple.objects.L1BarrelEcal = cms.VInputTag('pfClustersFromL1EGClusters' )
    process.ntuple.objects.L1BarrelCalo = cms.VInputTag('pfClustersFromCombinedCaloHCal:calibrated')
    process.ntuple.objects.L1HGCal   = cms.VInputTag('pfClustersFromHGC3DClusters')
    process.ntuple.objects.L1HFCalo  = cms.VInputTag('pfClustersFromCombinedCaloHF:calibrated')
    process.ntuple.objects.L1HGCalEM = cms.VInputTag('pfClustersFromHGC3DClustersEM', )

def addSeededConeJets(what="Puppi",src="l1ctLayer1:Puppi"):
    process.load('L1Trigger.Phase2L1ParticleFlow.L1SeedConePFJetProducer_cfi')
    scModule = process.L1SeedConePFJetProducer.clone(L1PFObjects = src)
    scEmuModule = process.L1SeedConePFJetEmulatorProducer.clone(L1PFObjects = src)
    setattr(process, 'sc'+what, scModule)
    setattr(process, 'scEmu'+what, scEmuModule)
    process.extraPFStuff.add(scModule, scEmuModule)
    setattr(process.l1pfjetTable.jets, 'sc'+what, cms.InputTag('sc'+what))
    setattr(process.l1pfjetTable.jets, 'scEmu'+what, cms.InputTag('scEmu'+what))
def addDeregSeededConeJets():
    process.load('L1Trigger.Phase2L1ParticleFlow.L1SeedConePFJetProducer_cfi')
    process.load('L1Trigger.Phase2L1ParticleFlow.DeregionizerProducer_cfi')
    process.scDeregPuppiJets = process.L1SeedConePFJetEmulatorProducer.clone(L1PFObjects = "DeregionizerProducer:Puppi")
    process.extraPFStuff.add(process.scDeregPuppiJets, process.DeregionizerProducer)
    setattr(process.l1pfjetTable.jets, 'scDeregPuppi', cms.InputTag('scDeregPuppiJets'))



def addRefs(calo=True,tk=True):
    process.load('L1Trigger.L1CaloTrigger.Phase1L1TJets_cff')
    process.Phase1L1TJetProducer.inputCollectionTag = cms.InputTag("l1ctLayer1", "Puppi") # make sure the process name is not pre-encoded
    process.extraPFStuff.add(process.Phase1L1TJetProducer,process.Phase1L1TJetCalibrator)
    process.l1pfjetTable.jets.RefPhase1PuppiJets = cms.InputTag("Phase1L1TJetCalibrator", "Phase1L1TJetFromPfCandidates")
    if calo:
        process.l1pfjetTable.jets.RefCaloJets = cms.InputTag("L1CaloJetProducer","L1CaloJetCollectionBXV")
    if tk:
        process.l1pfjetTable.jets.RefTwoLayerJets = cms.InputTag("TwoLayerJets", "L1TwoLayerJets")
        process.l1pfjetTable.jets.RefTwoLayerJets_sel = cms.string("pt > 5")
        process.l1pfmetTable.mets.RefL1TrackerEtMiss = cms.InputTag("L1TrackerEtMiss","trkMET")
        process.l1pfmetCentralTable.mets.RefL1TrackerEtMiss = cms.InputTag("L1TrackerEtMiss","trkMET")

def addTkPtCut(ptCut):
    process.l1ctLayer1BarrelTkPt3 = process.l1ctLayer1Barrel.clone(trkPtCut = ptCut)
    process.l1ctLayer1HGCalTkPt3 = process.l1ctLayer1HGCal.clone(trkPtCut = ptCut)
    process.l1ctLayer1TkPt3 = cms.EDProducer("L1TPFCandMultiMerger",
        pfProducers = cms.VInputTag(
            cms.InputTag("l1ctLayer1BarrelTkPt3"), 
            cms.InputTag("l1ctLayer1HGCalTkPt3"),
            cms.InputTag("l1ctLayer1HGCalNoTK"),
            cms.InputTag("l1ctLayer1HF")
            ),
        labelsToMerge = cms.vstring("TK", "PF", "Puppi"),
        regionalLabelsToMerge = cms.vstring(),
    )
    process.extraPFStuff.add(process.l1ctLayer1BarrelTkPt3, process.l1ctLayer1HGCalTkPt3, process.l1ctLayer1TkPt3)
    monitorPerf("L1PFTkPt3", "l1ctLayer1TkPt3:PF")
    monitorPerf("L1PuppiTkPt3", "l1ctLayer1TkPt3:Puppi")
    #if hasattr(process, "l1tkv5Stubs"):
    #    process.l1tkv5StubsTkPt3 = cms.EDFilter("L1TPFCandSelector", src = cms.InputTag("l1ctLayer1:TKVtx"), cut = cms.string("pfTrack.nStubs >= 5 && pt > 3"))
    #    process.extraPFStuff.add(process.l1tkv5StubsTkPt3)
    #    monitorPerf("L1TKV5TkPt3", "l1tkv5StubsTkPt3", makeRespSplit = False)

def addGenLep(pdgs=[11,13,22]):
    genLepTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
                src = cms.InputTag("genParticles"),
                doc = cms.string("gen leptons"),
                singleton = cms.bool(False), # the number of entries is variable
                extension = cms.bool(False), # this is the main table
                variables = cms.PSet(
                    pt  = Var("pt",  float,precision=8),
                    phi = Var("phi", float,precision=8),
                    eta  = Var("eta", float,precision=8),
                    vz   = Var("vz",  float,precision=8),
                    charge  = Var("charge", int, doc="charge id"),
                    prompt  = Var("2*statusFlags().isPrompt() + statusFlags().isDirectPromptTauDecayProduct()", int, doc="Particle status."),
                )
            )
    for pdgId in pdgs:
        if pdgId == 13:
            process.genMuTable = genLepTable.clone(
                        cut  = cms.string("abs(pdgId) == %d && status == 1 && pt > 2" % pdgId),
                        name = cms.string("GenMu"))
            process.extraPFStuff.add(process.genMuTable)
        elif pdgId == 11:
            process.genElTable = genLepTable.clone(
                        cut  = cms.string("abs(pdgId) == %d && status == 1 && pt > 2" % pdgId),
                        name = cms.string("GenEl"))
            process.extraPFStuff.add(process.genElTable)
        elif pdgId == 22:
            process.genPhTable = genLepTable.clone(
                        cut  = cms.string("abs(pdgId) == %d && status == 1 && pt > 8 && statusFlags().isPrompt()" % pdgId),
                        name = cms.string("GenPh"))
            process.extraPFStuff.add(process.genPhTable)

def addStaMu():
    process.staMuTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
                        src = cms.InputTag('simGmtStage2Digis'),
                        cut = cms.string(""),
                        name = cms.string("StaMu"),
                        doc = cms.string("reco leptons"),
                        singleton = cms.bool(False), # the number of entries is variable
                        extension = cms.bool(False), # this is the main table
                        variables = cms.PSet(
                            pt  = Var("pt",  float,precision=8),
                            phi = Var("phi", float,precision=8),
                            eta  = Var("eta", float,precision=8),
                            charge  = Var("charge", int, doc="charge id"),
                            quality  = Var("hwQual", int, doc="charge id"),
                        )
    )
    process.extraPFStuff.add(process.staMuTable)

def addPFLep(pdgs=[11,13,22],opts=["PF","Puppi"], postfix=""):
    for w in opts:
        pfLepTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
                        src = cms.InputTag("l1ctLayer1%s:%s"%(postfix,w)),
                        doc = cms.string("reco leptons"),
                        singleton = cms.bool(False), # the number of entries is variable
                        extension = cms.bool(False), # this is the main table
                        variables = cms.PSet(
                            pt  = Var("pt",  float,precision=8),
                            phi = Var("phi", float,precision=8),
                            eta  = Var("eta", float,precision=8),
                            vz   = Var("vz",  float,precision=8),
                            charge  = Var("charge", int, doc="charge id"),
                        )
                    )
        for pdgId in pdgs:
            if pdgId == 13:
                muTable = pfLepTable.clone(
                            cut  = cms.string("abs(pdgId) == %d && pt > 2" % pdgId),
                            name = cms.string(w+"Mu"+postfix))
                setattr(process, w+"Mu"+postfix+"Table", muTable)
                process.extraPFStuff.add(muTable)
                muTable.variables.quality = Var("? muon.isNonnull ? muon.hwQual : -1", int, doc="Quality")
            elif pdgId == 11:
                elTable = pfLepTable.clone(
                            cut  = cms.string("abs(pdgId) == %d && pt > 2" % pdgId),
                            name = cms.string(w+"El"+postfix))
                setattr(process, w+"El"+postfix+"Table", elTable)
                process.extraPFStuff.add(elTable)
                elTable.variables.hwEmID = Var("? pfCluster.isNonnull ? pfCluster.hwEmID : -1", int, doc="Quality")
                elTable.variables.pfEmID = Var("? pfCluster.isNonnull ? pfCluster.egVsPionMVAOut : -1", float, precision=8)
                elTable.variables.pfPuID = Var("? pfCluster.isNonnull ? pfCluster.egVsPUMVAOut : -1", float, precision=8)
                elTable.variables.clPt   = Var("? pfCluster.isNonnull ? pfCluster.pt : -1", float, precision=8)
                elTable.variables.clEmEt = Var("? pfCluster.isNonnull ? pfCluster.emEt : -1", float, precision=8)
            elif pdgId == 22:
                phTable = pfLepTable.clone(
                            cut  = cms.string("abs(pdgId) == %d && pt > 8" % pdgId),
                            name = cms.string(w+"Ph"+postfix))
                phTable.variables.hwEmID = Var("hwEmID", int)
                phTable.variables.pfPuID = Var("? pfCluster.isNonnull ? pfCluster.egVsPUMVAOut : -1", float, precision=8)
                elTable.variables.pfEmID = Var("? pfCluster.isNonnull ? pfCluster.egVsPionMVAOut : -1", float, precision=8)
                elTable.variables.clPt   = Var("? pfCluster.isNonnull ? pfCluster.pt : -1", float, precision=8)
                elTable.variables.clEmEt = Var("? pfCluster.isNonnull ? pfCluster.emEt : -1", float, precision=8)
                if w == "Puppi":
                    phTable.variables.puppiW = Var("puppiWeight", float, precision=8)
                setattr(process, w+"Ph"+postfix+"Table", phTable)
                process.extraPFStuff.add(phTable)
def addTkEG(postfix=""):
    for w in "EB","EE":
        tkEmTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
                        name = cms.string("TkEm"+w+postfix),
                        src = cms.InputTag("l1ctLayer1EG%s:L1TkEm%s" % (postfix, w)),
                        cut = cms.string(""),
                        doc = cms.string(""),
                        singleton = cms.bool(False), # the number of entries is variable
                        extension = cms.bool(False), # this is the main table
                        variables = cms.PSet(
                            pt  = Var("pt",  float,precision=8),
                            phi = Var("phi", float,precision=8),
                            eta  = Var("eta", float,precision=8),
                            charge  = Var("charge", int, doc="charge"),
                            emid    = Var("EGRef.hwQual", int, doc="id"),
                            tkIso   = Var("trkIsol", float, precision=8),
                            tkIsoV  = Var("trkIsolPV", float, precision=8),
                        )
                    )
        tkEleTable = tkEmTable.clone(
                        name = cms.string("TkEle"+w+postfix),
                        src = cms.InputTag("l1ctLayer1EG%s:L1TkEle%s" % (postfix, w)),
                    )
        tkEleTable.variables.charge = Var("charge", int, doc="charge")
        tkEleTable.variables.vz     = Var("trkzVtx",  float,precision=8)
        tkEleTable.variables.caloEta = Var("EGRef.eta", float,precision=8)
        tkEleTable.variables.caloPhi = Var("EGRef.phi", float,precision=8)
        setattr(process, "TkEm%s%sTable" % (w,postfix), tkEmTable)
        setattr(process, "TkEle%s%sTable" % (w,postfix), tkEleTable)
        process.extraPFStuff.add(tkEmTable,tkEleTable)

def goGun(calib=1):
    process.ntuple.isParticleGun = True
    respOnly()
    if calib: 
        addCalib()
def goMT(nthreads=2):
    process.options.numberOfThreads = cms.untracked.uint32(nthreads)
    process.options.numberOfStreams = cms.untracked.uint32(0)
def noPU():
    for X in "", "EM", "Raw", "EMRaw":
        pfc = getattr(process, 'pfClustersFromHGC3DClusters'+X, None)
        if not pfc: continue
        pfc.emVsPUID.wp = "-1.0"
def oldInputs():
    process.pfClustersFromCombinedCaloHCal.phase2barrelCaloTowers = [cms.InputTag("L1EGammaClusterEmuProducer",)]


def addEDMOutput():
    process.out = cms.OutputModule("PoolOutputModule",
                                   fileName = cms.untracked.string("debugPF.root"),
                                   SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring("p"))
                               )
    process.end = cms.EndPath(process.out)
    process.maxEvents.input = 10

def addCrops():
    process.l1ctLayer1HGCalIdealCrop = process.l1ctLayer1HGCal.clone()
    #process.l1ctLayer1HGCalIdealCrop.regionizerAlgoParameters.useAlsoVtxCoords = False
    process.l1ctLayer1HGCalIdealCrop.pfAlgoParameters.nTrack = 30
    process.l1ctLayer1HGCalIdealCrop.pfAlgoParameters.nCalo = 20
    process.l1ctLayer1HGCalIdealCrop.pfAlgoParameters.nMu = 4
    process.l1ctLayer1HGCalIdealCrop.puAlgoParameters.nFinalSort = 18
    process.l1ctLayer1HGCalRealCrop = process.l1ctLayer1HGCalIdealCrop.clone()
    process.l1ctLayer1HGCalRealCrop.regionizerAlgo = "Multififo"
    process.l1ctLayer1HGCalRealCrop.regionizerAlgoParameters = cms.PSet(
                useAlsoVtxCoords = cms.bool(True),
                nClocks = cms.uint32(54),
                nTrack = cms.uint32(30),
                nCalo = cms.uint32(20),
                nEmCalo = cms.uint32(0),
                nMu = cms.uint32(4),
            )
    process.l1ctLayer1HGCalRealCrop.pfAlgoParameters.nTrack = 30
    process.l1ctLayer1HGCalRealCrop.pfAlgoParameters.nCalo = 20
    process.l1ctLayer1HGCalRealCrop.pfAlgoParameters.nMu = 4
    process.l1ctLayer1HGCalRealCrop.pfAlgoParameters.nSelCalo = 20
    process.l1ctLayer1HGCalRealCrop.puAlgoParameters.nTrack = 30
    process.l1ctLayer1HGCalRealCrop.puAlgoParameters.nIn = 20
    process.l1ctLayer1HGCalRealCrop.puAlgoParameters.nOut = 20
    process.l1ctLayer1HGCalRealCrop.puAlgoParameters.nFinalSort = 18
    process.l1ctLayer1HGCalRealCrop.puAlgoParameters.finalSortAlgo = "Hybrid"
    for X in "IdealCrop", "RealCrop":
        setattr(process, "l1ctLayer1"+X, process.l1ctLayer1.clone(
            pfProducers = [  cms.InputTag("l1ctLayer1Barrel"), cms.InputTag("l1ctLayer1HGCal"+X), cms.InputTag("l1ctLayer1HGCalNoTK"), cms.InputTag("l1ctLayer1HF") ]))
        process.extraPFStuff.add( getattr(process, "l1ctLayer1HGCal"+X), getattr(process, "l1ctLayer1"+X) )
        monitorPerf("L1PF"+X,    "l1ctLayer1"+X+":PF")
        monitorPerf("L1Puppi"+X, "l1ctLayer1"+X+":Puppi")

def useTkInputEmulator(postfix="",bitwise=True,newTrackWord=True):
    process.pfTracksFromL1Tracks.redigitizeTrackWord = newTrackWord
    if postfix:
        prodhgcal = process.l1ctLayer1HGCal.clone()
        setattr(process, 'l1ctLayer1HGCal' + postfix, prodhgcal)
        merger = process.l1ctLayer1.clone(
            pfProducers = [  cms.InputTag("l1ctLayer1Barrel"), cms.InputTag("l1ctLayer1HGCal" + postfix), cms.InputTag("l1ctLayer1HGCalNoTK"), cms.InputTag("l1ctLayer1HF") ])
        setattr(process, 'l1ctLayer1' + postfix, merger)
        monitorPerf("L1PF"+ postfix,    "l1ctLayer1"+postfix+":PF")
        monitorPerf("L1Puppi"+ postfix, "l1ctLayer1"+postfix+":Puppi")
        egmerger = process.l1ctLayer1EG.clone()
        egmerger.tkElectrons[0].pfProducers = [ cms.InputTag("l1ctLayer1HGCal" + postfix, "L1TkEle") ]
        egmerger.tkEms[0].pfProducers = [ cms.InputTag("l1ctLayer1HGCal" + postfix, "L1TkEm"),  cms.InputTag("l1ctLayer1HGCalNoTK", "L1TkEm")  ]
        egmerger.tkEgs[0].pfProducers = [ cms.InputTag("l1ctLayer1HGCal" + postfix, "L1Eg"),  cms.InputTag("l1ctLayer1HGCalNoTK", "L1Eg")  ]
        setattr(process, 'l1ctLayer1EG' + postfix, egmerger)
        process.extraPFStuff.add(prodhgcal, merger, egmerger)
    else:
        prodhgcal = process.l1ctLayer1HGCal
    prodhgcal.trackInputConversionAlgo = "Emulator"
    prodhgcal.trackInputConversionParameters = cms.PSet(
            region = cms.string("endcap"),
            trackWordEncoding = cms.string("biased" if newTrackWord else "stepping"),
            bitwiseAccurate = cms.bool(bitwise),
            ptLUTBits = cms.uint32(11),
            etaLUTBits = cms.uint32(11),
            etaPreOffs = cms.int32(0),
            etaShift = cms.uint32(15-11),
            etaPostOffs = cms.int32(150),
            phiBits = cms.uint32(10),
            z0Bits = cms.uint32(12),
            dEtaHGCalBits = cms.uint32(10),
            dEtaHGCalZ0PreShift = cms.uint32(2),
            dEtaHGCalRInvPreShift = cms.uint32(6),
            dEtaHGCalLUTBits = cms.uint32(10),
            dEtaHGCalLUTShift = cms.uint32(2),
            dEtaHGCalFloatOffs = cms.double(0.0),
            dPhiHGCalBits = cms.uint32(4),
            dPhiHGCalZ0PreShift = cms.uint32(4),
            dPhiHGCalZ0PostShift = cms.uint32(6),
            dPhiHGCalRInvShift = cms.uint32(4),
            dPhiHGCalTanlInvShift = cms.uint32(22),
            dPhiHGCalTanlLUTBits = cms.uint32(10),
            dPhiHGCalFloatOffs = cms.double(0.0)
            )
def useMuInputEmulator(postfix=""):
    for reg in "Barrel", "HGCal", "HGCalNoTK":
        if postfix:
            prod = getattr(process, 'l1ctLayer1'+reg).clone()
            setattr(process, 'l1ctLayer1' + reg + postfix, prod)
            process.extraPFStuff.add(prod)
        else:
            prod = getattr(process, 'l1ctLayer1'+reg)
        prod.muonInputConversionAlgo = "Emulator"
        prod.muonInputConversionParameters = cms.PSet(
                z0Scale = cms.double(1.8),
                dxyScale = cms.double(3.0)
                )
    if postfix:
        merger = process.l1ctLayer1.clone(
            pfProducers = [  cms.InputTag("l1ctLayer1Barrel" + postfix), cms.InputTag("l1ctLayer1HGCal" + postfix), cms.InputTag("l1ctLayer1HGCalNoTK" + postfix), cms.InputTag("l1ctLayer1HF") ])
        setattr(process, 'l1ctLayer1' + postfix, merger)
        monitorPerf("L1PF"+ postfix,    "l1ctLayer1"+postfix+":PF")
        monitorPerf("L1Puppi"+ postfix, "l1ctLayer1"+postfix+":Puppi")
        process.extraPFStuff.add(merger)

if False:
    #process.source.fileNames  = [ '/store/cmst3/group/l1tr/gpetrucc/11_1_0/NewInputs110X/110121.done/TTbar_PU200/inputs110X_%d.root' % i for i in (1,)] #3,7,8,9) ]
    #process.source.fileNames  = [ '/store/cmst3/group/l1tr/gpetrucc/11_1_0/NewInputs110X/110121.done/DYToLL_PU200/inputs110X_%d.root' % i for i in (1,2,3)] #3,7,8,9) ]
    process.source.fileNames  = [ '/store/cmst3/group/l1tr/gpetrucc/11_1_0/NewInputs110X/110121.done/DYToLL_PU200/inputs110X_%d.root' % i for i in (1,)] #3,7,8,9) ]
    addOld()
    #goMT(4)
    addStaMu()
    addGenLep([13])
    addPFLep([13],["PF","OldPF"])

    if True:
        #process.source.eventsToProcess = cms.untracked.VEventRange("1:1376:240626","1:1376:240627","1:1376:240628","1:1376:240642","1:1376:240638")
        process.source.eventsToProcess = cms.untracked.VEventRange("1:1376:240656", "1:1376:240742", "1:1108:193756", "1:1108:193762", "1:1108:193772")
        #process.source.eventsToProcess = cms.untracked.VEventRange()
        #process.maxEvents.input = 10
        addEDMOutput()
        R='HGCal'
        getattr(process, 'l1ctLayer1'  +R).trkPtCut = 10
        getattr(process, 'l1pfProducer'+R).trkPtCut = 10
        getattr(process, 'l1pfProducer'+R).debug = 2
        getattr(process, 'l1ctLayer1'  +R).pfAlgoParameters.debug = True

def saveCands():
    process.l1pfcandTable = cms.EDProducer("L1PFCandTableProducer",
                                           commonSel = cms.string("pt > 0.0 && abs(eta) < 10.0"),
                                           cands = cms.PSet(
                                           ),
                                           moreVariables = cms.PSet(
                                               puppiWeight = cms.string("puppiWeight"),
                                               pdgId = cms.string("pdgId"),
                                               charge = cms.string("charge")
                                           ),
                                       )
    monitorPerf("L1PF", "l1ctLayer1:PF", saveCands=True)
    monitorPerf("L1Puppi", "l1ctLayer1:Puppi", saveCands=True)
    process.p += process.l1pfcandTable

def saveGenCands():
    process.gencandTable = cms.EDProducer("L1PFCandTableProducer",
                                           commonSel = cms.string("pt > 0.0"),
                                           cands = cms.PSet(
                                               Gen = cms.InputTag("genParticlesForMETAllVisible")
                                           ),
                                           moreVariables = cms.PSet(
                                               pdgId = cms.string("pdgId"),
                                               charge = cms.string("charge")
                                           ),
                                      )
    process.p += process.gencandTable

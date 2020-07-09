import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras
from math import sqrt

process = cms.Process("RESP", eras.Phase2C8_trigger)

process.load('Configuration.StandardSequences.Services_cff')
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True), allowUnscheduled = cms.untracked.bool(False) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1))
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:inputs106X.root'),
    inputCommands = cms.untracked.vstring("keep *", 
            "drop l1tPFClusters_*_*_*",
            "drop l1tPFTracks_*_*_*",
            "drop l1tPFCandidates_*_*_*")
)

process.load('Configuration.Geometry.GeometryExtended2026D41Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2026D41_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff') # needed to read HCal TPs
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('RecoMET.Configuration.GenMETParticles_cff')
process.load('RecoMET.METProducers.genMetTrue_cfi')

from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets
from RecoMET.METProducers.PFMET_cfi import pfMet

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

process.load("L1Trigger.Phase2L1ParticleFlow.l1ParticleFlow_cff")

process.extraPFStuff = cms.Task()

process.runPF = cms.Sequence( 
    process.l1ParticleFlow 
)

process.centralGen = cms.EDFilter("CandPtrSelector", src = cms.InputTag("genParticlesForMETAllVisible"), cut = cms.string("abs(eta) < 2.4"))
process.barrelGen = cms.EDFilter("CandPtrSelector", src = cms.InputTag("genParticlesForMETAllVisible"), cut = cms.string("abs(eta) < 1.5"))
process.genMetCentralTrue = process.genMetTrue.clone(src = cms.InputTag("centralGen"))
process.genMetBarrelTrue = process.genMetTrue.clone(src = cms.InputTag("barrelGen"))
process.extraPFStuff.add(
    process.genParticlesForMETAllVisible,
    process.centralGen,
    process.barrelGen,
    process.genMetCentralTrue,
    process.genMetBarrelTrue
)

def monitorPerf(label, tag, makeResp=True, makeRespSplit=True, makeJets=True, makeMET=True, makeCentralMET=True, makeBarrelMET=True,
                makeInputMultiplicities=False, makeOutputMultiplicities=False):
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
    if makeMET:
        _add('met'+label, pfMet.clone(src = tag, calculateSignificance = False))
        setattr(process.l1pfmetTable.mets, label, cms.InputTag('met'+label))
        if makeCentralMET:
            _add('central'+label, cms.EDFilter("CandPtrSelector", src = cms.InputTag(tag), cut = cms.string("abs(eta) < 2.4")))
            _add('met'+label+'Central', pfMet.clone(src = 'central'+label, calculateSignificance = False))
            setattr(process.l1pfmetCentralTable.mets, label, cms.InputTag('met'+label+'Central'))
        if makeBarrelMET:
            _add('barrel'+label, cms.EDFilter("CandPtrSelector", src = cms.InputTag(tag), cut = cms.string("abs(eta) < 1.5")))
            _add('met'+label+'Barrel', pfMet.clone(src = 'barrel'+label, calculateSignificance = False))
            setattr(process.l1pfmetBarrelTable.mets, label, cms.InputTag('met'+label+'Barrel'))
    if makeInputMultiplicities:
        D = tag.split(":")[0] # l1pfProducer[Barrel,HGCal,HF] usually
        I = tag.split(":")[1] # Calo, EmCalo, TK, or Mu, usually
        for X in ["tot","max"]:
            process.ntuple.copyUInts.append( "%s:%sNL1%s" % (D,X,I))
        process.ntuple.copyVecUInts.append( "%s:vecNL1%s" % (D,I))
    if makeOutputMultiplicities:
        D = tag.split(":")[0] # l1pfProducer[Barrel,HGCal,HF] usually
        P = tag.split(":")[1] # PF or Puppi, usually
        for O in [""] + "Charged Neutral ChargedHadron NeutralHadron Photon Electron Muon".split():
            for X in ["tot","max"]:
                process.ntuple.copyUInts.append( "%s:%sNL1%s%s" % (D,X,P,O))
            process.ntuple.copyVecUInts.append( "%s:vecNL1%s%s" % (D,P,O))

process.ntuple = cms.EDAnalyzer("ResponseNTuplizer",
    genJets = cms.InputTag("ak4GenJetsNoNu"),
    genParticles = cms.InputTag("genParticles"),
    isParticleGun = cms.bool(False),
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
process.l1pfmetBarrelTable  = process.l1pfmetTable.clone(genMet = "genMetBarrelTrue", flavour = "Barrel")

monitorPerf("L1Calo", "l1pfCandidates:Calo", makeRespSplit = False)
monitorPerf("L1TK", "l1pfCandidates:TK", makeRespSplit = False, makeJets=False, makeMET=False)
monitorPerf("L1TKV", "l1pfCandidates:TKVtx", makeRespSplit = False, makeJets=False, makeMET=False)
monitorPerf("L1PF", "l1pfCandidates:PF")
monitorPerf("L1Puppi", "l1pfCandidates:Puppi")

for D in ['Barrel','HF','HGCal','HGCalNoTK']:
    monitorPerf("L1%sCalo"%D,"l1pfProducer%s:Calo"%D, makeResp=False, makeRespSplit=False, makeJets=False, makeMET=False, 
               makeCentralMET=False, makeBarrelMET=False, makeInputMultiplicities=True)
    monitorPerf("L1%sEmCalo"%D,"l1pfProducer%s:EmCalo"%D, makeResp=False, makeRespSplit=False, makeJets=False, makeMET=False, 
               makeCentralMET=False, makeBarrelMET=False, makeInputMultiplicities=True)
    monitorPerf("L1%sTK"%D,"l1pfProducer%s:TK"%D, makeResp=False, makeRespSplit=False, makeJets=False, makeMET=False, 
               makeCentralMET=False, makeBarrelMET=False, makeInputMultiplicities=True)
    monitorPerf("L1%sMu"%D,"l1pfProducer%s:Mu"%D, makeResp=False, makeRespSplit=False, makeJets=False, makeMET=False, 
               makeCentralMET=False, makeBarrelMET=False, makeInputMultiplicities=True)

    monitorPerf("L1%sPF"%D,"l1pfProducer%s:PF"%D, makeResp=False, makeRespSplit=False, makeJets=False, makeMET=False, 
               makeCentralMET=False, makeBarrelMET=False, makeOutputMultiplicities=True)
    monitorPerf("L1%sPuppi"%D,"l1pfProducer%s:Puppi"%D, makeResp=False, makeRespSplit=False, makeJets=False, makeMET=False, 
               makeCentralMET=False, makeBarrelMET=False, makeOutputMultiplicities=True)

# define regions
def goRegional(postfix="", relativeCoordinates=False):
    overlap=0.25 # 0.3
    getattr(process, 'l1pfProducer'+postfix+'Barrel').regions = cms.VPSet(
        cms.PSet(
            etaBoundaries = cms.vdouble(-1.5, -0.5, 0.5, 1.5),
            etaExtra = cms.double(overlap),
            phiExtra = cms.double(overlap),
            phiSlices = cms.uint32(9)
        )
    )
    getattr(process, 'l1pfProducer'+postfix+'HGCalNoTK').regions = cms.VPSet(
        cms.PSet(
            etaBoundaries = cms.vdouble(-3, -2.5),
            etaExtra = cms.double(overlap),
            phiExtra = cms.double(overlap),
            phiSlices = cms.uint32(9)
        ),
        cms.PSet(
            etaBoundaries = cms.vdouble(2.5, 3),
            etaExtra = cms.double(overlap),
            phiExtra = cms.double(overlap),
            phiSlices = cms.uint32(9)
        )
    )
    getattr(process, 'l1pfProducer'+postfix+'HGCal').regions = cms.VPSet(
        cms.PSet(
            etaBoundaries = cms.vdouble(-2.5, -1.5),
            etaExtra = cms.double(overlap),
            phiExtra = cms.double(overlap),
            phiSlices = cms.uint32(9)
        ),
        cms.PSet(
            etaBoundaries = cms.vdouble(1.5, 2.5),
            etaExtra = cms.double(overlap),
            phiExtra = cms.double(overlap),
            phiSlices = cms.uint32(9)
        )
    )
    getattr(process, 'l1pfProducer'+postfix+'HF').regions = cms.VPSet(
        cms.PSet(
            etaBoundaries = cms.vdouble(-5, -4, -3),
            etaExtra = cms.double(overlap),
            phiExtra = cms.double(overlap),
            phiSlices = cms.uint32(9)
        ),
        cms.PSet(
            etaBoundaries = cms.vdouble(3, 4, 5),
            etaExtra = cms.double(overlap),
            phiExtra = cms.double(overlap),
            phiSlices = cms.uint32(9)
        )
    )
    for D in 'Barrel', 'HGCal', 'HGCalNoTK', 'HF':
        getattr(process, 'l1pfProducer'+postfix+D).useRelativeRegionalCoordinates = relativeCoordinates

process.runPF.associate(process.extraPFStuff)
# to check available tags:
#process.content = cms.EDAnalyzer("EventContentAnalyzer")
process.p = cms.Path(
        process.runPF + 
        process.ntuple + #process.content +
        process.l1pfjetTable + 
        process.l1pfmetTable + process.l1pfmetCentralTable + process.l1pfmetBarrelTable
        )
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
    process.ntuple.objects.MuGenAcc = cms.VInputTag(cms.InputTag("genInAcceptance"))
    process.ntuple.objects.MuGenAcc_sel = cms.string("abs(pdgId) == 13")
    process.extraPFStuff.add(process.genInAcceptance)
if False: # test also PF leptons
    process.ntuple.objects.L1PFMuon = cms.VInputTag("l1pfCandidates:PF",)
    process.ntuple.objects.L1PFMuon_sel = cms.string("abs(pdgId) == 13")
    process.ntuple.objects.L1PFElectron = cms.VInputTag("l1pfCandidates:PF",)
    process.ntuple.objects.L1PFElectron_sel = cms.string("abs(pdgId) == 11")
def respOnly():
    process.p.remove(process.l1pfjetTable)
    process.p.remove(process.l1pfmetTable)
    process.p.remove(process.l1pfmetCentralTable)
    process.p.remove(process.l1pfmetBarrelTable)
    process.end.remove(process.outnano)
def addCHS():
    process.l1PuppiCharged = cms.EDFilter("L1TPFCandSelector",
        src = cms.InputTag("l1pfCandidates:Puppi"),
        cut = cms.string("charge != 0"))
    process.l1PFNeutral = cms.EDFilter("L1TPFCandSelector",
        src = cms.InputTag("l1pfCandidates:PF"),
        cut = cms.string("charge == 0"))
    process.extraPFStuff.add(process.l1PuppiCharged, process.l1PFNeutral)
    monitorPerf("L1CHS", [ "l1PuppiCharged", "l1PFNeutral" ], makeRespSplit = False)
def addPFnoMu():
    process.l1pfProducerBarrelPFnoMu = process.l1pfProducerBarrel.clone()
    process.l1pfProducerBarrelPFnoMu.useStandaloneMuons       = cms.bool(False) 
    process.l1pfProducerBarrelPFnoMu.useTrackerMuons          = cms.bool(False)
    #
    process.l1pfProducerHGCalPFnoMu = process.l1pfProducerHGCal.clone()
    process.l1pfProducerHGCalPFnoMu.useStandaloneMuons       = cms.bool(False) 
    process.l1pfProducerHGCalPFnoMu.useTrackerMuons          = cms.bool(False)
    #
    process.l1pfProducerHFPFnoMu = process.l1pfProducerHF.clone()
    process.l1pfProducerHFPFnoMu.useStandaloneMuons       = cms.bool(False) 
    process.l1pfProducerHFPFnoMu.useTrackerMuons          = cms.bool(False)
    #
    process.l1pfCandidatesPFnoMu = process.l1pfCandidates.clone(
            pfProducers = ["l1pfProducerBarrelPFnoMu", "l1pfProducerHGCalPFnoMu", "l1pfProducerHFPFnoMu"], 
            labelsToMerge = [ "Puppi", "PF" ])
    process.extraPFStuff.add(process.l1pfProducerBarrelPFnoMu, process.l1pfProducerHGCalPFnoMu, process.l1pfProducerHFPFnoMu, process.l1pfCandidatesPFnoMu)
    monitorPerf("L1PuppinoMu", "l1pfCandidatesPFnoMu:Puppi")
    monitorPerf("L1PFnoMu", "l1pfCandidatesPFnoMu:PF")
def addPFtkMu():
    process.l1pfProducerBarrelPFtkMu = process.l1pfProducerBarrel.clone()
    process.l1pfProducerBarrelPFtkMu.useStandaloneMuons       = cms.bool(False) 
    process.l1pfProducerBarrelPFtkMu.useTrackerMuons          = cms.bool(True)
    #
    process.l1pfProducerHGCalPFtkMu = process.l1pfProducerHGCal.clone()
    process.l1pfProducerHGCalPFtkMu.useStandaloneMuons       = cms.bool(False) 
    process.l1pfProducerHGCalPFtkMu.useTrackerMuons          = cms.bool(True)
    #
    process.l1pfProducerHFPFtkMu = process.l1pfProducerHF.clone()
    process.l1pfProducerHFPFtkMu.useStandaloneMuons       = cms.bool(False) 
    process.l1pfProducerHFPFtkMu.useTrackerMuons          = cms.bool(True)
    #
    process.l1pfCandidatesPFtkMu = process.l1pfCandidates.clone(
            pfProducers = ["l1pfProducerBarrelPFtkMu", "l1pfProducerHGCalPFtkMu", "l1pfProducerHFPFtkMu"], 
            labelsToMerge = [ "Puppi", "PF" ])
    process.extraPFStuff.add(process.l1pfProducerBarrelPFtkMu, process.l1pfProducerHGCalPFtkMu, process.l1pfProducerHFPFtkMu, process.l1pfCandidatesPFtkMu)
    monitorPerf("L1PuppitkMu", "l1pfCandidatesPFtkMu:Puppi")
    monitorPerf("L1PFtkMu", "l1pfCandidatesPFtkMu:PF")
def addPuppiOld():
    process.l1pfProducerBarrelPuppiOld = process.l1pfProducerBarrel.clone(puAlgo = "Puppi")
    process.l1pfProducerBarrelPuppiOld.puppiEtaCuts       = cms.vdouble(1.5) 
    process.l1pfProducerBarrelPuppiOld.puppiPtCuts        = cms.vdouble(0.0)
    process.l1pfProducerBarrelPuppiOld.puppiPtCutsPhotons = cms.vdouble(0.0)
    #
    process.l1pfProducerHGCalPuppiOld = process.l1pfProducerHGCal.clone(puAlgo = "Puppi")
    process.l1pfProducerHGCalPuppiOld.puppiEtaCuts       = cms.vdouble(2.5, 2.85,  3.0) 
    process.l1pfProducerHGCalPuppiOld.puppiPtCuts        = cms.vdouble( 20,  40,  9999)
    process.l1pfProducerHGCalPuppiOld.puppiPtCutsPhotons = cms.vdouble( 20,  40,  9999)
    #
    process.l1pfProducerHFPuppiOld = process.l1pfProducerHF.clone(puAlgo = "Puppi")
    process.l1pfProducerHFPuppiOld.puppiEtaCuts       = cms.vdouble(5.5) 
    process.l1pfProducerHFPuppiOld.puppiPtCuts        = cms.vdouble( 30)
    process.l1pfProducerHFPuppiOld.puppiPtCutsPhotons = cms.vdouble( 30)
    #
    process.l1pfCandidatesPuppiOld = process.l1pfCandidates.clone(
            pfProducers = ["l1pfProducerBarrelPuppiOld", "l1pfProducerHGCalPuppiOld", "l1pfProducerHFPuppiOld"], 
            labelsToMerge = [ "Puppi" ])
    process.extraPFStuff.add(process.l1pfProducerBarrelPuppiOld, process.l1pfProducerHGCalPuppiOld, process.l1pfProducerHFPuppiOld, process.l1pfCandidatesPuppiOld)
    monitorPerf("L1PuppiOld", "l1pfCandidatesPuppiOld:Puppi")
def addTKs():
    process.l1tkv5Stubs = cms.EDFilter("L1TPFCandSelector", src = cms.InputTag("l1pfCandidates:TKVtx"), cut = cms.string("pfTrack.nStubs >= 5"))
    process.l1tkv6Stubs = cms.EDFilter("L1TPFCandSelector", src = cms.InputTag("l1pfCandidates:TKVtx"), cut = cms.string("pfTrack.nStubs >= 6"))
    process.extraPFStuff.add(process.l1tkv5Stubs, process.l1tkv6Stubs)
    monitorPerf("L1TKV5", "l1tkv5Stubs", makeRespSplit = False)
    monitorPerf("L1TKV6", "l1tkv6Stubs", makeRespSplit = False)
    monitorPerf("L1TK", "l1pfCandidates:TK", makeRespSplit = False)
    monitorPerf("L1TKV", "l1pfCandidates:TKVtx", makeRespSplit = False)
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

def addRegional(label="Regional", relativeCoordinates=False):
    # clone each part
    for D in ['Barrel','HF','HGCal','HGCalNoTK']:
        pfreg = getattr(process, 'l1pfProducer'+D).clone() 
        setattr(process, 'l1pfProducer'+label+D, pfreg)
        process.extraPFStuff.add(pfreg)
    # merging
    pfcomb = cms.EDProducer("L1TPFCandMultiMerger",
            pfProducers = cms.VInputTag(
                cms.InputTag("l1pfProducer"+label+"Barrel"), 
                cms.InputTag("l1pfProducer"+label+"HGCal"),
                cms.InputTag("l1pfProducer"+label+"HGCalNoTK"),
                cms.InputTag("l1pfProducer"+label+"HF")
                ),
            labelsToMerge = cms.vstring("Calo", "TK", "TKVtx", "PF", "Puppi"),
    )
    setattr(process, 'l1pfCandidates'+label, pfcomb)
    process.extraPFStuff.add(pfcomb)
    # monitoring
    monitorPerf("L1Calo"+label, "l1pfCandidates"+label+":Calo", makeRespSplit = False)
    monitorPerf("L1TK"+label, "l1pfCandidates"+label+":TK", makeRespSplit = False, makeJets=False, makeMET=False)
    monitorPerf("L1TKV"+label, "l1pfCandidates"+label+":TKVtx", makeRespSplit = False, makeJets=False, makeMET=False)
    monitorPerf("L1PF"+label, "l1pfCandidates"+label+":PF")
    monitorPerf("L1Puppi"+label, "l1pfCandidates"+label+":Puppi")
    # and finally go regional
    goRegional(label, relativeCoordinates=relativeCoordinates)
def firmwareLike():
    ## Make a copy called "CMSSW", and then configure the setup to be close ro firmware
    for det in "Barrel", "HGCal", "HGCalNoTK", "HF":
        l1pf = getattr(process, 'l1pfProducer'+det).clone()
        setattr(process, 'l1pfProducerCMSSW'+det, l1pf)
        process.extraPFStuff.add(l1pf)
    process.l1pfCandidatesCMSSW = cms.EDProducer("L1TPFCandMultiMerger",
            pfProducers = cms.VInputTag(
                cms.InputTag("l1pfProducerCMSSWBarrel"), 
                cms.InputTag("l1pfProducerCMSSWHGCal"),
                cms.InputTag("l1pfProducerCMSSWHGCalNoTK"),
                cms.InputTag("l1pfProducerCMSSWHF")
                ),
            labelsToMerge = cms.vstring("PF","Puppi"),
    )
    process.extraPFStuff.add(process.l1pfCandidatesCMSSW)
    monitorPerf("L1CMSSWPF", "l1pfCandidatesCMSSW:PF")        
    monitorPerf("L1CMSSWPuppi", "l1pfCandidatesCMSSW:Puppi")        
    process.l1pfProducerBarrel.linking.ecalPriority = False # FIXME (was never implemented)
    process.l1pfProducerBarrel.linking.trackEmMayUseCaloMomenta = False # FIXME (was never implemented)
    process.l1pfProducerBarrel.linking.emCaloUseAlsoCaloSigma = False # FIXME (was never implemented, but could be)
    process.l1pfProducerBarrel.linking.emCaloSubtractionPtSlope = 1.0 # FIXME (could be implemented)
    for det in "Barrel", "HGCal", "HGCalNoTK", "HF":
        l1pf = getattr(process, 'l1pfProducer'+det)
        l1pf.linking.trackMuMatch = "drBestByPtDiff" # FIXME (could be implemented, but expensive)
        l1pf.linking.trackCaloLinkMetric = "bestByDR2Pt2" # FIXME (could be implemented, but very expensive - needs sqrt)
def addBitwise(label="Bitwise",alsoPuppi="linpuppi"):
    if not hasattr(process, "l1pfCandidatesRegional"):
        addRegional(relativeCoordinates=True)
    barrel = process.l1pfProducerRegionalBarrel.clone(pfAlgo = "BitwisePFAlgo",
            bitwiseAlgo = cms.string("pfalgo3"),
            bitwiseConfig = cms.PSet(
                NTRACK = cms.uint32(99),
                NEMCALO = cms.uint32(99),
                NCALO = cms.uint32(99),
                NMU = cms.uint32(9),
                NPHOTON = cms.uint32(99),
                NSELCALO = cms.uint32(99),
                NALLNEUTRAL = cms.uint32(199),
                DR2MAX_TK_MU = cms.uint32(2101),
                DR2MAX_TK_EM = cms.uint32(84),
                DR2MAX_EM_CALO = cms.uint32(525),
                DR2MAX_TK_CALO = cms.uint32(1182),
                TK_MAXINVPT_LOOSE = cms.uint32(40),
                TK_MAXINVPT_TIGHT = cms.uint32(80),
            ))
    hgcal = process.l1pfProducerRegionalHGCal.clone(pfAlgo = "BitwisePFAlgo",
            bitwiseAlgo = cms.string("pfalgo2hgc"),
            bitwiseConfig = cms.PSet(
                NTRACK = cms.uint32(99),
                NCALO = cms.uint32(99),
                NMU = cms.uint32(9),
                NSELCALO = cms.uint32(99),
                DR2MAX_TK_MU = cms.uint32(2101),
                DR2MAX_TK_CALO = cms.uint32(525),
                TK_MAXINVPT_LOOSE = cms.uint32(40),
                TK_MAXINVPT_TIGHT = cms.uint32(80),
            ))
    setattr(process, "l1pfProducer%sBarrel" % label, barrel)
    setattr(process, "l1pfProducer%sHGCal" % label, hgcal)
    process.extraPFStuff.add(barrel, hgcal)
    if alsoPuppi:
        barrel.puAlgo = "BitwisePuppiAlgo"
        #barrel.debugPuppi = cms.untracked.int32(1)
        barrel.bitwisePUAlgo = cms.string("linpuppi_flt" if alsoPuppi == "flt" else "linpuppi")
        barrel.bitwisePUConfig = cms.PSet(
                nTrack = cms.uint32(99),
                nIn = cms.uint32(99),
                nOut = cms.uint32(99),
                dR2Max = cms.uint32(4727),
                dR2Min = cms.uint32(257),
                ptMax  = cms.uint32(200),
                dzCut  = cms.uint32(10),
                absEtaBins = cms.vint32(),
                ptSlopeNe  = cms.vdouble(0.3),
                ptSlopePh  = cms.vdouble(0.3),
                ptZeroNe   = cms.vdouble(4.0),
                ptZeroPh   = cms.vdouble(2.5),
                alphaSlope = cms.vdouble(0.7),
                alphaZero  = cms.vdouble(6.0),
                alphaCrop  = cms.vdouble(4.0),
                priorNe    = cms.vdouble(5.0),
                priorPh    = cms.vdouble(1.0),
                ptCut      = cms.vuint32(4),
                )
        hgcal.puAlgo = "BitwisePuppiAlgo"
        #hgcal.debugPuppi = cms.untracked.int32(1)
        hgcal.bitwisePUAlgo = cms.string("linpuppi_flt" if alsoPuppi == "flt" else "linpuppi")
        hgcal.bitwisePUConfig = cms.PSet(
                nTrack = cms.uint32(99),
                nIn = cms.uint32(99),
                nOut = cms.uint32(99),
                dR2Max = cms.uint32(4727),
                dR2Min = cms.uint32(84),
                ptMax  = cms.uint32(200),
                dzCut  = cms.uint32(40),
                absEtaBins = cms.vint32(0),
                ptSlopeNe  = cms.vdouble(0.3,0.3),
                ptSlopePh  = cms.vdouble(0.4,0.4),
                ptZeroNe   = cms.vdouble(5.0,7.0),
                ptZeroPh   = cms.vdouble(3.0,4.0),
                alphaSlope = cms.vdouble(1.5,1.5),
                alphaZero  = cms.vdouble(6.0,6.0),
                alphaCrop  = cms.vdouble(3.0,3.0),
                priorNe    = cms.vdouble(5.0,5.0),
                priorPh    = cms.vdouble(1.5,1.5),
                ptCut      = cms.vuint32(4, 8),
                )
        hgcalNoTK = process.l1pfProducerRegionalHGCalNoTK.clone(puAlgo = "BitwisePuppiAlgo")
        #hgcalNoTK.debugPuppi = cms.untracked.int32(1)
        hgcalNoTK.bitwisePUAlgo = cms.string("fwdlinpuppi_flt" if alsoPuppi == "flt" else "fwdlinpuppi")
        hgcalNoTK.bitwisePUConfig = cms.PSet(
                nTrack = cms.uint32(99),
                nIn = cms.uint32(99),
                nOut = cms.uint32(99),
                dR2Max = cms.uint32(4727),
                dR2Min = cms.uint32(84),
                ptMax  = cms.uint32(200),
                dzCut  = cms.uint32(40),
                absEtaBins = cms.vint32(),
                ptSlopeNe  = cms.vdouble(0.3),
                ptSlopePh  = cms.vdouble(0.4),
                ptZeroNe   = cms.vdouble(9.0),
                ptZeroPh   = cms.vdouble(5.0),
                alphaSlope = cms.vdouble(2.2),
                alphaZero  = cms.vdouble(9.0),
                alphaCrop  = cms.vdouble(4.0),
                priorNe    = cms.vdouble(7.0),
                priorPh    = cms.vdouble(5.0),
                ptCut      = cms.vuint32(16),
                )
        hf = process.l1pfProducerRegionalHF.clone(puAlgo = "BitwisePuppiAlgo")
        #hf.debugPuppi = cms.untracked.int32(1)
        hf.bitwisePUAlgo = cms.string("fwdlinpuppi_flt" if alsoPuppi == "flt" else "fwdlinpuppi")
        hf.bitwisePUConfig = cms.PSet(
                nTrack = cms.uint32(99),
                nIn = cms.uint32(99),
                nOut = cms.uint32(99),
                dR2Max = cms.uint32(4727),
                dR2Min = cms.uint32(525),
                ptMax  = cms.uint32(400),
                dzCut  = cms.uint32(40),
                absEtaBins = cms.vint32(),
                ptSlopeNe  = cms.vdouble(0.25),
                ptSlopePh  = cms.vdouble(0.25),
                ptZeroNe   = cms.vdouble(14. ),
                ptZeroPh   = cms.vdouble(14. ),
                alphaSlope = cms.vdouble(0.6 ),
                alphaZero  = cms.vdouble(9.0 ),
                alphaCrop  = cms.vdouble(4.0 ),
                priorNe    = cms.vdouble(6.0 ),
                priorPh    = cms.vdouble(6.0 ),
                ptCut      = cms.vuint32(40),
                )
        setattr(process, "l1pfProducer%sHGCalNoTK" % label, hgcalNoTK)
        setattr(process, "l1pfProducer%sHF" % label, hf)
        process.extraPFStuff.add(hgcalNoTK, hf)
        pfcands = cms.EDProducer("L1TPFCandMultiMerger",
                pfProducers = cms.VInputTag(
                    cms.InputTag("l1pfProducer%sBarrel" % label), 
                    cms.InputTag("l1pfProducer%sHGCal" % label),
                    cms.InputTag("l1pfProducer%sHGCalNoTK" % label),
                    cms.InputTag("l1pfProducer%sHF" % label)
                    ),
                labelsToMerge = cms.vstring("PF", "Puppi"),
        )
    else:
        pfcands = cms.EDProducer("L1TPFCandMultiMerger",
                pfProducers = cms.VInputTag(
                    cms.InputTag("l1pfProducer%sBarrel" % label), 
                    cms.InputTag("l1pfProducer%sHGCal" % label),
                    cms.InputTag("l1pfProducerRegionalHGCalNoTK"),
                    cms.InputTag("l1pfProducerRegionalHF")
                    ),
                labelsToMerge = cms.vstring("PF", "Puppi"),
        ) 
    setattr(process, "l1pfCandidates%s" % label, pfcands)
    process.extraPFStuff.add(pfcands)
    monitorPerf("L1PF%s" % label, "l1pfCandidates%s:PF" % label)
    monitorPerf("L1Puppi%s" % label, "l1pfCandidates%s:Puppi"% label)


def addRefs(calo=True,tk=True):
    process.load('L1Trigger.L1CaloTrigger.Phase1L1TJets_cff')
    process.Phase1L1TJetProducer.inputCollectionTag = cms.InputTag("l1pfCandidates", "Puppi") # make sure the process name is not pre-encoded
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
    process.l1pfProducerBarrelTkPt3 = process.l1pfProducerBarrel.clone(trkPtCut = ptCut)
    process.l1pfProducerHGCalTkPt3 = process.l1pfProducerHGCal.clone(trkPtCut = ptCut)
    process.l1pfCandidatesTkPt3 = cms.EDProducer("L1TPFCandMultiMerger",
        pfProducers = cms.VInputTag(
            cms.InputTag("l1pfProducerBarrelTkPt3"), 
            cms.InputTag("l1pfProducerHGCalTkPt3"),
            cms.InputTag("l1pfProducerHGCalNoTK"),
            cms.InputTag("l1pfProducerHF")
            ),
        labelsToMerge = cms.vstring("PF", "Puppi"),
    )
    process.extraPFStuff.add(process.l1pfProducerBarrelTkPt3, process.l1pfProducerHGCalTkPt3, process.l1pfCandidatesTkPt3)
    monitorPerf("L1PFTkPt3", "l1pfCandidatesTkPt3:PF")
    monitorPerf("L1PuppiTkPt3", "l1pfCandidatesTkPt3:Puppi")
    if hasattr(process, "l1tkv5Stubs"):
        process.l1tkv5StubsTkPt3 = cms.EDFilter("L1TPFCandSelector", src = cms.InputTag("l1pfCandidates:TKVtx"), cut = cms.string("pfTrack.nStubs >= 5 && pt > 3"))
        process.extraPFStuff.add(process.l1tkv5StubsTkPt3)
        monitorPerf("L1TKV5TkPt3", "l1tkv5StubsTkPt3", makeRespSplit = False)
 
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

def doDumpFile(basename="TTbar_PU200"):
    goRegional(relativeCoordinates=True)
    for det in "Barrel", "HGCal", "HGCalNoTK", "HF":
        l1pf = getattr(process, 'l1pfProducer'+det)
        l1pf.dumpFileName = cms.untracked.string(basename+"_"+det+".dump")
        l1pf.genOrigin = cms.InputTag("genParticles","xyz0")
    process.maxEvents.input = 100



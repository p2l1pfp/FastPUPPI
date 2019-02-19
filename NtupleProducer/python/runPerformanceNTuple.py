import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras
from math import exp

process = cms.Process("RESP", eras.Phase2C4_trigger)

process.load('Configuration.StandardSequences.Services_cff')
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True), allowUnscheduled = cms.untracked.bool(False) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100))
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:/eos/cms/store/cmst3/user/gpetrucc/l1tr/105X/NewInputs104X/150219/TTbar_PU0/inputs104X_TTbar_PU0_job15.root'),
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
    skipBadFiles = cms.untracked.bool(True)
)

process.load('Configuration.Geometry.GeometryExtended2023D35Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2023D35_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff') # needed to read HCal TPs
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('RecoMET.Configuration.GenMETParticles_cff')
process.load('RecoMET.METProducers.genMetTrue_cfi')

from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets
from RecoMET.METProducers.PFMET_cfi import pfMet
pfMet.calculateSignificance = False

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '103X_upgrade2023_realistic_v2', '') 

process.load("L1Trigger.Phase2L1ParticleFlow.l1ParticleFlow_split_cff")

process.extraPFStuff = cms.Task()

process.runPF = cms.Sequence( 
    process.l1ParticleFlow_proper # excludes the prerequisites (3D clusters and L1EG clusters)
)

process.centralGen = cms.EDFilter("CandPtrSelector", src = cms.InputTag("genParticlesForMETAllVisible"), cut = cms.string("abs(eta) < 2.4"))
process.genMetCentralTrue = process.genMetTrue.clone(src = cms.InputTag("centralGen"))
process.extraPFStuff.add(
    process.genParticlesForMETAllVisible,
    process.centralGen,
    process.genMetCentralTrue
)

def monitorPerf(label, tag, makeResp=True, makeRespSplit=True, makeJets=True, makeMET=True, makeCentralMET=True):
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
    if makeJets:
        _add('ak4'+label, ak4PFJets.clone(src = tag))
        setattr(process.l1pfjetTable.jets, label, cms.InputTag('ak4'+label))
    if makeMET:
        _add('met'+label, pfMet.clone(src = tag))
        setattr(process.l1pfmetTable.mets, label, cms.InputTag('met'+label))
        if makeCentralMET:
            _add('central'+label, cms.EDFilter("CandPtrSelector", src = cms.InputTag(tag), cut = cms.string("abs(eta) < 2.4")))
            _add('met'+label+'Central', pfMet.clone(src = 'central'+label))
            setattr(process.l1pfmetTable.mets, label+'Central', cms.InputTag('met'+label+'Central'))

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
)

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
    genMetCentral = cms.InputTag("genMetCentralTrue"),
    mets = cms.PSet(
    ),
)


monitorPerf("L1Calo", "l1pfCandidates:Calo", makeRespSplit = False)
monitorPerf("L1TK", "l1pfCandidates:TK", makeRespSplit = False)
monitorPerf("L1TKV", "l1pfCandidates:TKVtx", makeRespSplit = False)
monitorPerf("L1PF", "l1pfCandidates:PF")
monitorPerf("L1Puppi", "l1pfCandidates:Puppi")
if True:
    process.l1PuppiForMET = cms.EDFilter("L1TPFCandSelector", 
            src = cms.InputTag("l1pfCandidates:Puppi"), 
            cut = cms.string("charge != 0 || abs(eta) < 1.5 || (pt > 20 && abs(eta) < 2.5) || (pt > 40 && 2.5 <= abs(eta) <= 2.85) || (pt > 30 && abs(eta) > 3.0)"))
    process.extraPFStuff.add(process.l1PuppiForMET)
    monitorPerf("L1PuppiForMET", "l1PuppiForMET")

if True:
    process.load("L1Trigger.Phase2L1ParticleFlow.l1ParticleFlow_cff")
    process.pfClustersFromHGC3DClustersEM.src = "hgcalBackEndLayer2ProducerSTC:HGCalBackendLayer2Processor3DClustering"
    process.extraPFStuff.add(
        process.pfClustersFromHGC3DClustersEM, 
        process.pfClustersFromCombinedCalo, 
        process.l1pfProducer
    )
    monitorPerf("L1OldCalo", "l1pfProducer:Calo", makeRespSplit = False)
    monitorPerf("L1OldPF", "l1pfProducer:PF")
    monitorPerf("L1OldPuppi", "l1pfProducer:Puppi")
    if True:
        process.l1OldPuppiForMET = cms.EDFilter("L1TPFCandSelector", 
                src = cms.InputTag("l1pfProducer:Puppi"), 
                cut = cms.string("charge != 0 || abs(eta) < 1.5 || (pt > 20 && abs(eta) < 2.5) || (pt > 40 && 2.5 <= abs(eta) <= 2.85) || (pt > 30 && abs(eta) > 3.0)"))
        process.extraPFStuff.add(process.l1OldPuppiForMET)
        monitorPerf("L1OldPuppiForMET", "l1OldPuppiForMET")

if True: # track quantities with different quality cuts
    process.l1tkv5Stubs = cms.EDFilter("L1TPFCandSelector", src = cms.InputTag("l1pfCandidates:TKVtx"), cut = cms.string("pfTrack.nStubs >= 5"))
    process.l1tkv6Stubs = cms.EDFilter("L1TPFCandSelector", src = cms.InputTag("l1pfCandidates:TKVtx"), cut = cms.string("pfTrack.nStubs >= 6"))
    process.extraPFStuff.add(process.l1tkv5Stubs, process.l1tkv6Stubs)
    monitorPerf("L1TKV5", "l1tkv5Stubs", makeRespSplit = False)
    monitorPerf("L1TKV6", "l1tkv6Stubs", makeRespSplit = False)

process.runPF.associate(process.extraPFStuff)
process.p = cms.Path(process.runPF + process.ntuple + process.l1pfjetTable + process.l1pfmetTable)
process.TFileService = cms.Service("TFileService", fileName = cms.string("perfTuple.root"))

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
if True: # test also PF leptons
    process.ntuple.objects.L1PFMuon = cms.VInputTag("l1pfCandidates:PF",)
    process.ntuple.objects.L1PFMuon_sel = cms.string("abs(pdgId) == 13")
    process.ntuple.objects.L1PFElectron = cms.VInputTag("l1pfCandidates:PF",)
    process.ntuple.objects.L1PFElectron_sel = cms.string("abs(pdgId) == 11")
def respOnly():
    process.p.remove(process.l1pfjetTable)
    process.p.remove(process.l1pfmetTable)
    process.end.remove(process.outnano)
def addCalib():
    process.pfClustersFromL1EGClustersRaw    = process.pfClustersFromL1EGClusters.clone(corrector = "")
    process.pfClustersFromHGC3DClustersRaw   = process.pfClustersFromHGC3DClusters.clone(corrector = "")
    process.pfClustersFromHGC3DClustersEMRaw = process.pfClustersFromHGC3DClustersRaw.clone(emOnly = True, etMin = 0.)
    process.extraPFStuff.add(
            process.pfClustersFromL1EGClustersRaw, 
            process.pfClustersFromHGC3DClustersRaw, 
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
    if hasattr(process.ntuple.objects, 'L1OldPF'):
        process.pfClustersFromHGC3DClustersOldEMRaw = process.pfClustersFromHGC3DClustersEM.clone(corrector = "")
        process.extraPFStuff.add(process.pfClustersFromHGC3DClustersOldEMRaw)
        process.ntuple.objects.L1OldRawHGCalEM = cms.VInputTag('pfClustersFromHGC3DClustersOldEMRaw', )
        process.ntuple.objects.L1OldHGCalEM = cms.VInputTag('pfClustersFromHGC3DClustersEM', )
        process.ntuple.objects.L1OldRawCalo = cms.VInputTag('pfClustersFromCombinedCalo:uncalibrated')
        process.ntuple.objects.L1OldRawCaloEM = cms.VInputTag('pfClustersFromCombinedCalo:emUncalibrated')
        process.ntuple.objects.L1OldRawEcal = cms.VInputTag('pfClustersFromL1EGClustersRaw', 'pfClustersFromHGC3DClustersEMRaw')
        process.ntuple.objects.L1OldEcal = cms.VInputTag(cms.InputTag('l1pfProducer','EmCalo'))
    if False: # extra calibrations for TC clusters
        process.hgc3DClustersNoNoiseTC = process.hgc3DClustersNoNoise.clone(src = "hgcalBackEndLayer2ProducerTC:HGCalBackendLayer2Processor3DClustering")
        process.pfClustersFromHGC3DClustersRawTC = process.pfClustersFromHGC3DClustersRaw.clone(src = "hgc3DClustersNoNoiseTC")
        process.pfClustersFromHGC3DClustersTC    = process.pfClustersFromHGC3DClustersRawTC.clone(corrector =  "L1Trigger/Phase2L1ParticleFlow/data/hadcorr_HGCal3D_TC.root")
        process.pfClustersFromHGC3DClustersEMRawTC = process.pfClustersFromHGC3DClustersRawTC.clone(emOnly = True, etMin = 0.)
        process.extraPFStuff.add(process.hgc3DClustersNoNoiseTC, process.pfClustersFromHGC3DClustersRawTC, process.pfClustersFromHGC3DClustersTC, process.pfClustersFromHGC3DClustersEMRawTC)
        process.ntuple.objects.L1RawHGCalEMTC = cms.VInputTag('pfClustersFromHGC3DClustersEMRawTC')
        process.ntuple.objects.L1RawHGCalTC = cms.VInputTag('pfClustersFromHGC3DClustersRawTC')
        process.ntuple.objects.L1HGCalTC    = cms.VInputTag('pfClustersFromHGC3DClustersTC') 

def goGun(calib=1):
    process.ntuple.isParticleGun = True
    respOnly()
    if calib: 
        addCalib()
def goMT(nthreads=2):
    process.options.numberOfThreads = cms.untracked.uint32(nthreads)
    process.options.numberOfStreams = cms.untracked.uint32(0)
def goOld():
    process.pfClustersFromL1EGClusters.corrector  = "L1Trigger/Phase2L1ParticleFlow/data/emcorr_barrel_93X.root"
    process.pfClustersFromHGC3DClustersEM.corrector =  "L1Trigger/Phase2L1ParticleFlow/data/emcorr_hgc_old3d_93X.root"
    process.pfClustersFromCombinedCalo.hadCorrector =  "L1Trigger/Phase2L1ParticleFlow/data/hadcorr_93X.root"
    process.pfClustersFromHGC3DClusters.corrector =  "L1Trigger/Phase2L1ParticleFlow/data/hadcorr_HGCal3D_STC_93X.root"
    process.pfClustersFromCombinedCaloHCal.hadCorrector =  "L1Trigger/Phase2L1ParticleFlow/data/hadcorr_barrel_93X.root"
    process.pfClustersFromCombinedCaloHF.hadCorrector =  "L1Trigger/Phase2L1ParticleFlow/data/hfcorr_93X.root"


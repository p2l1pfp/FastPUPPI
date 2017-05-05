import FWCore.ParameterSet.Config as cms

process = cms.Process("RESP")

process.load('Configuration.StandardSequences.Services_cff')
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True), allowUnscheduled = cms.untracked.bool(False) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1))
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:l1pf_out.root')
)
process.source.duplicateCheckMode = cms.untracked.string("noDuplicateCheck")

process.ntuple = cms.EDAnalyzer("ResponseNTuplizer",
    genJets = cms.InputTag("ak4GenJetsNoNu"),
    genParticles = cms.InputTag("genParticles"),
    isParticleGun = cms.bool(False),
    objects = cms.PSet(
        # -- offline inputs --
        Ecal = cms.VInputTag('l1tPFEcalProducerFromOfflineRechits:towers','l1tPFHGCalEEProducerFromOfflineRechits:towers', 'l1tPFHFProducerFromOfflineRechits:towers'),
        Hcal = cms.VInputTag('l1tPFHcalProducerFromOfflineRechits:towers','l1tPFHGCalFHProducerFromOfflineRechits:towers', 'l1tPFHGCalBHProducerFromOfflineRechits:towers', 'l1tPFHFProducerFromOfflineRechits:towers'),
        Calo = cms.VInputTag('l1tPFEcalProducerFromOfflineRechits:towers','l1tPFHGCalEEProducerFromOfflineRechits:towers', 'l1tPFHcalProducerFromOfflineRechits:towers', 'l1tPFHGCalFHProducerFromOfflineRechits:towers', 'l1tPFHGCalBHProducerFromOfflineRechits:towers', 'l1tPFHFProducerFromOfflineRechits:towers'),
        TK   = cms.VInputTag('l1tPFTkProducersFromOfflineTracksStrips'),
        # -- TP inputs --
        #TPEcal = cms.VInputTag('l1tPFEcalProducerFromTPDigis:towers','l1tPFHGCalProducerFromTriggerCells:towersEE',),
        #TPHcal = cms.VInputTag('l1tPFHcalProducerFromTPDigis','l1tPFHGCalProducerFromTriggerCells:towersFHBH',),
        #TPTK   = cms.VInputTag('l1tPFTkProducersFromL1Tracks',),
        # -- processed --
        L1RawCalo = cms.VInputTag("InfoOut:RawCalo",),
        L1Calo = cms.VInputTag("InfoOut:Calo",),
        L1TK = cms.VInputTag("InfoOut:TK",),
        L1PF = cms.VInputTag("InfoOut:PF",),
        L1Puppi = cms.VInputTag("InfoOut:Puppi",),
        # -- processed (integer math) --
        L1ICalo = cms.VInputTag("InfoOut:L1Calo",),
        L1ITK = cms.VInputTag("InfoOut:L1TK",),
        L1IPF = cms.VInputTag("InfoOut:L1PF",),
        L1IPuppi = cms.VInputTag("InfoOut:L1Puppi",),
       ## -- clustered --
       #L1ak4RawCalo = cms.VInputTag("ak4L1RawCalo",),
       #L1ak4Calo = cms.VInputTag("ak4L1Calo",),
       #L1ak4TK = cms.VInputTag("ak4L1TK",),
       #L1ak4PF = cms.VInputTag("ak4L1PF",),
       #L1ak4Puppi = cms.VInputTag("ak4L1Puppi",),
    ),
    copyUInts = cms.VInputTag(
        "InfoOut:totNL1Calo", "InfoOut:totNL1TK", "InfoOut:totNL1Mu", "InfoOut:totNL1PF", "InfoOut:totNL1Puppi",
        "InfoOut:maxNL1Calo", "InfoOut:maxNL1TK", "InfoOut:maxNL1Mu", "InfoOut:maxNL1PF", "InfoOut:maxNL1Puppi",
    )
)
process.p = cms.Path(process.ntuple)
process.TFileService = cms.Service("TFileService", fileName = cms.string("respTupleNew.root"))


if True:
    process.load('FastPUPPI.NtupleProducer.caloNtupleProducer_cfi')
    process.load('FastPUPPI.NtupleProducer.ntupleProducer_cfi')
    process.load('FastPUPPI.NtupleProducer.l1tPFMuProducerFromL1Mu_cfi')
    process.CaloInfoOut.outputName = ""; # turn off Ntuples
    process.InfoOut.outputName = ""; # turn off Ntuples
    process.p = cms.Path(process.CaloInfoOut + process.InfoOut + process.ntuple)
def goGun():
    process.ntuple.isParticleGun = True
def goSpring17():
    del process.ntuple.objects.TK
    process.ntuple.objects.TPEcal = cms.VInputTag('l1tPFEcalProducerFromTPDigis:towers', 'l1tPFHGCalProducerFromTriggerCells:towersEE',)
    process.ntuple.objects.TPHcal = cms.VInputTag('l1tPFHcalProducerFromTPDigis', 'l1tPFHGCalProducerFromTriggerCells:towersFHBH',)
    process.ntuple.objects.TPTK   = cms.VInputTag('l1tPFTkProducersFromL1Tracks',)
    if hasattr(process, 'InfoOut'):
        process.InfoOut.L1TrackTag = 'l1tPFTkProducersFromL1Tracks'
        process.CaloInfoOut.EcalTPTags = [ 'l1tPFEcalProducerFromTPDigis:towers', 'l1tPFHGCalProducerFromTriggerCells:towersEE' ]
        process.CaloInfoOut.HcalTPTags = [ 'l1tPFHcalProducerFromTPDigis', 'l1tPFHGCalProducerFromTriggerCells:towersFHBH', ]
        process.CaloInfoOut.caloClusterer.linker.useCorrectedEcal = False
def goRegional(inParallel=False):
    regions = cms.VPSet(
            cms.PSet(
                etaBoundaries = cms.vdouble(-5.5,-4,-3),
                phiSlices = cms.uint32(4),
                etaExtra = cms.double(0.2),
                phiExtra = cms.double(0.2),
            ),
            cms.PSet(
                etaBoundaries = cms.vdouble(-3,-1.5,0,1.5,3),
                phiSlices = cms.uint32(6),
                etaExtra = cms.double(0.2),
                phiExtra = cms.double(0.2),
            ),
            cms.PSet(
                etaBoundaries = cms.vdouble(3,4,5.5),
                phiSlices = cms.uint32(4),
                etaExtra = cms.double(0.2),
                phiExtra = cms.double(0.2),
            ),
    )
    if inParallel:
        process.InfoOutReg = process.InfoOut.clone(regions = regions)
        process.p = cms.Path(process.InfoOut + process.InfoOutReg + process.ntuple)
    else:
        process.InfoOut.regions = regions
if False:
    process.out = cms.OutputModule("PoolOutputModule",
            fileName = cms.untracked.string("l1pf_remade.root"),
    )
    process.e = cms.EndPath(process.out)
if False:
    process.MessageLogger.cerr.FwkReport.reportEvery = 1
    process.maxEvents.input = 3
    process.InfoOut.debug = 2
    if False:
        process.filter = cms.EDFilter("CandViewSelector",
            src = cms.InputTag("genParticles"),
            cut = cms.string("pt > 10 && abs(eta) < 1.5"),
            filter = cms.bool(True),
        )
        process.p = cms.Path(process.filter + process.InfoOut)

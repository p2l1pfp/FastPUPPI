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
        #TK   = cms.VInputTag('l1tPFTkProducersFromOfflineTracksStrips'),
        # -- TP inputs --
        TPEcal = cms.VInputTag('l1tPFEcalProducerFromTPDigis:towers','l1tPFHGCalProducerFromTriggerCells:towersEE',),
        TPHcal = cms.VInputTag('l1tPFHcalProducerFromTPDigis','l1tPFHGCalProducerFromTriggerCells:towersFHBH',),
        TPCalo = cms.VInputTag('l1tPFEcalProducerFromTPDigis:towers', 'l1tPFHGCalProducerFromTriggerCells:towersEE', 'l1tPFHcalProducerFromTPDigis', 'l1tPFHGCalProducerFromTriggerCells:towersFHBH',),
        TPTK   = cms.VInputTag('l1tPFTkProducersFromL1Tracks',),
        # -- processed --
        L1RawEcal = cms.VInputTag(cms.InputTag('CaloInfoOut','emUncalibrated')),
        L1Ecal = cms.VInputTag(cms.InputTag('CaloInfoOut','emCalibrated')),
        L1RawCalo = cms.VInputTag(cms.InputTag('CaloInfoOut','uncalibrated')),
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
    process.CaloInfoOut.outputName = ""; # turn off Ntuples
    process.InfoOut.outputName = ""; # turn off Ntuples
    process.p = cms.Path(process.CaloInfoOut + process.InfoOut + process.ntuple)
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
    process.p = cms.Path(process.genInAcceptance + process.p._seq)
def goGun():
    process.ntuple.isParticleGun = True
def tmpCalib():
    process.CaloInfoOut.caloClusterer.linker.useCorrectedEcal = True
    process.CaloInfoOut.simpleCorrEm = cms.PSet(
                etaBins = cms.vdouble( 0.500,  1.000,  1.500,  2.000,  2.500,  3.000),
                offset  = cms.vdouble(-1.402, -1.733, -2.007, -0.983, -1.140, -1.362),
                scale   = cms.vdouble( 0.977,  0.976,  0.960,  0.915,  0.949,  0.986),
                )
    process.CaloInfoOut.simpleCorrHad = cms.PSet(
            etaBins = cms.vdouble( 0.500,  0.500,  0.500,  1.000,  1.000,  1.000,  1.500,  1.500,  1.500,  2.000,  2.000,  2.000,  2.500,  2.500,  2.500,  3.000,  3.000,  3.000,  3.500,  4.000,  4.500,  5.000),
            emfBins = cms.vdouble( 0.125,  0.500,  0.875,  0.125,  0.500,  0.875,  0.125,  0.500,  0.875,  0.125,  0.500,  0.875,  0.125,  0.500,  0.875,  0.125,  0.500,  0.875,  1.100,  1.100,  1.100,  1.100),
            offset  = cms.vdouble(-2.765, -1.101, -3.387, -3.069, -0.775, -2.698, -5.154,  0.823, -1.693, -2.871, -0.408, -1.276, -2.205, -0.923, -2.050, -3.338,  0.284, -1.839, -3.910, -3.679, -3.361, -4.131),
            scale   = cms.vdouble( 0.951,  0.719,  0.721,  0.977,  0.702,  0.722,  0.915,  0.651,  0.647,  0.586,  0.671,  0.722,  0.608,  0.670,  0.732,  0.544,  0.578,  0.674,  1.157,  1.154,  1.060,  0.744),
            )
def tmpResol():
    process.InfoOut.simpleResolHad = cms.PSet(
            etaBins = cms.vdouble( 1.300,  1.700,  2.800,  3.200,  4.000,  5.000),
            offset  = cms.vdouble( 3.522,  0.078,  2.071,  1.708,  1.148, -0.265),
            scale   = cms.vdouble( 0.124,  0.494,  0.183,  0.257,  0.162,  0.428),
            kind    = cms.string('calo'),
            )
    process.InfoOut.simpleResolEm = cms.PSet(
            etaBins = cms.vdouble( 1.300,  1.700,  2.800,  3.200,  4.000,  5.000),
            offset  = cms.vdouble( 0.849,  0.626,  0.157, -1.305,  0.607, -3.985),
            scale   = cms.vdouble( 0.016,  0.097,  0.043,  0.305,  0.142,  0.626),
            kind    = cms.string('calo'),
            )
    process.InfoOut.simpleResolTrk  = cms.PSet(
            etaBins = cms.vdouble( 0.800,  1.200,  1.500,  2.000,  2.500),
            offset  = cms.vdouble( 0.006,  0.010,  0.010,  0.019,  0.027),
            scale   = cms.vdouble( 0.303,  0.465,  1.003,  1.219,  1.518),
            kind    = cms.string('track'),
            )
def newLink():
    process.InfoOut.linking = cms.PSet(
            trackCaloDR = cms.double(0.15),
            trackCaloNSigmaLow = cms.double(2.0),
            trackCaloNSigmaHigh = cms.double(2.0),
            useTrackCaloSigma = cms.bool(True),
            rescaleUnmatchedTrack = cms.bool(False),
            maxInvisiblePt = cms.double(10.0),
            )
def useClusters():
        process.ntuple.objects.TPEcal = cms.VInputTag('l1tPFEcalProducerFromTPDigis:crystals', 'l1tPFHGCalProducerFrom3DTPs',)
        process.ntuple.objects.TPHcal = cms.VInputTag('l1tPFHcalProducerFromTPDigis', 'l1tPFHGCalProducerFrom3DTPs',)
        process.ntuple.objects.TPCalo = cms.VInputTag('l1tPFEcalProducerFromTPDigis:crystals', 'l1tPFHGCalProducerFrom3DTPs', 'l1tPFHcalProducerFromTPDigis' )
        if hasattr(process, 'InfoOut'):
            process.CaloInfoOut.EcalTPTags = [ 'l1tPFEcalProducerFromTPDigis:crystals' ]
            process.CaloInfoOut.HcalTPTags = [ 'l1tPFHcalProducerFromTPDigis' ]
            process.CaloInfoOut.caloClusterer.linker.useCorrectedEcal = False
            process.InfoOut.CaloClusterTags = [ cms.InputTag('CaloInfoOut','uncalibrated'), cms.InputTag('l1tPFHGCalProducerFrom3DTPs') ]
            process.InfoOut.correctCaloEnergies = False # to become True when calibration will be available
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
        process.p = cms.Path(process.CaloInfoOut + process.InfoOut + process.InfoOutReg + process.ntuple)
    else:
        process.InfoOut.regions = regions
if False:
    goGun(); tmpCalib(); tmpResol();
    #process.source.fileNames = ['/store/cmst3/user/gpetrucc/l1phase2/Spring17D/010517/inputs_17D_SinglePion0_NoPU_job42.root']
    process.source.fileNames = ['/store/cmst3/user/gpetrucc/l1phase2/Spring17D/010517/inputs_17D_SinglePion_NoPU_job42.root']
    #process.source.fileNames = ['/store/cmst3/user/gpetrucc/l1phase2/Spring17D/010517/inputs_17D_TTbar_NoPU_job1.root']
    process.out = cms.OutputModule("PoolOutputModule",
            fileName = cms.untracked.string("l1pf_remade.root"),
    )
    process.e = cms.EndPath(process.out)
    process.maxEvents.input = 50
    process.MessageLogger.cerr.FwkReport.reportEvery = 1
    if False:
        #process.source.fileNames = [ 'file:l1pf_remade.root' ]
        process.TFileService.fileName = cms.string("respTupleNew_1.root")
        process.out.fileName = cms.untracked.string("l1pf_remade_1.root")
        process.source.eventsToProcess = cms.untracked.VEventRange("1:6315:307280")
    process.InfoOut.debug = 1

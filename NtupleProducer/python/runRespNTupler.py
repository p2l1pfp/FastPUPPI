import FWCore.ParameterSet.Config as cms

process = cms.Process("RESP")

process.load('Configuration.StandardSequences.Services_cff')
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True), allowUnscheduled = cms.untracked.bool(False) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1))
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:/eos/cms/store/cmst3/user/gpetrucc/l1phase2/Spring17D/200517/inputs_17D_SinglePion_PU0_job42.root')
)
process.source.duplicateCheckMode = cms.untracked.string("noDuplicateCheck")


process.load('FastPUPPI.NtupleProducer.caloNtupleProducer_cfi')
process.load('FastPUPPI.NtupleProducer.ntupleProducer_cfi')
process.CaloInfoOut.outputName = ""; # turn off Ntuples
process.InfoOut.outputName = ""; # turn off Ntuples

process.ntuple = cms.EDAnalyzer("ResponseNTuplizer",
    genJets = cms.InputTag("ak4GenJetsNoNu"),
    genParticles = cms.InputTag("genParticles"),
    isParticleGun = cms.bool(False),
    doRandom = cms.bool(False),
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
        L1TKV = cms.VInputTag("InfoOut:TKVtx",),
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
process.p = cms.Path(process.CaloInfoOut + process.InfoOut + process.ntuple)
process.TFileService = cms.Service("TFileService", fileName = cms.string("respTupleNew.root"))

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
    process.chGenInAcceptance = cms.EDFilter("GenParticleSelector",
        src = cms.InputTag("genParticles"),
        cut = cms.string("status == 1 && (abs(eta) < 2.5 && pt > 2 && charge != 0)")
    )
    process.phGenInAcceptance = cms.EDFilter("GenParticleSelector",
        src = cms.InputTag("genParticles"),
        cut = cms.string("status == 1 && pt > 1 && pdgId == 22")
    )
    process.ntuple.objects.ChGenAcc = cms.VInputTag(cms.InputTag("chGenInAcceptance"))
    process.ntuple.objects.PhGenAcc = cms.VInputTag(cms.InputTag("phGenInAcceptance"))
    process.p = cms.Path(process.genInAcceptance + process.chGenInAcceptance + process.phGenInAcceptance + process.p._seq)
    process.ntuple.objects.L1PFCharged = cms.VInputTag("InfoOut:PF",)
    process.ntuple.objects.L1PFCharged_sel = cms.string("charge != 0")
    process.ntuple.objects.L1PFPhoton = cms.VInputTag("InfoOut:PF",)
    process.ntuple.objects.L1PFPhoton_sel = cms.string("pdgId == 22")
    process.ntuple.objects.L1IPFCharged = cms.VInputTag("InfoOut:L1PF",)
    process.ntuple.objects.L1IPFCharged_sel = cms.string("charge != 0")
    process.ntuple.objects.L1IPFPhoton = cms.VInputTag("InfoOut:L1PF",)
    process.ntuple.objects.L1IPFPhoton_sel = cms.string("pdgId == 22")
def goGun():
    process.ntuple.isParticleGun = True
def goRandom():
    process.ntuple.doRandom = True
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
def haveFun():
    process.l1tPFHGCalProducerFrom3DTPsEM = cms.EDProducer('GetEMPart',
            src = cms.InputTag('l1tPFHGCalProducerFrom3DTPs'),
    )
    process.p.replace(process.CaloInfoOut, process.l1tPFHGCalProducerFrom3DTPsEM + process.CaloInfoOut)
    process.ntuple.objects.TPEcal  = [ cms.InputTag('l1tPFHGCalProducerFrom3DTPsEM'),cms.InputTag('l1tPFEcalProducerFromL1EGCrystalClusters'), ]
    process.CaloInfoOut.EcalTPTags = [ cms.InputTag('l1tPFHGCalProducerFrom3DTPsEM'),cms.InputTag('l1tPFEcalProducerFromL1EGCrystalClusters'), ]
    process.InfoOut.CaloClusterTags = [ cms.InputTag('CaloInfoOut','uncalibrated'), ]
    process.InfoOut.EmClusterTags   = [ cms.InputTag('l1tPFHGCalProducerFrom3DTPsEM'),cms.InputTag('l1tPFEcalProducerFromL1EGCrystalClusters'), ]
    process.InfoOut.linking.altAlgo = cms.string("PFAlgo3")
    process.InfoOut.linking.trackCaloNSigmaLow  =  2.0
    process.InfoOut.linking.trackCaloNSigmaHigh =  1.2
    process.InfoOut.linking.trackCaloDR = cms.double(0.12) 
    process.InfoOut.linking.trackEmDR   = cms.double(0.04) # 1 Ecal crystal size is 0.02, and ~2 cm in HGCal is ~0.007
    process.InfoOut.linking.emCaloDR    = cms.double(0.10) # 1 Hcal tower size is ~0.09
    process.InfoOut.linking.trackEmPtMinFrac = cms.double(0.5) # Calo object must have an EM Et at least half of that of the EM cluster to allow linking
                                                               # currently it's likely to be always >= 1, but it may be different with different calibrations
    process.InfoOut.linking.caloReLink  = cms.bool(False)
    process.InfoOut.linking.caloReLinkDR = cms.double(0.3)
    process.InfoOut.linking.caloReLinkThreshold = cms.double(0.5)
    process.InfoOut.linking.maxInvisiblePt      = 10.0
    process.InfoOut.linking.tightTrackMinStubs = cms.uint32(6)
    process.InfoOut.linking.tightTrackMaxChi2  = cms.double(50)
    process.InfoOut.linking.tightTrackMaxInvisiblePt = cms.double(20)
    process.InfoOut.linking.sumTkCaloErr2 = cms.bool(True) # add up track calo errors in quadrature instead of linearly
    process.InfoOut.linking.ecalPriority  = cms.bool(True)
    process.ntuple.objects.AltEcal = cms.VInputTag("InfoOut:L1EmCalo",)
    process.ntuple.objects.AltPF = cms.VInputTag("InfoOut:AltPF",)
    process.ntuple.objects.AltPFCharged = cms.VInputTag("InfoOut:AltPF",)
    process.ntuple.objects.AltPFPhoton  = cms.VInputTag("InfoOut:AltPF",)
    process.ntuple.objects.AltPFCharged_sel = cms.string("charge != 0")
    process.ntuple.objects.AltPFPhoton_sel  = cms.string("pdgId == 22")
    process.ntuple.objects.AltPuppi = cms.VInputTag("InfoOut:AltPuppi",)
    process.ntuple.objects.AltPuppiCharged = cms.VInputTag("InfoOut:AltPuppi",)
    process.ntuple.objects.AltPuppiCharged_sel = cms.string("charge != 0")
    process.ntuple.objects.AltPFChargedScaled = cms.VInputTag("InfoOut:AltPF",)
    process.ntuple.objects.AltPFChargedScaled_sel = cms.string("charge != 0 && status == 2")
    process.ntuple.objects.AltPFDiscTrack = cms.VInputTag("InfoOut:AltPFDiscarded",)
    process.ntuple.objects.AltPFDiscCaloT = cms.VInputTag("InfoOut:AltPFDiscarded",)
    process.ntuple.objects.AltPFDiscCaloG = cms.VInputTag("InfoOut:AltPFDiscarded",)
    process.ntuple.objects.AltPFDiscEm    = cms.VInputTag("InfoOut:AltPFDiscarded",)
    process.ntuple.objects.AltPFDiscTrack_sel = cms.string("charge != 0 && status == 1")
    process.ntuple.objects.AltPFDiscCaloT_sel = cms.string("charge == 0 && status == 0")
    process.ntuple.objects.AltPFDiscCaloG_sel = cms.string("charge == 0 && status == 1")
    process.ntuple.objects.AltPFDiscEm_sel    = cms.string("charge == 0 && status == 2")
    for M in process.InfoOut, process.CaloInfoOut:
        if hasattr(M, 'simpleCorrEm'):  del M.simpleCorrEm
        if hasattr(M, 'simpleCorrHad'): del M.simpleCorrHad
        #M.correctorEmfMax  = 0.875
        #M.correctorEmfBins = 9
        M.correctorEmfMax  = 1.0
        M.correctorEmfBins = 11
    process.InfoOut.correctCaloEnergies = True
    process.InfoOut.rawEmCaloPtMin = cms.double(0.5)
    process.InfoOut.rawCaloPtMin   = cms.double(1.0)
def addBH():
    process.ntuple.objects.TPHcal += [ 'l1tPFHGCalBHProducerFromOfflineRechits:towers' ] 
    process.ntuple.objects.TPCalo += [ 'l1tPFHGCalBHProducerFromOfflineRechits:towers' ] 
    process.CaloInfoOut.HcalTPTags += [ 'l1tPFHGCalBHProducerFromOfflineRechits:towers' ]
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
    process.InfoOut.fillTrackTree = cms.untracked.int32(1)
    process.source.fileNames = ['file:/eos/cms/store/cmst3/user/gpetrucc/l1phase2/Spring17D/200517/inputs_17D_TTbar_PU0_job1.root' ]
    process.p.remove(process.ntuple)
    process.TFileService.fileName = cms.string("trackTupleNew.root")
if False:
    process.CaloInfoOutBackup = process.CaloInfoOut.clone()
    process.InfoOutBackup = process.InfoOut.clone()
    process.p.replace(process.CaloInfoOut, process.CaloInfoOutBackup + process.CaloInfoOut)
    process.p.replace(process.InfoOut, process.InfoOutBackup + process.InfoOut)
    #goGun(); 
    haveFun(); 
    addBH();
    #tmpCalib(); tmpResol(); newLink();
    #process.source.fileNames = ['file:/eos/cms/store/cmst3/user/gpetrucc/l1phase2/Spring17D/200517/inputs_17D_SinglePion0_PU0_job42.root']
    #process.source.fileNames = ['/store/cmst3/user/gpetrucc/l1phase2/Spring17D/200517/inputs_17D_SinglePion_PU0_job42.root']
    process.source.fileNames = ['file:/eos/cms/store/cmst3/user/gpetrucc/l1phase2/Spring17D/200517/inputs_17D_TTbar_PU0_job1.root']
    #process.source.fileNames = ['/store/cmst3/user/gpetrucc/l1phase2/Spring17D/200517/inputs_17D_SingleTau_PU0_job42.root']
    process.out = cms.OutputModule("PoolOutputModule",
            fileName = cms.untracked.string("l1pf_remade.root"),
    )
    process.e = cms.EndPath(process.out)
    #process.source.skipEvents = cms.untracked.uint32(10)
    process.maxEvents.input = 10
    process.MessageLogger.cerr.FwkReport.reportEvery = 1
    process.InfoOut.altDebug = cms.untracked.int32(1)
    #process.InfoOutNRL = process.InfoOut.clone()
    #process.InfoOutNRL.linking.caloReLink  = cms.bool(False)
    #process.p.replace(process.InfoOut, process.InfoOutNRL + process.InfoOut)
    if True:
        #process.CaloInfoOut.debug = cms.untracked.int32(1)
        process.TFileService.fileName = cms.string("respTupleNew_1.root")
        process.out.fileName = cms.untracked.string("l1pf_remade_1.root")
        process.source.eventsToProcess = cms.untracked.VEventRange("1:71:2999",)
        process.InfoOut.debugEta = cms.untracked.double(-1.1)
        process.InfoOut.debugPhi = cms.untracked.double(+0.7)
        process.InfoOut.debugR   = cms.untracked.double(0.5)
    if False:
        process.InfoOut.regions = cms.VPSet(
                cms.PSet(
                    etaBoundaries = cms.vdouble(-5.5,-3),
                    phiSlices = cms.uint32(1),
                    etaExtra = cms.double(0.2),
                    phiExtra = cms.double(0.2),
                ),
                cms.PSet(
                    etaBoundaries = cms.vdouble(-3,-1.3,1.3,3),
                    phiSlices = cms.uint32(3),
                    etaExtra = cms.double(0.2),
                    phiExtra = cms.double(0.2),
                ),
                cms.PSet(
                    etaBoundaries = cms.vdouble(3,5.5),
                    phiSlices = cms.uint32(1),
                    etaExtra = cms.double(0.2),
                    phiExtra = cms.double(0.2),
                ),
        )


import FWCore.ParameterSet.Config as cms

process = cms.Process("RESP")

process.load('Configuration.StandardSequences.Services_cff')
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True), allowUnscheduled = cms.untracked.bool(False) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1))
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:/eos/cms/store/cmst3/user/jngadiub/L1PFInputs/SingleNeutrino_PU140/inputs_17D_SingleNeutrino_PU140_job1.root'),
)
process.source.duplicateCheckMode = cms.untracked.string("noDuplicateCheck")


process.load('FastPUPPI.NtupleProducer.l1tPFHGCalProducerFrom3DTPsEM_cfi')
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
        #TPEcal = cms.VInputTag('l1tPFEcalProducerFromTPDigis:towers','l1tPFHGCalProducerFromTriggerCells:towersEE',),
        TPEcal = cms.VInputTag('l1tPFHGCalProducerFrom3DTPsEM', 'l1tPFEcalProducerFromL1EGCrystalClusters', ),
        TPHcal = cms.VInputTag('l1tPFHcalProducerFromTPDigis', 'l1tPFHGCalProducerFromTriggerCells:towersFHBH', 'l1tPFHGCalBHProducerFromOfflineRechits:towers'),
        TPCalo = cms.VInputTag('l1tPFHGCalProducerFrom3DTPsEM', 'l1tPFEcalProducerFromL1EGCrystalClusters', 'l1tPFHcalProducerFromTPDigis', 'l1tPFHGCalProducerFromTriggerCells:towersFHBH', 'l1tPFHGCalBHProducerFromOfflineRechits:towers'),
        TPTK   = cms.VInputTag('l1tPFTkProducersFromL1Tracks',),
        # -- processed --
        #L1RawEcal = cms.VInputTag(cms.InputTag('CaloInfoOut','emUncalibrated')),
        #L1RawCalo = cms.VInputTag(cms.InputTag('CaloInfoOut','uncalibrated')),
        #L1Ecal = cms.VInputTag(cms.InputTag('CaloInfoOut','emCalibrated')),
        L1RawEcal = cms.VInputTag(cms.InputTag('InfoOut','RawEmCalo')),
        L1RawCalo = cms.VInputTag(cms.InputTag('InfoOut','RawCalo')),
        L1Ecal = cms.VInputTag(cms.InputTag('InfoOut','EmCalo')),
        L1Calo = cms.VInputTag("InfoOut:Calo",),
        L1TK = cms.VInputTag("InfoOut:TK",),
        L1TKV = cms.VInputTag("InfoOut:TKVtx",),
        L1PF = cms.VInputTag("InfoOut:PF",),
        L1Puppi = cms.VInputTag("InfoOut:Puppi",),
    ),
    copyUInts = cms.VInputTag(),
)
for X in "tot","max":
    for I in "Calo EmCalo TK Mu".split(): process.ntuple.copyUInts.append( "InfoOut:%sNL1%s" % (X,I))
    for O in [""] + "Charged Neutral ChargedHadron NeutralHadron Photon Electron Muon".split():
        process.ntuple.copyUInts.append( "InfoOut:%sNL1PF%s" % (X,O))
        process.ntuple.copyUInts.append( "InfoOut:%sNL1Puppi%s" % (X,O))
process.p = cms.Path(process.l1tPFHGCalProducerFrom3DTPsEM + process.CaloInfoOut + process.InfoOut + process.ntuple)
process.TFileService = cms.Service("TFileService", fileName = cms.string("respTuple.root"))

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
def goGun():
    process.ntuple.isParticleGun = True
def goRandom():
    process.ntuple.doRandom = True
if False:
    process.ntuple.objects.PFChargedScaled = cms.VInputTag("InfoOut:PF",)
    process.ntuple.objects.PFChargedScaled_sel = cms.string("charge != 0 && status == 2")
    process.ntuple.objects.PFDiscTrack = cms.VInputTag("InfoOut:PFDiscarded",)
    process.ntuple.objects.PFDiscCaloT = cms.VInputTag("InfoOut:PFDiscarded",)
    process.ntuple.objects.PFDiscCaloG = cms.VInputTag("InfoOut:PFDiscarded",)
    process.ntuple.objects.PFDiscEm    = cms.VInputTag("InfoOut:PFDiscarded",)
    process.ntuple.objects.PFDiscTrack_sel = cms.string("charge != 0 && status == 1")
    process.ntuple.objects.PFDiscCaloT_sel = cms.string("charge == 0 && status == 0")
    process.ntuple.objects.PFDiscCaloG_sel = cms.string("charge == 0 && status == 1")
    process.ntuple.objects.PFDiscEm_sel    = cms.string("charge == 0 && status == 2")
def goRegional(inParallel=False,mode="atCalo"):
    regions = cms.VPSet(
            cms.PSet(
                etaBoundaries = cms.vdouble(-5.5,-4,-3),
                phiSlices = cms.uint32(4),
                etaExtra = cms.double(0.25),
                phiExtra = cms.double(0.25),
            ),
            cms.PSet(
                etaBoundaries = cms.vdouble(-3,-1.5,-0.5,0.5,1.5,3),
                phiSlices = cms.uint32(6),
                etaExtra = cms.double(0.25),
                phiExtra = cms.double(0.25),
            ),
            cms.PSet(
                etaBoundaries = cms.vdouble(3,4,5.5),
                phiSlices = cms.uint32(4),
                etaExtra = cms.double(0.25),
                phiExtra = cms.double(0.25),
            ),
    )
    if inParallel:
        process.InfoOutReg = process.InfoOut.clone(regions = regions)
        process.p = cms.Path(process.CaloInfoOut + process.InfoOut + process.InfoOutReg + process.ntuple)
        process.InfoOut.useRelativeRegionalCoordinates = cms.bool(False)
    else:
        process.InfoOut.regions = regions
        process.InfoOut.useRelativeRegionalCoordinates = cms.bool(True)
        process.InfoOut.trackRegionMode = cms.string(mode)
def gbr(neta,nphi,etaex=0.3,phiex=0.2,mode="any"):
    regions = cms.VPSet(
            cms.PSet(
                etaBoundaries = cms.vdouble(*[(-1.5+3*i/neta) for i in xrange(neta+1)]),
                phiSlices = cms.uint32(nphi),
                etaExtra = cms.double(etaex),
                phiExtra = cms.double(phiex),
            ),
    )
    process.InfoOut.regions = regions
    process.InfoOut.useRelativeRegionalCoordinates = cms.bool(True)
    process.InfoOut.trackRegionMode = cms.string(mode)
def gbrExt(nphi,mode="atCalo"): ## take the barrel calo but include the full tracker
    process.CaloInfoOut.EcalTPTags = [ cms.InputTag('l1tPFEcalProducerFromL1EGCrystalClusters') ]
    process.CaloInfoOut.HcalTPTags = [ cms.InputTag('l1tPFHcalProducerFromTPDigis') ]
    process.InfoOut.EmClusterTags = [ cms.InputTag('l1tPFEcalProducerFromL1EGCrystalClusters') ]
    regions = cms.VPSet(
            cms.PSet(
                etaBoundaries = cms.vdouble([-3,3]),
                phiSlices = cms.uint32(nphi),
                etaExtra = cms.double(0),
                phiExtra = cms.double(0),
            ),
    )
    process.InfoOut.regions = regions
    process.InfoOut.useRelativeRegionalCoordinates = cms.bool(True)
    process.InfoOut.trackRegionMode = cms.string(mode)
if False:
    process.InfoOut.fillTrackTree = cms.untracked.int32(1)
    process.source.fileNames = ['file:/eos/cms/store/cmst3/user/gpetrucc/l1phase2/Spring17D/200517/inputs_17D_TTbar_PU0_job1.root' ]
    process.p.remove(process.ntuple)
    process.TFileService.fileName = cms.string("trackTupleNew.root")
if False: # run bitwise PF for Vivado
    goRegional()
    process.InfoOut.useRelativeRegionalCoordinates = cms.bool(True)
    process.InfoOut.linking.algo = "BitwisePF"
    process.source.fileNames = ['file:/eos/cms/store/cmst3/user/jngadiub/L1PFInputs/TTbar_PU140/inputs_17D_TTbar_PU140_job1.root'] 
if False: # prepare dump file for Vivado (PF IP core)
    goRegional()
    process.InfoOut.useRelativeRegionalCoordinates = cms.bool(True)
    process.InfoOut.regionDumpFileName = cms.untracked.string("regions_TTbar_PU140.dump")
    process.source.fileNames = ['file:/eos/cms/store/cmst3/user/jngadiub/L1PFInputs/TTbar_PU140/inputs_17D_TTbar_PU140_job1.root'] #, 'file:/eos/cms/store/cmst3/user/jngadiub/L1PFInputs/TTbar_PU140/inputs_17D_TTbar_PU140_job2.root', 'file:/eos/cms/store/cmst3/user/jngadiub/L1PFInputs/TTbar_PU140/inputs_17D_TTbar_PU140_job3.root', 'file:/eos/cms/store/cmst3/user/jngadiub/L1PFInputs/TTbar_PU140/inputs_17D_TTbar_PU140_job4.root', 'file:/eos/cms/store/cmst3/user/jngadiub/L1PFInputs/TTbar_PU140/inputs_17D_TTbar_PU140_job5.root', 'file:/eos/cms/store/cmst3/user/jngadiub/L1PFInputs/TTbar_PU140/inputs_17D_TTbar_PU140_job6.root', ]
if False: # prepare dump file for Vivado (Regionizer)
    gbr(1,12,0,0,"atCalo")
    process.InfoOut.regionDumpFileName = cms.untracked.string("barrel_sectors_1x12_TTbar_PU140.dump")
    process.source.fileNames = ['file:/eos/cms/store/cmst3/user/jngadiub/L1PFInputs/TTbar_PU140/inputs_17D_TTbar_PU140_job1.root'] #, 'file:/eos/cms/store/cmst3/user/jngadiub/L1PFInputs/TTbar_PU140/inputs_17D_TTbar_PU140_job2.root', 'file:/eos/cms/store/cmst3/user/jngadiub/L1PFInputs/TTbar_PU140/inputs_17D_TTbar_PU140_job3.root', 'file:/eos/cms/store/cmst3/user/jngadiub/L1PFInputs/TTbar_PU140/inputs_17D_TTbar_PU140_job4.root', 'file:/eos/cms/store/cmst3/user/jngadiub/L1PFInputs/TTbar_PU140/inputs_17D_TTbar_PU140_job5.root', 'file:/eos/cms/store/cmst3/user/jngadiub/L1PFInputs/TTbar_PU140/inputs_17D_TTbar_PU140_job6.root', ]
if False: # prepare dump file for Vivado (Vertexing)
    gbrExt(12,"atCalo")
    process.InfoOut.regionDumpFileName = cms.untracked.string("barrel_alltracks_sectors_1x12_TTbar_PU140.dump")
    process.source.fileNames = ['file:/eos/cms/store/cmst3/user/jngadiub/L1PFInputs/TTbar_PU140/inputs_17D_TTbar_PU140_job1.root'] #, 'file:/eos/cms/store/cmst3/user/jngadiub/L1PFInputs/TTbar_PU140/inputs_17D_TTbar_PU140_job2.root', 'file:/eos/cms/store/cmst3/user/jngadiub/L1PFInputs/TTbar_PU140/inputs_17D_TTbar_PU140_job3.root', 'file:/eos/cms/store/cmst3/user/jngadiub/L1PFInputs/TTbar_PU140/inputs_17D_TTbar_PU140_job4.root', 'file:/eos/cms/store/cmst3/user/jngadiub/L1PFInputs/TTbar_PU140/inputs_17D_TTbar_PU140_job5.root', 'file:/eos/cms/store/cmst3/user/jngadiub/L1PFInputs/TTbar_PU140/inputs_17D_TTbar_PU140_job6.root', ]

def comp4():
    process.InfoOutGlobal = process.InfoOut.clone()
    goRegional(mode="atCalo")
    process.InfoOutRegional = process.InfoOut.clone()
    process.InfoOut.linking.ecalPriority = False 
    process.InfoOut.linking.trackMuMatch = "drBestByPtDiff" 
    process.InfoOut.linking.trackCaloLinkMetric = "bestByDR2Pt2"
    process.InfoOut.linking.trackEmUseAlsoTrackSigma = False 
    process.InfoOut.linking.emCaloUseAlsoCaloSigma = False
    process.InfoOutSimpler = process.InfoOut.clone()
    process.InfoOut.linking.algo = "BitwisePF"
    process.p.replace(process.InfoOut, process.InfoOutGlobal + process.InfoOutRegional + process.InfoOutSimpler + process.InfoOut)
    process.ntuple.objects.L1PFGlobal = cms.VInputTag("InfoOutGlobal:PF",)
    process.ntuple.objects.L1PFRegional = cms.VInputTag("InfoOutRegional:PF",)
    process.ntuple.objects.L1PFSimpler = cms.VInputTag("InfoOutSimpler:PF",)
if False:
    if 1:
        goRegional(mode="atCalo")
        process.InfoOut.useRelativeRegionalCoordinates = cms.bool(True)
        process.InfoOut.linking.ecalPriority = False 
        process.InfoOut.linking.trackMuMatch = "drBestByPtDiff" 
        process.InfoOut.linking.trackCaloLinkMetric = "bestByDR2Pt2"
        process.InfoOut.linking.trackEmUseAlsoTrackSigma = False 
        process.InfoOut.linking.emCaloUseAlsoCaloSigma = False
        process.InfoOut.altDebug = cms.untracked.int32(1)
        #process.CaloInfoOutBackup = process.CaloInfoOut.clone()
        process.InfoOutBackup = process.InfoOut.clone()
        #process.p.replace(process.CaloInfoOut, process.CaloInfoOutBackup + process.CaloInfoOut)
        process.p.replace(process.InfoOut, process.InfoOutBackup + process.InfoOut)
        process.InfoOut.linking.algo = "BitwisePF"
        process.InfoOut.bitwiseDebug = cms.untracked.int32(1)
        process.ntuple.objects.L1PFBackup = cms.VInputTag("InfoOutBackup:PF",)
    else:
        comp4();
    #process.source.fileNames = ['file:/eos/cms/store/cmst3/user/jngadiub/L1PFInputs/TTbar_PU140/inputs_17D_TTbar_PU140_job1.root'] 
    process.source.fileNames = ['file:/eos/cms/store/cmst3/user/jngadiub/L1PFInputs/TTbar_PU0/inputs_17D_TTbar_PU0_job51.root'] 
    process.out = cms.OutputModule("PoolOutputModule",
            fileName = cms.untracked.string("l1pf_remade.root"),
    )
    process.e = cms.EndPath(process.out)
    #process.source.skipEvents = cms.untracked.uint32(10)
    process.maxEvents.input = 500
    process.MessageLogger.cerr.FwkReport.reportEvery = 1
    #process.InfoOut.debug = cms.untracked.int32(1)
    process.source.eventsToProcess = cms.untracked.VEventRange("1:505:24891",)
    if 0:
        process.InfoOut.altDebug = cms.untracked.int32(1)
        #process.CaloInfoOut.debug = cms.untracked.int32(1)
        process.TFileService.fileName = cms.string("respTupleNew_1.root")
        process.out.fileName = cms.untracked.string("l1pf_remade_1.root")
        process.source.eventsToProcess = cms.untracked.VEventRange("1:979:48261",)
    if 1:
        for X in "InfoOut", "InfoOutBackup",:
            P = getattr(process,X,None)
            if P:
                P.debugEta = cms.untracked.double(+0.58)
                P.debugPhi = cms.untracked.double(+1.87)
                P.debugR   = cms.untracked.double(0.6)

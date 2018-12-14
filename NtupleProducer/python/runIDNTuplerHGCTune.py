import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras

process = cms.Process("ID", eras.Phase2_trigger)

process.load('Configuration.StandardSequences.Services_cff')
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10))
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:inputs.root'),
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck")
)

process.load('Configuration.Geometry.GeometryExtended2023D17Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '100X_upgrade2023_realistic_v1', '')

process.load("L1Trigger.L1THGCal.hgcalTriggerPrimitives_cff")
process.hgcalConcentratorProducerSTC = process.hgcalConcentratorProducer.clone()
process.hgcalConcentratorProducerSTC.ProcessorParameters.Method = 'superTriggerCellSelect'
process.hgcalConcentratorProducerSTC.ProcessorParameters.triggercell_threshold_silicon = 0
process.hgcalConcentratorProducerSTC.ProcessorParameters.TCThreshold_fC = 0

process.hgcalBackEndLayer1ProducerSTC =  process.hgcalBackEndLayer1Producer.clone()
process.hgcalBackEndLayer1ProducerSTC.InputTriggerCells = cms.InputTag("hgcalConcentratorProducerSTC","HGCalConcentratorProcessorSelection")
process.hgcalBackEndLayer1ProducerSTC.ProcessorParameters.C2d_parameters.clusterType = cms.string('dummyC2d')
process.hgcalBackEndLayer1ProducerSTC.ProcessorParameters.C2d_parameters.threshold_scintillator = cms.double(0)
process.hgcalBackEndLayer1ProducerSTC.ProcessorParameters.C2d_parameters.threshold_silicon      = cms.double(0)

process.hgcalBackEndLayer2ProducerSTC = process.hgcalBackEndLayer2Producer.clone()
process.hgcalBackEndLayer2ProducerSTC.ProcessorParameters.C3d_parameters.dR_multicluster = 0.03
process.hgcalBackEndLayer2ProducerSTC.ProcessorParameters.C3d_parameters.type_multicluster = 'HistoMaxC3d'
process.hgcalBackEndLayer2ProducerSTC.ProcessorParameters.C3d_parameters.threshold_histo_multicluster = 20
process.hgcalBackEndLayer2ProducerSTC.InputCluster = cms.InputTag("hgcalBackEndLayer1ProducerSTC","HGCalBackendLayer1Processor2DClustering")
process.hgcalBackEndLayer2ProducerSTC.ProcessorParameters.C3d_parameters.dR_multicluster_byLayer = cms.vdouble(
            [0] + [0.010]*7 + [0.020]*7 + [0.030]*7 + [0.040]*7 +   [0.040]*6 + [0.050]*6  +  [0.050]*12 )

ntuple = cms.EDAnalyzer("IDNTuplizer",
    src = cms.InputTag("FILLME"),
    cut = cms.string("pt > 1"),
    genParticles = cms.InputTag("genParticles"),
    propagateToCalo = cms.bool(True),
    drMax = cms.double(0.1),
    minRecoPtOverGenPt = cms.double(0.3),
    onlyMatched = cms.bool(True),
    variables = cms.PSet(
	pt = cms.string("pt"),
	emPt = cms.string("? hOverE >= 0 ? pt/(1+hOverE) : 0"),
 	eta = cms.string("eta"),
 	phi = cms.string("phi"),
 	hOverE = cms.string("hOverE"),
        emf = cms.string("? hOverE >= 0 ? 1/(1+hOverE) : 0"),
 	showerLength = cms.string("showerLength"),
 	coreShowerLength = cms.string("coreShowerLength"),
 	firstLayer = cms.string("firstLayer"),
 	maxLayer = cms.string("maxLayer"),
 	eMax = cms.string("eMax"),
 	eMaxOverE = cms.string("eMax/energy"),
 	sigmaEtaEtaMax = cms.string("sigmaEtaEtaMax"),
 	sigmaPhiPhiMax = cms.string("sigmaPhiPhiMax"),
 	sigmaEtaEtaTot = cms.string("sigmaEtaEtaTot"),
 	sigmaPhiPhiTot = cms.string("sigmaPhiPhiTot"),
 	sigmaZZ = cms.string("sigmaZZ"),
 	sigmaRRTot = cms.string("sigmaRRTot"),
 	sigmaRRMax = cms.string("sigmaRRMax"),
 	sigmaRRMean = cms.string("sigmaRRMean"),
    ),
)

process.ntupleSTC = ntuple.clone( src = cms.InputTag("hgcalBackEndLayer2ProducerSTC","HGCalBackendLayer2Processor3DClustering") )
  
modules = [
    # process.hgcalConcentratorProducerSTC, # already in the inputs
    process.hgcalBackEndLayer1ProducerSTC,
    process.hgcalBackEndLayer2ProducerSTC,
    process.ntupleSTC
]

def newClustering(postfix,
        concentratorAlgo="thresholdSelect", # 'thresholdSelect' or 'superTriggerCellSelect'
        concentratorThreshold=2,
        layer1Algo="dummyC2d",
        layer1Threshold=0,
        layer2Algo="HistoMaxC3d",
        layer2Threshold=0,
        layer2dR=0.03,
        reuseConc=None, reuseL1=None):
    global modules
    if reuseConc is None:
        conc = process.hgcalConcentratorProducer.clone()
        conc.ProcessorParameters.Method = concentratorAlgo
        conc.ProcessorParameters.triggercell_threshold_silicon = concentratorThreshold
        if concentratorThreshold == 0: conc.ProcessorParameters.TCThreshold_fC = 0
        setattr(process, "hgcalConcentratorProducer"+postfix, conc)
        reuseConc = postfix
        modules.append(conc)
    #
    if reuseL1 is None:
        bel1 = process.hgcalBackEndLayer1Producer.clone()
        bel1.InputTriggerCells = cms.InputTag("hgcalConcentratorProducer"+reuseConc,"HGCalConcentratorProcessorSelection")
        bel1.ProcessorParameters.C2d_parameters.clusterType = cms.string(layer1Algo)
        bel1.ProcessorParameters.C2d_parameters.threshold_scintillator = cms.double(layer1Threshold)
        bel1.ProcessorParameters.C2d_parameters.threshold_silicon      = cms.double(layer1Threshold)
        setattr(process, "hgcalBackEndLayer1Producer"+postfix, bel1)
        modules.append(bel1)
        reuseL1 = postfix
    #
    bel2 = process.hgcalBackEndLayer2Producer.clone()
    bel2.InputCluster = cms.InputTag("hgcalBackEndLayer1Producer"+reuseL1,"HGCalBackendLayer1Processor2DClustering")
    bel2.ProcessorParameters.C3d_parameters.type_multicluster = layer2Algo
    bel2.ProcessorParameters.C3d_parameters.threshold_histo_multicluster = layer2Threshold
    if type(layer2dR) == list:
        bel2.ProcessorParameters.C3d_parameters.dR_multicluster_byLayer = cms.vdouble(layer2dR)
        bel2.ProcessorParameters.C3d_parameters.dR_multicluster = 0.03 if type(layer2dR) == list else layer2dR
    else:
        bel2.ProcessorParameters.C3d_parameters.dR_multicluster = 0.03 if type(layer2dR) == list else layer2dR
    setattr(process, "hgcalBackEndLayer2Producer"+postfix, bel2)
    modules.append(bel2)
    #
    nt = ntuple.clone()
    nt.src = cms.InputTag("hgcalBackEndLayer2Producer"+postfix,"HGCalBackendLayer2Processor3DClustering")
    setattr(process, "ntuple"+postfix, nt)
    modules.append(nt)


newClustering("STCdR03", reuseConc="STC", reuseL1="STC", layer2Threshold=20, layer2dR = 0.03)

#newClustering("STC070", concentratorAlgo="superTriggerCellSelect", concentratorThreshold=0.7, layer1Threshold=0, layer2Threshold=20, 
#                    layer2dR = ([0] + [0.010]*7 + [0.020]*7 + [0.030]*7 + [0.040]*7 +   [0.040]*6 + [0.050]*6  +  [0.050]*12))
#newClustering("STC073",  reuseConc="STC070", layer1Threshold=3, layer2Threshold=20, 
#                    layer2dR = ([0] + [0.010]*7 + [0.020]*7 + [0.030]*7 + [0.040]*7 +   [0.040]*6 + [0.050]*6  +  [0.050]*12))
newClustering("TC", reuseConc="", layer1Threshold=0, layer2Threshold=20, layer2dR = ([0] + [0.010]*7 + [0.020]*7 + [0.030]*7 + [0.040]*7 +   [0.040]*6 + [0.050]*6  +  [0.050]*12))
newClustering("TCdR03", reuseConc="", reuseL1="TC", layer2Threshold=20, layer2dR = 0.03)

process.p = cms.Path(sum(modules[1:], modules[0]))
process.TFileService = cms.Service("TFileService", fileName = cms.string("idTupleNew.root"))

def goRandom():
    for aname in process.analyzers_().keys():
        if aname.startswith("ntuple"):
            ana = getattr(process,aname)
            if hasattr(ana, 'onlyMatched'):
                ana.onlyMatched = False
def hgcAcc(pdgId,ptmin=2):
    process.acceptance = cms.EDFilter("GenParticleSelector",
        src = cms.InputTag("genParticles"),
        cut = cms.string("abs(pdgId) == %d && pt > %g && abs(eta) > 1.6 && abs(eta) < 2.6" % (pdgId,ptmin)),
        filter = cms.bool(True),
    )
    process.p.insert(0, process.acceptance)
def xdup():
    process.options.numberOfThreads = cms.untracked.uint32(2)
    process.options.numberOfStreams = cms.untracked.uint32(0)


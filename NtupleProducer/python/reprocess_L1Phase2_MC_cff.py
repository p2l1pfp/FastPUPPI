import FWCore.ParameterSet.Config as cms

HcalHardcodeGeometryEP = cms.ESProducer("HcalHardcodeGeometryEP",
    UseOldLoader = cms.bool(False)
)


HcalTPGCoderULUT = cms.ESProducer("HcalTPGCoderULUT",
    FGLUTs = cms.FileInPath('CalibCalorimetry/HcalTPGAlgos/data/HBHE_FG_LUT.dat'),
    LUTGenerationMode = cms.bool(True),
    MaskBit = cms.int32(32768),
    RCalibFile = cms.FileInPath('CalibCalorimetry/HcalTPGAlgos/data/RecHit-TPG-calib.dat'),
    inputLUTs = cms.FileInPath('CalibCalorimetry/HcalTPGAlgos/data/inputLUTcoder_physics.dat'),
    read_Ascii_LUTs = cms.bool(False),
    read_FG_LUTs = cms.bool(False),
    read_XML_LUTs = cms.bool(False)
)


HcalTrigTowerGeometryESProducer = cms.ESProducer("HcalTrigTowerGeometryESProducer")

CaloGeometryBuilder = cms.ESProducer("CaloGeometryBuilder",
    SelectedCalos = cms.vstring('HCAL', 
        'ZDC', 
        'EcalBarrel', 
        'TOWER', 
        'HGCalEESensitive', 
        'HGCalHESiliconSensitive')
)


CaloTPGTranscoder = cms.ESProducer("CaloTPGTranscoderULUTs",
    HFTPScaleShift = cms.PSet(
        NCT = cms.int32(1),
        RCT = cms.int32(3)
    ),
    LUTfactor = cms.vint32(1, 2, 5, 0),
    RCTLSB = cms.double(0.25),
    ZS = cms.vint32(4, 2, 1, 0),
    hcalLUT1 = cms.FileInPath('CalibCalorimetry/CaloTPG/data/outputLUTtranscoder_physics.dat'),
    hcalLUT2 = cms.FileInPath('CalibCalorimetry/CaloTPG/data/TPGcalcDecompress2.txt'),
    ietaLowerBound = cms.vint32(1, 18, 27, 29),
    ietaUpperBound = cms.vint32(17, 26, 28, 32),
    nominal_gain = cms.double(0.177),
    read_Ascii_Compression_LUTs = cms.bool(False),
    read_Ascii_RCT_LUTs = cms.bool(False)
)


CaloTopologyBuilder = cms.ESProducer("CaloTopologyBuilder")


CaloTowerHardcodeGeometryEP = cms.ESProducer("CaloTowerHardcodeGeometryEP")


CaloTowerTopologyEP = cms.ESProducer("CaloTowerTopologyEP")

from L1Trigger.L1THGCal.hgcalTriggerPrimitives_cff import *
hgcl1tpg_step = cms.Sequence(hgcalTriggerPrimitives)

from SimCalorimetry.EcalEBTrigPrimProducers.ecalEBTriggerPrimitiveDigis_cff import *
EcalEBtp_step = cms.Sequence(simEcalEBTriggerPrimitiveDigis)

from SimCalorimetry.HcalTrigPrimProducers.hcalTTPDigis_cff import *
HcalTPsimulation_step = cms.Sequence(hcalTTPSequence)

from Configuration.StandardSequences.SimL1Emulator_cff import *
L1simulation_step = cms.Sequence(SimL1Emulator)

from L1Trigger.TrackTrigger.TrackTrigger_cff import *
from L1Trigger.TrackTrigger.TTStubAlgorithmRegister_cfi import *
TTClusterStub = cms.Sequence(TrackTriggerClustersStubs)

from L1Trigger.TrackFindingTracklet.L1TrackletTracks_cff import *
#TTClusterAssociatorFromPixelDigis.digiSimLinks          = cms.InputTag( "simSiPixelDigis","Tracker" )
L1TrackTrigger_step = cms.Sequence(L1TrackletTracks)

reprocess_L1Phase2_MC = cms.Sequence(
    hgcl1tpg_step + EcalEBtp_step + TTClusterStub + L1TrackTrigger_step + HcalTPsimulation_step
)

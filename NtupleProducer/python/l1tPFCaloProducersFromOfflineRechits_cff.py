import FWCore.ParameterSet.Config as cms

from FastPUPPI.NtupleProducer.l1tPFEcalProducerFromOfflineRechits_cfi import *
from FastPUPPI.NtupleProducer.l1tPFHcalProducerFromOfflineRechits_cfi import *
from FastPUPPI.NtupleProducer.l1tPFHFProducerFromOfflineRechits_cfi import *
from FastPUPPI.NtupleProducer.l1tPFHGCalProducerFromOfflineRechits_cfi import *
from FastPUPPI.NtupleProducer.l1tPFHGCalBHProducerFromOfflineRechits_cfi import *

l1tPFCaloProducersFromOfflineRechits = cms.Sequence(
    l1tPFEcalProducerFromOfflineRechits +
    l1tPFHcalProducerFromOfflineRechits +
    l1tPFHFProducerFromOfflineRechits +
    l1tPFHGCalEEProducerFromOfflineRechits +
    l1tPFHGCalFHProducerFromOfflineRechits +
    l1tPFHGCalBHProducerFromOfflineRechits 
)


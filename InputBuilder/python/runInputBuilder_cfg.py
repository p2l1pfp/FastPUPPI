import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:/afs/cern.ch/work/n/ntran/private/Correlator/eos/cms/store/mc/TTI2023Upg14D/Zmumu_TuneZ2star_14TeV_Eta4-pythia6/GEN-SIM-DIGI-RAW/PU140bx25_PH2_1K_FB_V3-v2/00000/08B6D2DB-D5E6-E311-9356-002618943975.root'
    )
)

process.demo = cms.EDAnalyzer('InputBuilder',
	L1TrackTag  = cms.InputTag('TTTracksFromPixelDigis','Level1TTTracks'),
	EcalTPTag   = cms.InputTag('simEcalTriggerPrimitiveDigis'),
	HcalTPTag   = cms.InputTag('simHcalTriggerPrimitiveDigis'),
)


process.p = cms.Path(process.demo)

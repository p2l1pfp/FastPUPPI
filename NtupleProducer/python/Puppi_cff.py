import FWCore.ParameterSet.Config as cms

puppiCentral = cms.VPSet(
                 cms.PSet(
                  algoId           = cms.int32(5),  #0 is default Puppi
                  useCharged       = cms.bool(True),
                  applyLowPUCorr   = cms.bool(False),
                  combOpt          = cms.int32(0),
                  cone             = cms.double(0.4),
                  rmsPtMin         = cms.double(0.5),
                  rmsScaleFactor   = cms.double(1.0)
                 )
                )

puppiForward = cms.VPSet(
                cms.PSet(
                 algoId         = cms.int32(5),  #0 is default Puppi
                 useCharged     = cms.bool(False),
                 applyLowPUCorr = cms.bool(False),
                 combOpt        = cms.int32(0),
                 cone           = cms.double(0.4),
                 rmsPtMin       = cms.double(0.5),
                 rmsScaleFactor = cms.double(1.0)
                 )
                )

puppi = cms.PSet(#"PuppiProducer",
                       puppiDiagnostics = cms.bool(False),
                       puppiForLeptons = cms.bool(False),
                       UseDeltaZCut   = cms.bool(True),
                       DeltaZCut      = cms.double(0.3),
                       candName       = cms.InputTag('particleFlow'),
                       vertexName     = cms.InputTag('offlinePrimaryVertices'),
                       #candName      = cms.string('packedPFCandidates'),
                       #vertexName     = cms.string('offlineSlimmedPrimaryVertices'),
                       applyCHS       = cms.bool  (False),
                       invertPuppi    = cms.bool  (False),
                       useExp         = cms.bool  (False),
                       MinPuppiWeight = cms.double(0.01),
                       useExistingWeights = cms.bool(False),
                       useWeightsNoLep    = cms.bool(False),
                       clonePackedCands   = cms.bool(False), # should only be set to True for MiniAOD
                       vtxNdofCut     = cms.int32(4),
                       vtxZCut        = cms.double(24),
                       algos          = cms.VPSet( 
                        cms.PSet( 
                         etaMin = cms.vdouble(0.),
                         etaMax = cms.vdouble(2.5),
                         ptMin  = cms.vdouble(0.),
                         MinNeutralPt   = cms.vdouble(4.0),
                         MinNeutralPtSlope   = cms.vdouble(0.0),
                         RMSEtaSF = cms.vdouble(1.0),
                         MedEtaSF = cms.vdouble(1.0),
                         EtaMaxExtrap = cms.double(2.0),
                         puppiAlgos = puppiCentral
                        ),
                        cms.PSet( 
                         etaMin              = cms.vdouble( 2.5),
                         etaMax              = cms.vdouble( 3.0),
                         ptMin               = cms.vdouble( 0.),
                         MinNeutralPt        = cms.vdouble( 4.0),
                         MinNeutralPtSlope   = cms.vdouble(0.00),
                         RMSEtaSF            = cms.vdouble(1.20),
                         MedEtaSF            = cms.vdouble(0.90),
                         EtaMaxExtrap        = cms.double( 2.0),
                         puppiAlgos = puppiForward
                        ),
                        cms.PSet( 
                         etaMin              = cms.vdouble( 3.0),
                         etaMax              = cms.vdouble( 10.0),
                         ptMin               = cms.vdouble( 0.),
                         MinNeutralPt        = cms.vdouble( 4.0),
                         MinNeutralPtSlope   = cms.vdouble(0.00),
                         RMSEtaSF            = cms.vdouble(0.95),
                         MedEtaSF            = cms.vdouble(0.75),
                         EtaMaxExtrap        = cms.double (2.0),
                         puppiAlgos = puppiForward
                        ),
                      )
)

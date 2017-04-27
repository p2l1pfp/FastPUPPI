import FWCore.ParameterSet.Config as cms

genMuonsFromGun = cms.EDFilter("CandViewSelector",
    src = cms.InputTag("genParticles"),
    cut = cms.string("abs(pdgId) == 13 && isPromptFinalState && numberOfMothers == 0"),
)

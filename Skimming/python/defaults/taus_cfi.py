import FWCore.ParameterSet.Config as cms

defaultSelector = cms.EDFilter("PATTauSelector", 
	   filter = cms.bool(True),
           src = cms.InputTag("cleanLayer1Taus"),
           cut = cms.string('abs( eta ) <= 2.5 & pt >= 7')#GeV
)

patMatchedTauSelector = cms.EDFilter("PATTauMatchedSelector", 
   src = cms.InputTag("cleanLayer1Taus"),
   pdgId = cms.int32(15),
   status = cms.uint32(0),
   autoCharge = cms.bool(True)
)

patJetMatchedTauSelector = cms.EDFilter("PATTauJetMatchedSelector", 
   src = cms.InputTag("cleanLayer1Taus"),
)

from SuSyAachen.Skimming.jetSelectors_cfi import candViewCountFilter 
tauCountSelector = candViewCountFilter.clone(
           filter = cms.bool(True),
           src = cms.InputTag("cleanLayer1Taus"),
           minNumber = 1
)

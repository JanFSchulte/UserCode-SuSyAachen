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

from SuSyAachen.Skimming.tauSelection_cfi import scIsolatedTaus
tauScIsolationSelector = scIsolatedTaus.clone(
    filter = cms.bool(True),
     src = cms.InputTag("selectedPatTausPF"),
    otherTauSource = cms.InputTag("selectedPatTausPFShrinkingCone"),
    otherTauId = cms.string("byIsolation"),
    dRMin = cms.double(0.15)
)


tauMuonCleaner = cms.EDProducer('tauMuonCleaner',
                                     filter = cms.bool(True),
                                     src = cms.InputTag("selectedPatTausPF"),
                                     leptSrc = cms.InputTag("basicMuons"),
                                     dRJetLepton = cms.double(0.15)
                                     )

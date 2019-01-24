import FWCore.ParameterSet.Config as cms

defaultSelector = cms.EDFilter("PATJetSelector", 
       filter = cms.bool(True),
           src = cms.InputTag("cleanLayer1Jets"),
           cut = cms.string('abs( eta ) <= 2.0 & pt >= 10')#GeV
)

from SuSyAachen.Skimming.jetSelectors_cfi import patPFJetIDFilter
patPFJetIDSelector =  patPFJetIDFilter.clone(
           filter = cms.bool(True),
           src = cms.InputTag("cleanLayer1JetsAK4"),
           version = cms.string('FIRSTDATA'),
           quality = cms.string('LOOSE')
)


from SuSyAachen.Skimming.jetSelectors_cfi import jetMuonCleaner
muonCleanJets = jetMuonCleaner.clone(
    src = cms.InputTag("selectedPatJetsPF")
)

from SuSyAachen.Skimming.jetSelectors_cfi import jetElectronCleaner
electronCleanJets = jetElectronCleaner.clone(
    src = cms.InputTag("selectedPatJetsPF")
)




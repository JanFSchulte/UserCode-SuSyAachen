import FWCore.ParameterSet.Config as cms

defaultSelector = cms.EDFilter("PATJetSelector", 
       filter = cms.bool(True),
           src = cms.InputTag("cleanLayer1Jets"),
           cut = cms.string('abs( eta ) <= 2.0 & pt >= 10')#GeV
)

from SuSyAachen.Skimming.jetSelectors_cfi import patJetCountFilter as patJetCountFilterOrig
patJetCountFilter = patJetCountFilterOrig.clone(
           filter = cms.bool(True),
           src = cms.InputTag("cleanLayer1JetsAK4"),
           minNumber = 1
)

genJetSelector = cms.EDFilter("GenJetSelector", 
           filter = cms.bool(True),
           src = cms.InputTag("ak4GenJets"),
           cut = cms.string("pt > 100")
)

from SuSyAachen.Skimming.jetSelectors_cfi import patJetFlagFilter
patJetFlagSelector =  patJetFlagFilter.clone(
           filter = cms.bool(True),
           src = cms.InputTag("cleanLayer1JetsAK4"),
           cut = cms.string("pat::Flags::Overlap::Electrons")
)

from SuSyAachen.Skimming.jetSelectors_cfi import patPFJetIDFilter
patPFJetIDSelector =  patPFJetIDFilter.clone(
           filter = cms.bool(True),
           src = cms.InputTag("cleanLayer1JetsAK4"),
           version = cms.string('FIRSTDATA'),
           quality = cms.string('LOOSE')
)

from SuSyAachen.Skimming.jetSelectors_cfi import candViewCountFilter 
genJetCountSelector = candViewCountFilter.clone(
           filter = cms.bool(True),
           src = cms.InputTag("ak4GenJets"),
           minNumber = 1
)

from SuSyAachen.Skimming.jetSelectors_cfi import jetMuonCleaner
muonCleanJets = jetMuonCleaner.clone(
    src = cms.InputTag("selectedPatJetsPF")
)

from SuSyAachen.Skimming.jetSelectors_cfi import jetElectronCleaner
electronCleanJets = jetElectronCleaner.clone(
    src = cms.InputTag("selectedPatJetsPF")
)


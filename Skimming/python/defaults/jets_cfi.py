import FWCore.ParameterSet.Config as cms

defaultSelector = cms.EDFilter("PATJetSelector", 
	   filter = cms.bool(True),
           src = cms.InputTag("cleanLayer1Jets"),
           cut = cms.string('abs( eta ) <= 2.0 & pt >= 10')#GeV
)

from SuSyAachen.Skimming.jetSelectors_cfi import patJetCountFilter as patJetCountFilterOrig
patJetCountFilter = patJetCountFilterOrig.clone(
           filter = cms.bool(True),
           src = cms.InputTag("cleanLayer1JetsAK5"),
           minNumber = 1
)

genJetSelector = cms.EDFilter("GenJetSelector", 
           filter = cms.bool(True),
           src = cms.InputTag("ak5GenJets"),
           cut = cms.string("pt > 100")
)

from SuSyAachen.Skimming.jetSelectors_cfi import patJetFlagFilter
patJetFlagSelector =  patJetFlagFilter.clone(
           filter = cms.bool(True),
           src = cms.InputTag("cleanLayer1JetsAK5"),
           cut = cms.string("pat::Flags::Overlap::Electrons")
)

from SuSyAachen.Skimming.jetSelectors_cfi import patPFJetIDFilter
patPFJetIDSelector =  patPFJetIDFilter.clone(
           filter = cms.bool(True),
           src = cms.InputTag("cleanLayer1JetsAK4"),
           version = cms.string('FIRSTDATA'),
           quality = cms.string('LOOSE')
)

from SuSyAachen.Skimming.jetSelectors_cfi import patJetIDFilter
patJetIDSelector =  patJetIDFilter.clone(
           filter = cms.bool(True),
           src = cms.InputTag("cleanLayer1JetsAK4"),
           version = cms.string('PURE09'),
           quality = cms.string('LOOSE')
)


from SuSyAachen.Skimming.jetSelectors_cfi import candViewCountFilter 
genJetCountSelector = candViewCountFilter.clone(
           filter = cms.bool(True),
           src = cms.InputTag("ak4GenJets"),
           minNumber = 1
)

from SuSyAachen.Skimming.htFilter_cfi import htFilter
htJetFilter = htFilter.clone(
            src = cms.InputTag("selectedPatJetsPF"),
            minHT = cms.double(100.00)
)

from SuSyAachen.Skimming.htFilter_cfi import metSqrtHtFilter
metSqrtHtJetFilter = metSqrtHtFilter.clone(
            src = cms.InputTag("selectedPatJetsPF"),
            minCut = cms.double(13.5)
)

from SuSyAachen.Skimming.htFilter_cfi import rawHtFilter
rawHtJetFilter = rawHtFilter.clone(
            src = cms.InputTag("selectedPatJetsPF"),
            minHT = cms.double(100.00)
)

from SuSyAachen.Skimming.mhtFilter_cfi import mhtFilter
mhtJetFilter = mhtFilter.clone(
            src = cms.InputTag("selectedPatJetsPF"),
            minMHT = cms.double(100.00)
)

from SuSyAachen.Skimming.jetSelectors_cfi import jetMuonCleaner
muonCleanJets = jetMuonCleaner.clone(
    src = cms.InputTag("selectedPatJetsPF")
)

from SuSyAachen.Skimming.jetSelectors_cfi import jetElectronCleaner
electronCleanJets = jetElectronCleaner.clone(
    src = cms.InputTag("selectedPatJetsPF")
)

from SuSyAachen.Skimming.jetSelectors_cfi import jetTauCleaner
tauCleanJets = jetTauCleaner.clone(
    src = cms.InputTag("selectedPatJetsPF")
)

from SuSyAachen.Skimming.jetSelectors_cfi import resCorrectedJetProducer
resCorrectedJets = resCorrectedJetProducer.clone(
    src = cms.InputTag("selectedPatJetsPF")
)

from SuSyAachen.Skimming.jetSelectors_cfi import unCorrectedJetProducer
unCorrectedJets = unCorrectedJetProducer.clone(
    src = cms.InputTag("selectedPatJetsPF")
)

from SuSyAachen.Skimming.jetSelectors_cfi import fastJetUnCorrectedJetProducer
fastJetUnCorrectedJets = fastJetUnCorrectedJetProducer.clone(
    src = cms.InputTag("selectedPatJetsPF")
)


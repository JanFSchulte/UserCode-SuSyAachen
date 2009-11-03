import FWCore.ParameterSet.Config as cms

filterJets = cms.bool(False)

basicJets = cms.EDProducer("PATJetSelector", filter = filterJets,
  src = cms.InputTag("cleanLayer1JetsAK5"),
  cut = cms.string("pt > 20 & abs( eta ) < 2.5")
)

METFilter =  cms.EDProducer("PATMETSelector", filter = cms.bool(True),
  src = cms.InputTag("layer1METsAK5"),
  cut = cms.string("pt > 100")
)

from SuSyAachen.Skimming.jetSelectors_cfi import patJetCountFilter
jetFilterHigh = patJetCountFilter.clone( filter = filterJets,
    src = "basicJets",
    minNumber = 1,
    cut = "pt > 100"
)

jetFilterLow = patJetCountFilter.clone( filter = filterJets,
    src = "basicJets",
    minNumber = 3,
    cut = "pt > 50"
)

seqSelectJetMET = cms.Sequence( METFilter + basicJets * ( jetFilterHigh + jetFilterLow) )

import FWCore.ParameterSet.Config as cms

filterJets = cms.bool(False)
filterMET = cms.bool(True)

basicJets = cms.EDFilter("PATJetSelector", filter = filterJets,
  src = cms.InputTag("cleanLayer1JetsAK5"),
  cut = cms.string("pt > 50 & abs( eta ) < 2.5")
)

METFilter =  cms.EDFilter("PATMETSelector", filter = filterMET,
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

seqSelectJetMET = cms.Sequence( METFilter  +basicJets * ( cms.ignore(jetFilterHigh) + cms.ignore(jetFilterLow)) )

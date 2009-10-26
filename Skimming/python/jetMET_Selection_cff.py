import FWCore.ParameterSet.Config as cms

cleanMuons = cms.EDProducer("CandSelector",
                       src = cms.InputTag("cleanLayer1Muons"),
                       cut = cms.string("pt > 10 & abs( eta ) < 2.5")
                       )

seqSelectJetMET = cms.Sequence(cleanMuons)

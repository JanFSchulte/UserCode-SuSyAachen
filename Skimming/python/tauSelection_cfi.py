import FWCore.ParameterSet.Config as cms
scIsolatedTaus = cms.EDFilter("PATTauToTauMatchSelector", 
                                  src = cms.InputTag("selectedPatTausPF"),
                                  otherTauSource = cms.InputTag("selectedPatTausPFShrinkingCone"),
                                  otherTauId = cms.string("byIsolation"),
                                  dRMin = cms.double(0.15)
                                
                                )

vicinityTaus = cms.EDFilter("PATTauVicinitySelector",
                              src = cms.InputTag("selectedPatTausPF"),
                              pfCandSource = cms.InputTag("particleFlow","","RECO"),
                              dRMax = cms.double(0.5)
                              )

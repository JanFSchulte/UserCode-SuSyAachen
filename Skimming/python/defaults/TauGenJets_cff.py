import FWCore.ParameterSet.Config as cms
def TauGenJets(process):
    from PhysicsTools.JetMCAlgos.TauGenJets_cfi import tauGenJets
    process.tauGenJets = tauGenJets
    process.hadronicGenTauJets = cms.EDFilter("TauGenJetDecayModeSelector",
      src = cms.InputTag("tauGenJets"),
      select = cms.vstring('oneProng0Pi0', 'oneProng1Pi0', 'oneProng2Pi0', 'oneProngOther',
                           'threeProng0Pi0', 'threeProng1Pi0', 'threeProngOther', 'rare'),
      filter = cms.bool(False)
    )
    process.seqTauGenJets = cms.Sequence(process.tauGenJets*process.hadronicGenTauJets)

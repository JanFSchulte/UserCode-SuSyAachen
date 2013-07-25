import FWCore.ParameterSet.Config as cms

process = cms.Process("rho2011Producer")


from RecoJets.JetProducers.kt4PFJets_cfi import *
kt6PFJetsForIsolation = kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
kt6PFJetsForIsolation.Rho_EtaMax = cms.double(2.5)

process.rho2011Sequence = cms.Sequence(kt6PFJetsForIsolation)
process.p = cms.Path(process.rho2011Sequence)


process.outpath = cms.EndPath(process.out)
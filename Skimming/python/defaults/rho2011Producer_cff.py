import FWCore.ParameterSet.Config as cms

def rho2011Producer(process):


	process.load("RecoJets.JetProducers.kt4PFJets_cfi")
	process.kt6PFJetsForIsolation = process.kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
	process.kt6PFJetsForIsolation.Rho_EtaMax = cms.double(2.5)
	
	process.seqrho2011Producer = cms.Sequence(process.kt6PFJetsForIsolation)
	process.rho2011Path = cms.Path(process.seqrho2011Producer)
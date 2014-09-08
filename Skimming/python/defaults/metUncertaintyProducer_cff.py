import FWCore.ParameterSet.Config as cms

def metUncertaintyProducer(process,taskname):


	from PhysicsTools.PatUtils.tools.metUncertaintyTools import runMEtUncertainties
	runMEtUncertainties(process,jetCollection="selectedPatJetsAK5PF",electronCollection="IsoElectrons",muonCollection="IsoMuons",dRjetCleaning=0.4,tauCollection=None,addToPatDefaultSequence=False)
	process.metUncertaintyPath = cms.Path(process.metUncertaintySequence)
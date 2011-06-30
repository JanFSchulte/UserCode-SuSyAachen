import FWCore.ParameterSet.Config as cms
def RedoTaus(process):
	from copy import deepcopy
	process.load("Configuration.StandardSequences.Geometry_cff")
#	process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#	from Configuration.PyReleaseValidation.autoCond import autoCond
#	process.GlobalTag.globaltag = cms.string( autoCond[ 'mc' ] )
	process.load("Configuration.StandardSequences.MagneticField_cff")
	process.load("RecoJets.Configuration.RecoPFJets_cff")
	process.load('RecoTauTag.Configuration.RecoPFTauTag_cff')
	process.load("PhysicsTools.PatAlgos.producersLayer1.tauProducer_cff")

	process.patDefaultSequence = deepcopy(process.makePatTaus)
#	from PhysicsTools.PatAlgos.tools.tauTools import switchToPFTauHPS
#	switchToPFTauHPS(process)
	
#	process.dump = cms.EDAnalyzer('EventContentAnalyzer') 
#	isoConeFormula = '0.3'
#	process.shrinkingConePFTauProducer.builders[0].isoConeChargedHadrons = isoConeFormula
#	process.shrinkingConePFTauProducer.builders[0].isoConePiZeros = isoConeFormula
#	process.shrinkingConePFTauProducer.builders[0].isoConeNeutralHadrons = isoConeFormula
#	process.patTaus.userIsolation.pfAllParticles.deltaR = eval(isoConeFormula)
#	process.patTaus.userIsolation.pfNeutralHadron.deltaR = eval(isoConeFormula)
#	process.patTaus.userIsolation.pfChargedHadron.deltaR = eval(isoConeFormula)
#	process.patTaus.userIsolation.pfGamma.deltaR = eval(isoConeFormula)

		
	process.seqRedoTaus = cms.Sequence(process.ak5PFJets * 
					   process.PFTau *
					   process.patDefaultSequence
#					   process.makePatTaus #* process.dump
					   )

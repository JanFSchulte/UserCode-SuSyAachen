import FWCore.ParameterSet.Config as cms

def electronMVAProducer(process):


	process.load("RecoEgamma.ElectronIdentification.ElectronMVAValueMapProducer_cfi")

	
	process.seqelectronMVAProducer = cms.Sequence(process.electronMVAValueMapProducer)
	process.electronMVAPath = cms.Path(process.seqelectronMVAProducer)
